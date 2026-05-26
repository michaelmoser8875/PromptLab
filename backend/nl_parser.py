"""Phase 1 - Natural language to SMILES identifier.

Turns a free-text reaction prompt into structured molecular data:

    {
        "reactants": [{"name": ..., "smiles": ..., "source": ...}, ...],
        "products":  [...],
        "context":   "..."
    }

PubChem is queried first for each compound name; if that misses, Claude is
asked for a SMILES as a fallback. Every SMILES is validated and canonicalized
with RDKit before it leaves this module.
"""
from __future__ import annotations

import json
import sys
import urllib.parse
from dataclasses import asdict, dataclass
from typing import Optional

import anthropic
import requests
from rdkit import Chem
from rdkit import RDLogger

# RDKit prints noisy stderr warnings on invalid SMILES; we treat those as
# resolution failures rather than user-visible errors.
RDLogger.DisableLog("rdApp.*")


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_TIMEOUT = 10.0
MODEL = "claude-opus-4-7"


@dataclass
class Compound:
    name: str
    smiles: str
    source: str  # "pubchem" | "llm"


@dataclass
class ParsedReaction:
    reactants: list[Compound]
    products: list[Compound]
    context: str

    def to_dict(self) -> dict:
        return {
            "reactants": [asdict(c) for c in self.reactants],
            "products": [asdict(c) for c in self.products],
            "context": self.context,
        }


class ResolutionError(Exception):
    """Raised when a compound name cannot be resolved to a valid SMILES."""


# ---------------------------------------------------------------------------
# Compound name extraction (LLM)
# ---------------------------------------------------------------------------

_EXTRACTION_SYSTEM = """\
You extract chemical reaction components from natural-language prompts.

Given a user prompt that describes (or asks about) a chemical reaction,
identify three things:

- reactants: the compounds being mixed or reacted. Always extract these when
  the prompt names any chemicals.
- products: only compounds the user EXPLICITLY names as products. Usually
  empty - leave it as an empty list if the user only mentions the inputs.
- context: a short note on what kind of reaction is implied. Examples:
  "acid-base neutralization", "esterification", "SN2 substitution",
  "combustion in air", "precipitation". Empty string if no specific reaction
  type is implied.

For each compound, return the canonical English name in lowercase
(e.g. "methanol", "acetic acid", "sodium hydroxide"). Prefer common names
over IUPAC where both exist - they resolve more reliably against external
databases.

If the prompt does not describe a chemical reaction at all, return empty
reactants and products lists with an empty context.
"""

_EXTRACTION_SCHEMA = {
    "type": "object",
    "properties": {
        "reactants": {"type": "array", "items": {"type": "string"}},
        "products": {"type": "array", "items": {"type": "string"}},
        "context": {"type": "string"},
    },
    "required": ["reactants", "products", "context"],
    "additionalProperties": False,
}


def extract_compounds(prompt: str, client: anthropic.Anthropic) -> dict:
    """Return {'reactants': [str], 'products': [str], 'context': str}."""
    response = client.messages.create(
        model=MODEL,
        max_tokens=2000,
        system=_EXTRACTION_SYSTEM,
        messages=[{"role": "user", "content": prompt}],
        output_config={
            "format": {
                "type": "json_schema",
                "schema": _EXTRACTION_SCHEMA,
            },
        },
    )
    text = next(b.text for b in response.content if b.type == "text")
    return json.loads(text)


# ---------------------------------------------------------------------------
# Name to SMILES resolution
# ---------------------------------------------------------------------------

# PubChem renamed CanonicalSMILES to SMILES in 2024; try both so this works
# against both old and new deployments.
_PUBCHEM_PROPS = ("SMILES", "CanonicalSMILES", "IsomericSMILES")


def _pubchem_lookup(name: str) -> Optional[str]:
    """Return a SMILES string from PubChem, or None on miss/failure."""
    encoded = urllib.parse.quote(name, safe="")
    for prop in _PUBCHEM_PROPS:
        url = f"{PUBCHEM_BASE}/compound/name/{encoded}/property/{prop}/JSON"
        try:
            r = requests.get(url, timeout=PUBCHEM_TIMEOUT)
        except requests.RequestException:
            continue
        if r.status_code != 200:
            continue
        try:
            props = r.json()["PropertyTable"]["Properties"]
        except (KeyError, ValueError):
            continue
        if not props:
            continue
        smiles = props[0].get(prop)
        if smiles:
            return smiles
    return None


_LLM_SMILES_SYSTEM = (
    "You return the canonical SMILES string for a named chemical compound. "
    "Reply with the SMILES only - no explanation, no surrounding quotes, no "
    "code fences, no prose. If you do not know the structure with confidence, "
    "reply with exactly: UNKNOWN"
)


def _llm_smiles(name: str, client: anthropic.Anthropic) -> Optional[str]:
    """Last-resort fallback: ask Claude for the SMILES of a named compound."""
    response = client.messages.create(
        model=MODEL,
        max_tokens=200,
        system=_LLM_SMILES_SYSTEM,
        messages=[{"role": "user", "content": name}],
    )
    text = "".join(b.text for b in response.content if b.type == "text").strip()
    if not text or text.upper() == "UNKNOWN":
        return None
    return text


def _canonicalize(smiles: str) -> Optional[str]:
    """Validate via RDKit and return the canonical form, or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def resolve(name: str, client: anthropic.Anthropic) -> Compound:
    """Resolve a compound name to a Compound with a validated SMILES.

    PubChem is tried first; the LLM is a fallback. RDKit canonicalizes both.
    Raises ResolutionError if neither path produces valid SMILES.
    """
    raw = _pubchem_lookup(name)
    if raw:
        canonical = _canonicalize(raw)
        if canonical:
            return Compound(name=name, smiles=canonical, source="pubchem")

    raw = _llm_smiles(name, client)
    if raw:
        canonical = _canonicalize(raw)
        if canonical:
            return Compound(name=name, smiles=canonical, source="llm")

    raise ResolutionError(
        f"Could not resolve {name!r} to a valid SMILES via PubChem or LLM."
    )


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def parse_reaction(
    prompt: str,
    client: Optional[anthropic.Anthropic] = None,
) -> ParsedReaction:
    """Parse a natural-language reaction prompt into structured molecular data.

    Raises ResolutionError if any extracted compound cannot be resolved.
    """
    if client is None:
        client = anthropic.Anthropic()

    extracted = extract_compounds(prompt, client)
    reactants = [resolve(n, client) for n in extracted["reactants"]]
    products = [resolve(n, client) for n in extracted["products"]]
    return ParsedReaction(
        reactants=reactants,
        products=products,
        context=extracted.get("context", ""),
    )


def _main() -> None:
    prompt = " ".join(sys.argv[1:]).strip()
    if not prompt:
        prompt = "what happens when methanol reacts with acetic acid?"
        print(f"(no prompt given; using default: {prompt!r})", file=sys.stderr)
    result = parse_reaction(prompt)
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    _main()
