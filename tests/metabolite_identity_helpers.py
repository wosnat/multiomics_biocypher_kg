"""Shared builders for metabolite-identity unit tests (no network, no MNX)."""
from multiomics_kg.utils.chebi_utils import ChebiIndex


def term(name, *, charge=None, alt=(), kegg=(), conj=(), taut=(), enan=(), is_a=(),
         gformula=None, inchikey=None):
    return {"name": name, "charge": charge, "alt_ids": list(alt), "kegg": list(kegg),
            "conjugates": list(conj), "tautomers": list(taut), "enantiomers": list(enan),
            "is_a": list(is_a), "generalized_formula": gformula, "inchikey": inchikey}


def make_chebi(**terms) -> ChebiIndex:
    """make_chebi(CHEBI_17632=term(...)) — underscores become colons."""
    return ChebiIndex({k.replace("_", ":"): v for k, v in terms.items()}, release="256")


def cpd(name, *, formula=None, mass=None, inchikey=None, smiles=None, chebi_id=None,
        mnxm_id=None, hmdb_id=None, ev=("metabolism",), pathways=()):
    return {"name": name, "formula": formula, "mass": mass, "inchikey": inchikey,
            "smiles": smiles, "chebi_id": chebi_id, "mnxm_id": mnxm_id, "hmdb_id": hmdb_id,
            "evidence_sources": list(ev), "pathways": list(pathways)}
