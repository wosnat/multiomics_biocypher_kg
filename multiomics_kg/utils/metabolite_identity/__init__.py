"""Metabolite identity pipeline for prepare-data step 6 (spec 2026-10-13)."""
from multiomics_kg.utils.metabolite_identity.names import (  # noqa: F401
    PERIODIC_TABLE, is_formula_like, is_generic_formula, is_polymer_formula,
    is_stereo_qualified, name_key, normalize_hmdb, parse_elements, strip_anomer_tokens)
from multiomics_kg.utils.metabolite_identity.canonical import (  # noqa: F401
    MERGE_BASES, IdentityCollisionError, IdentityOverrides, IdentityResult, MergeLink, canonicalise, collect_merge_links,
    entity_terms, load_identity_overrides, rewrite_to_canonical)
from multiomics_kg.utils.metabolite_identity.structure import (  # noqa: F401
    STRUCTURE_KINDS, apply_structure, classify_structure)
from multiomics_kg.utils.metabolite_identity.crossrefs import (  # noqa: F401
    build_identity_fields, is_isotopologue_name, pick_chebi_id, shared_chebi_terms)
from multiomics_kg.utils.metabolite_identity.stereo import (  # noqa: F401
    finalize_stereo, is_anomer_pair, stereo_families, stereo_links, unspecified_members)
from multiomics_kg.utils.metabolite_identity.paper import (  # noqa: F401
    STEREO_ASSIGNMENTS, Decision, Override, PaperRow, PaperTable, ResolutionResult, StereoContext,
    collapse_anomers, dominant_form, harvest_paper_tables, load_paper_overrides, resolve_paper_rows)
from multiomics_kg.utils.metabolite_identity.roles import (  # noqa: F401
    build_role_hierarchy, build_subtree_roles, partition_role_strings, role_node_id)
from multiomics_kg.utils.metabolite_identity.classes import (  # noqa: F401
    COMPOUND_BRITE_TREES, LEVEL_KINDS, download_compound_brite, parse_compound_brite,
    prune_compound_brite)
from multiomics_kg.utils.metabolite_identity.lumping import (  # noqa: F401
    LUMPING_MIN_LEVEL, LUMPING_MIN_METABOLITES, SEED_LUMPING_IDS, substrate_breadth,
    threshold_lumping, update_lumping_pin)
from multiomics_kg.utils.metabolite_identity.outputs import (  # noqa: F401
    build_identity_report, build_mapping_v2, remap_rows, repoint_reason, validate_build,
    write_all_atomic, write_resolved_tables)
