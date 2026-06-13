#!/usr/bin/env python3
"""
Convert a CyanoRak GenBank export (.gbff) into the (.gff + .gbk) pair that
``build_gene_mapping`` consumes, bridging CyanoRak's own locus tags to the
canonical NCBI locus tags via a protein-sequence match table.

Why this exists
---------------
A handful of CyanoRak strains (e.g. Prochlorococcus MIT1327) are **not** on the
public CyanoRak server export, so ``prepare_data`` step 0.2 cannot download the
usual ``Pro_<strain>.gff`` / ``.gbk`` pair. Instead the CyanoRak team provides a
single annotated ``.gbff``. That file carries every field the KG wants
(CK cluster numbers, CyanoRak/TIGR roles, COG/eggNOG, Pfam, GO, EC) in its CDS
``/note`` qualifiers — but its locus tags (``ProMIT1327_NNNNN``) do not match
the NCBI assembly's (``PMIT1327_NNNNN`` / ``PMIT1327_RSNNNNN``), and the
numbering is not gene-for-gene aligned (NCBI PGAP re-annotated independently,
so number-based alignment mismatches ~20% of genes by protein).

This script:
  1. Reads a bridge CSV (``<cyanorak_locus_tag>, locus_tag``) produced by
     ``scripts/map_img_to_ncbi_proteins.py`` (diamond protein match).
  2. Emits ``<organism>.gff`` in CyanoRak GFF format — one CDS row per mapped
     gene, with the annotation columns ``build_gene_mapping`` keeps
     (``_CYANORAK_COLS``), parsed out of the gbff notes. ``ID`` = the CyanoRak
     ORF Id (the join key the gbk maps). Code/id list fields are
     comma-separated (matching the real CyanoRak export, e.g. ``D.1.2,J``);
     description fields keep their literal spaces and percent-encode structural
     characters.
  3. Emits ``<organism>.gbk`` = the gbff filtered to mapped CDS, with each
     ``/locus_tag`` rewritten to the bridged **NCBI** locus tag. So
     ``_get_cyanorak_id_map_from_gbk`` yields ``ORF Id -> NCBI locus_tag`` and
     the standard locus_tag merge in ``load_gff_from_ncbi_and_cyanorak``
     succeeds (no unsafe coordinate fallback — multi-contig drafts repeat
     coordinates across contigs).

Genes with no bridge mapping are dropped from BOTH outputs so the merge never
manufactures phantom CyanoRak-only genes. When two CyanoRak ORFs bridge to the
same NCBI locus tag, the longer protein wins (deterministic).

Reproduce (MIT1327)
-------------------
    # 1. extract gbff translations keyed by the CyanoRak locus tag
    #    (>ProMIT1327_NNNNN) into cyanorak/cyanorak_proteins.faa
    # 2. diamond bridge -> cyanorak/cyanorak_to_ncbi_locus_tag.csv:
    uv run python scripts/map_img_to_ncbi_proteins.py \\
        --img-faa  cache/data/Prochlorococcus/genomes/MIT1327/cyanorak/cyanorak_proteins.faa \\
        --ncbi-faa cache/data/Prochlorococcus/genomes/MIT1327/protein.faa \\
        --gene-mapping cache/data/Prochlorococcus/genomes/MIT1327/gene_mapping.csv \\
        --output   cache/data/Prochlorococcus/genomes/MIT1327/cyanorak/cyanorak_to_ncbi_locus_tag.csv \\
        --source-id-col cyanorak_locus_tag
    # 3. this converter:
    uv run python -m multiomics_kg.download.convert_cyanorak_gbff_to_gff \\
        --gbff "data/Prochlorococcus/papers_and_supp/Soussan 2025/Pro_MIT1327_LVHR00000000.gbff" \\
        --bridge cache/data/Prochlorococcus/genomes/MIT1327/cyanorak/cyanorak_to_ncbi_locus_tag.csv \\
        --gene-mapping cache/data/Prochlorococcus/genomes/MIT1327/gene_mapping.csv \\
        --organism Pro_MIT1327 \\
        --out-dir  cache/data/Prochlorococcus/genomes/MIT1327/cyanorak
    # 4. rebuild gene_mapping + annotations (cyanorak now merges by locus_tag):
    bash scripts/prepare_data.sh --strains MIT1327 --steps 0 1 2 3 --skip-cyanorak --force
"""

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO

GFF_SOURCE = "cyanorak"

# CyanoRak GFF attribute order written into column 9 (subset of _CYANORAK_COLS
# that we can source from the gbff notes).
GFF_COLS = [
    "ID", "Name", "product", "cluster_number",
    "cyanorak_Role", "cyanorak_Role_description",
    "tIGR_Role", "tIGR_Role_description",
    "Ontology_term", "ontology_term_description",
    "eggNOG", "eggNOG_description", "kegg", "kegg_description",
    "protein_domains",
]


def _enc(value: str) -> str:
    """GFF3-encode a value the way CyanoRak exports do: structural chars
    percent-encoded, spaces kept literal (matches Pro_MED4.gff)."""
    if value is None:
        return ""
    out = value.replace("%", "%25")
    out = out.replace(";", "%3B").replace("=", "%3D")
    out = out.replace(",", "%2C").replace("&", "%26")
    out = out.replace("\t", "%09").replace("\n", " ").replace("\r", " ")
    return out


def _notes(feat) -> list[str]:
    return list(feat.qualifiers.get("note", []))


def _note_values(notes: list[str], prefix: str) -> list[str]:
    """All note values whose text starts with `prefix` (prefix stripped)."""
    return [n[len(prefix):].strip() for n in notes if n.startswith(prefix)]


def _strip_code(desc: str) -> str:
    """'157=Unknown function / General' -> 'Unknown function / General'.

    CyanoRak role/EC/GO def. notes are 'CODE=description'; the GFF description
    columns hold only the description (see Pro_MED4.gff)."""
    return desc.split("=", 1)[1].strip() if "=" in desc else desc.strip()


def _code_list(notes: list[str], prefix: str) -> list[str]:
    """Code/id tokens for a note prefix: a gene may carry several notes and/or
    a single space-separated note; CyanoRak GFF joins them with commas."""
    out: list[str] = []
    for v in _note_values(notes, prefix):
        out.extend(v.split())
    return out


def load_bridge(path: Path, source_col: str) -> dict[str, str]:
    """Load <cyanorak_locus_tag> -> NCBI locus_tag from the diamond bridge CSV."""
    bridge: dict[str, str] = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames or []
        src = source_col if source_col in cols else cols[0]
        for row in reader:
            cyan = (row.get(src) or "").strip()
            ncbi = (row.get("locus_tag") or "").strip()
            if cyan and ncbi:
                bridge[cyan] = ncbi
    return bridge


def load_joinable_locus_tags(gene_mapping: Path) -> set[str]:
    """NCBI locus tags the CyanoRak merge can actually join on.

    ``load_gff_from_ncbi_and_cyanorak`` merges on ``old_locus_tag`` (the NCBI
    ``locus_tag`` column is the *old* tag, NaN when a gene has none). A handful
    of RefSeq genes have no ``old_locus_tag`` (PGAP-added, RS tag only); a
    CyanoRak ORF bridged to such an RS tag cannot merge and would surface as a
    phantom CyanoRak-only duplicate. We therefore restrict bridging to genes
    that expose an ``old_locus_tags`` value."""
    joinable: set[str] = set()
    with open(gene_mapping) as f:
        for row in csv.DictReader(f):
            for o in (row.get("old_locus_tags") or "").replace(";", ",").split(","):
                if o.strip():
                    joinable.add(o.strip())
    return joinable


def build_gff_attributes(feat) -> dict[str, str]:
    """Map gbff CDS /note qualifiers to CyanoRak GFF attribute values."""
    notes = _notes(feat)
    quals = feat.qualifiers
    attrs: dict[str, str] = {}

    # ID — the CyanoRak ORF Id (join key matched by the gbk id-map)
    orf = _note_values(notes, "cyanorak ORF Id:")
    if orf:
        attrs["ID"] = orf[0]

    # Name — gene name (skip the 'None' sentinel CyanoRak uses)
    gene = quals.get("gene", [None])[0]
    if gene and gene != "None":
        attrs["Name"] = _enc(gene)
    else:
        prev = _note_values(notes, "Previous name:")
        if prev:
            attrs["Name"] = _enc(prev[0])

    # product
    prod = quals.get("product", [None])[0]
    if prod:
        attrs["product"] = _enc(prod)

    # cluster_number — the CK ortholog cluster (headline CyanoRak value)
    clu = _note_values(notes, "cyanorak cluster number:")
    if clu:
        attrs["cluster_number"] = clu[0]

    # CyanoRak role codes (comma-joined) + description(s)
    roles = _code_list(notes, "Cyanorak role:")
    if roles:
        attrs["cyanorak_Role"] = ",".join(roles)
    role_defs = [_enc(_strip_code(v)) for v in _note_values(notes, "Cyanorak role def.:")]
    if role_defs:
        attrs["cyanorak_Role_description"] = ",".join(role_defs)

    # TIGR role codes + description(s)
    tigr = _code_list(notes, "TIGR role:")
    if tigr:
        attrs["tIGR_Role"] = ",".join(tigr)
    tigr_defs = [_enc(_strip_code(v)) for v in _note_values(notes, "TIGR role def.:")]
    if tigr_defs:
        attrs["tIGR_Role_description"] = ",".join(tigr_defs)

    # GO terms (MF + BP + CC + Unclassified); gbff stores bare ids
    go_ids: list[str] = []
    go_descs: list[str] = []
    for ns in ("Molecular function", "Biological process",
               "Cellular component", "Unclassified"):
        for tok in _code_list(notes, f"GO terms {ns}:"):
            go_ids.append(f"GO:{tok}")
        go_descs.extend(_strip_code(v) for v in _note_values(notes, f"GO terms {ns} def.:"))
    if go_ids:
        attrs["Ontology_term"] = ",".join(go_ids)
    if go_descs:
        attrs["ontology_term_description"] = ",".join(_enc(d) for d in go_descs)

    # eggNOG-ish orthology ids (COG / bactNOG / cyaNOG / CyOG)
    egg: list[str] = []
    egg_desc: list[str] = []
    for key in ("COG", "bactNOG", "cyaNOG", "CyOG"):
        egg.extend(_code_list(notes, f"{key}:"))
        egg_desc.extend(_strip_code(v) for v in _note_values(notes, f"{key} def.:"))
    if egg:
        attrs["eggNOG"] = ",".join(egg)
    if egg_desc:
        attrs["eggNOG_description"] = ",".join(_enc(d) for d in egg_desc)

    # CyanoRak 'kegg' column actually stores EC numbers (see config comment)
    ec = _code_list(notes, "EC_number:")
    if ec:
        attrs["kegg"] = ",".join(ec)
    ec_desc = _note_values(notes, "EC_number def.:")
    if ec_desc:
        attrs["kegg_description"] = ",".join(_enc(_strip_code(d)) for d in ec_desc)

    # protein_domains — Pfam + TIGRFam + Interpro (all space-separated in gbff)
    domains: list[str] = []
    for prefix in ("pFam:", "TIGRFam:", "Interpro Domain:"):
        domains.extend(_code_list(notes, prefix))
    if domains:
        attrs["protein_domains"] = ",".join(domains)

    return attrs


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--gbff", type=Path, required=True,
                   help="CyanoRak GenBank export (.gbff)")
    p.add_argument("--bridge", type=Path, required=True,
                   help="diamond bridge CSV (<cyanorak_locus_tag>, locus_tag)")
    p.add_argument("--organism", required=True,
                   help="CyanoRak organism short name, e.g. Pro_MIT1327")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output cyanorak/ dir (e.g. cache/.../<Strain>/cyanorak)")
    p.add_argument("--gene-mapping", type=Path, required=True,
                   help="NCBI-only gene_mapping.csv (restricts bridging to "
                        "joinable old_locus_tags)")
    p.add_argument("--source-id-col", default="cyanorak_locus_tag",
                   help="source-id column name in the bridge CSV")
    args = p.parse_args()

    for f in (args.gbff, args.bridge, args.gene_mapping):
        if not f.exists():
            print(f"ERROR: {f} not found", file=sys.stderr)
            return 1

    bridge = load_bridge(args.bridge, args.source_id_col)
    joinable = load_joinable_locus_tags(args.gene_mapping)
    dropped_rs = {c: t for c, t in bridge.items() if t not in joinable}
    bridge = {c: t for c, t in bridge.items() if t in joinable}
    print(f"bridge entries (cyanorak -> NCBI locus_tag): {len(bridge)} joinable "
          f"(+{len(dropped_rs)} dropped: target has no old_locus_tag, unmergeable)")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    gff_path = args.out_dir / f"{args.organism}.gff"
    gbk_path = args.out_dir / f"{args.organism}.gbk"

    records = list(SeqIO.parse(str(args.gbff), "genbank"))
    print(f"gbff records (contigs): {len(records)}")

    # Pass 1: candidate CDS per NCBI locus_tag, longest-protein-wins.
    # value = (ncbi_lt, feat, record, attrs, protein_len)
    by_ncbi: dict[str, tuple] = {}
    n_cds = n_unmapped = n_dup = 0
    for rec in records:
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            n_cds += 1
            cyan_lt = feat.qualifiers.get("locus_tag", [None])[0]
            ncbi_lt = bridge.get(cyan_lt) if cyan_lt else None
            if not ncbi_lt:
                n_unmapped += 1
                continue
            attrs = build_gff_attributes(feat)
            if "ID" not in attrs:
                n_unmapped += 1
                continue
            plen = len(feat.qualifiers.get("translation", [""])[0])
            prev = by_ncbi.get(ncbi_lt)
            if prev is None:
                by_ncbi[ncbi_lt] = (ncbi_lt, feat, rec, attrs, plen)
            else:
                n_dup += 1
                if plen > prev[4]:
                    by_ncbi[ncbi_lt] = (ncbi_lt, feat, rec, attrs, plen)

    mapped = sorted(by_ncbi.values(), key=lambda t: t[0])
    print(f"CDS total: {n_cds}; mapped (deduped to NCBI locus_tag): {len(mapped)}; "
          f"unmapped CDS dropped: {n_unmapped}; collisions (longer won): {n_dup}")

    # ── write GFF ──
    with open(gff_path, "w") as gff:
        gff.write("##gff-version 3\n")
        for ncbi_lt, feat, rec, attrs, _plen in mapped:
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-"
            attr_str = ";".join(
                f"{k}={attrs[k]}" for k in GFF_COLS if attrs.get(k)
            )
            gff.write("\t".join([
                rec.id, GFF_SOURCE, "CDS", str(start), str(end),
                ".", strand, "0", attr_str,
            ]) + "\n")
    print(f"wrote {gff_path}")

    # ── write GBK (relabeled, winning mapped CDS only) ──
    winners = {id(t[1]) for t in mapped}
    ncbi_for_feat = {id(t[1]): t[0] for t in mapped}
    out_records = []
    for rec in records:
        new_feats = []
        for feat in rec.features:
            if feat.type == "source":
                new_feats.append(feat)
            elif feat.type == "CDS" and id(feat) in winners:
                feat.qualifiers["locus_tag"] = [ncbi_for_feat[id(feat)]]
                new_feats.append(feat)
        has_cds = any(f.type == "CDS" for f in new_feats)
        if has_cds:
            rec.features = new_feats
            out_records.append(rec)
    with open(gbk_path, "w") as gbk:
        SeqIO.write(out_records, gbk, "genbank")
    n_cds_out = sum(1 for r in out_records for f in r.features if f.type == "CDS")
    print(f"wrote {gbk_path} ({len(out_records)} records, {n_cds_out} CDS)")

    if n_cds_out != len(mapped):
        print(f"WARNING: gbk CDS count {n_cds_out} != mapped {len(mapped)}",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
