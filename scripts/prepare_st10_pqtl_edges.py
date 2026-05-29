#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CODING_CONSEQUENCES = {
    "missense_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "protein_altering_variant",
    "inframe_insertion",
    "inframe_deletion",
    "synonymous_variant",
}

REGULATORY_CONSEQUENCES = {
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
    "splice_region_variant",
    "regulatory_region_variant",
    "TF_binding_site_variant",
}


def split_consequences(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {x.strip() for x in str(value).split(",") if x.strip()}


def classify_consequence(value: object) -> str:
    terms = split_consequences(value)
    if terms & CODING_CONSEQUENCES:
        return "coding_or_splice"
    if terms & REGULATORY_CONSEQUENCES:
        return "regulatory_or_intronic"
    if terms:
        return "other_annotated"
    return "unknown"


def classify_regulatory_mode(row: pd.Series) -> str:
    consequence_class = row["consequence_class"]
    same_gene = str(row.get("annotated_gene", "")).upper() == str(row.get("gene_symbol", "")).upper()
    cis_trans = str(row.get("cis_trans", "")).lower()
    if cis_trans == "cis" and same_gene and consequence_class == "coding_or_splice":
        return "cis_target_coding_or_splice"
    if cis_trans == "cis" and same_gene:
        return "cis_target_regulatory_or_linked"
    if cis_trans == "cis":
        return "cis_locus_linked"
    if same_gene and consequence_class == "coding_or_splice":
        return "trans_annotated_gene_coding_or_splice"
    if consequence_class == "coding_or_splice":
        return "trans_other_gene_coding_or_splice"
    if consequence_class == "regulatory_or_intronic":
        return "trans_regulatory_or_intronic"
    return "trans_other_or_unknown"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, default=Path("文献资料补充/Nature-PQTL-补充材料.xlsx"))
    parser.add_argument("--matched-rsid-dir", type=Path, default=Path("st10_core_rsid/matched_ids"))
    parser.add_argument("--outdir", type=Path, default=Path("model_data/st10_core"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    matched = []
    for path in sorted(args.matched_rsid_dir.glob("chr*.matched_rsid.txt")):
        chrom = path.name.split(".")[0].replace("chr", "")
        with path.open() as fh:
            for line in fh:
                rsid = line.strip()
                if rsid:
                    matched.append((rsid, chrom))
    matched_df = pd.DataFrame(matched, columns=["rsid", "matched_chr"]).drop_duplicates()

    st10 = pd.read_excel(args.xlsx, sheet_name="ST10", header=4, dtype=object)
    st3 = pd.read_excel(args.xlsx, sheet_name="ST3", header=2, dtype=object)

    keep = [
        "Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)",
        "CHROM",
        "GENPOS (hg38)",
        "UKBPPP ProteinID",
        "Assay Target",
        "Target UniProt",
        "rsID",
        "A1FREQ",
        "BETA",
        "SE",
        "log10(p)",
        "cis/trans",
        "Bioinfomatic annotated gene",
        "Annotated gene consequence",
        "IMPACT",
    ]
    st10 = st10[keep].rename(
        columns={
            "Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)": "variant_id_hg37",
            "CHROM": "chrom",
            "GENPOS (hg38)": "pos_hg38",
            "UKBPPP ProteinID": "protein_id",
            "Assay Target": "assay_target",
            "Target UniProt": "target_uniprot",
            "rsID": "rsid",
            "A1FREQ": "a1freq",
            "BETA": "beta",
            "SE": "se",
            "log10(p)": "log10p",
            "cis/trans": "cis_trans",
            "Bioinfomatic annotated gene": "annotated_gene",
            "Annotated gene consequence": "consequence",
            "IMPACT": "impact",
        }
    )
    st10 = st10[st10["rsid"].notna() & (st10["rsid"].astype(str) != "-")].copy()
    st10["rsid"] = st10["rsid"].astype(str)
    st10 = st10.merge(matched_df, on="rsid", how="inner")

    protein_anno = st3[
        ["UKBPPP ProteinID", "Protein panel", "Gene symbol", "UniProt", "Gene CHROM", "Gene start", "Gene end"]
    ].drop_duplicates()
    protein_anno = protein_anno.rename(
        columns={
            "UKBPPP ProteinID": "protein_id",
            "Protein panel": "protein_panel",
            "Gene symbol": "gene_symbol",
            "UniProt": "uniprot",
            "Gene CHROM": "gene_chrom",
            "Gene start": "gene_start",
            "Gene end": "gene_end",
        }
    )
    edges = st10.merge(protein_anno, on="protein_id", how="left")
    for col in ["beta", "se", "log10p", "a1freq"]:
        edges[col] = pd.to_numeric(edges[col], errors="coerce")

    edges["abs_beta"] = edges["beta"].abs()
    edges["beta_logp_weight"] = edges["beta"] * edges["log10p"].fillna(0).pow(0.5)
    edges["is_cis"] = edges["cis_trans"].eq("cis").astype(int)
    edges["is_trans"] = edges["cis_trans"].eq("trans").astype(int)
    edges["consequence_class"] = edges["consequence"].map(classify_consequence)
    edges["annotated_gene_matches_target"] = (
        edges["annotated_gene"].fillna("").astype(str).str.upper()
        == edges["gene_symbol"].fillna("").astype(str).str.upper()
    ).astype(int)
    edges["regulatory_mode"] = edges.apply(classify_regulatory_mode, axis=1)
    snp_n_proteins = edges.groupby("rsid")["protein_id"].transform("nunique")
    snp_n_panels = edges.groupby("rsid")["protein_panel"].transform("nunique")
    edges["snp_n_target_proteins"] = snp_n_proteins.astype(int)
    edges["snp_n_target_panels"] = snp_n_panels.astype(int)
    edges["snp_is_pleiotropic_pqtl"] = (snp_n_proteins >= 5).astype(int)
    edges = edges.sort_values(["protein_id", "cis_trans", "rsid"]).reset_index(drop=True)

    protein_summary = (
        edges.groupby("protein_id")
        .agg(
            assay_target=("assay_target", "first"),
            protein_panel=("protein_panel", "first"),
            gene_symbol=("gene_symbol", "first"),
            uniprot=("uniprot", "first"),
            n_edges=("rsid", "size"),
            n_snps=("rsid", "nunique"),
            n_cis=("is_cis", "sum"),
            n_trans=("is_trans", "sum"),
            max_log10p=("log10p", "max"),
            sum_abs_beta=("abs_beta", "sum"),
            n_regulatory_modes=("regulatory_mode", "nunique"),
            regulatory_modes=("regulatory_mode", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
        .sort_values(["protein_panel", "assay_target", "protein_id"])
    )

    snp_summary = (
        edges.groupby("rsid")
        .agg(
            chrom=("matched_chr", "first"),
            n_proteins=("protein_id", "nunique"),
            n_edges=("protein_id", "size"),
            max_abs_beta=("abs_beta", "max"),
            max_log10p=("log10p", "max"),
            cis_trans_values=("cis_trans", lambda x: ";".join(sorted(set(map(str, x))))),
            regulatory_modes=("regulatory_mode", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
        .sort_values(["chrom", "rsid"])
    )

    edges.to_csv(args.outdir / "st10_pqtl_edges.csv", index=False)
    protein_summary.to_csv(args.outdir / "st10_protein_manifest.csv", index=False)
    snp_summary.to_csv(args.outdir / "st10_snp_manifest.csv", index=False)
    matched_df.to_csv(args.outdir / "matched_rsid_by_chr.csv", index=False)

    with (args.outdir / "st10_pqtl_edge_summary.txt").open("w") as out:
        out.write(f"Matched rsIDs: {matched_df['rsid'].nunique()}\n")
        out.write(f"ST10 matched edges: {len(edges)}\n")
        out.write(f"Proteins with matched edges: {edges['protein_id'].nunique()}\n")
        out.write(f"Unique rsIDs in edges: {edges['rsid'].nunique()}\n")
        out.write("\nEdges by cis/trans:\n")
        out.write(edges["cis_trans"].value_counts(dropna=False).to_string())
        out.write("\n\nProteins by panel:\n")
        out.write(protein_summary["protein_panel"].value_counts(dropna=False).to_string())
        out.write("\n")


if __name__ == "__main__":
    main()
