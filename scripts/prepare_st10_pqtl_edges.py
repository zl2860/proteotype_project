#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


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
