#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


GTEX_TISSUE_TO_CANCER_SITE = {
    "Bladder": ["bladder"],
    "Brain_Cerebellum": ["brain_cns"],
    "Brain_Cortex": ["brain_cns"],
    "Brain_Frontal_Cortex_BA9": ["brain_cns"],
    "Brain_Nucleus_accumbens_basal_ganglia": ["brain_cns"],
    "Breast_Mammary_Tissue": ["female_breast"],
    "Colon_Sigmoid": ["colorectal"],
    "Colon_Transverse": ["colorectal"],
    "Esophagus_Gastroesophageal_Junction": ["esophagus"],
    "Esophagus_Mucosa": ["esophagus"],
    "Esophagus_Muscularis": ["esophagus"],
    "Kidney_Cortex": ["kidney_renal_pelvis"],
    "Liver": ["liver_intrahepatic_bile_ducts"],
    "Lung": ["lung_trachea_bronchus"],
    "Ovary": ["ovary"],
    "Pancreas": ["pancreas"],
    "Prostate": ["prostate"],
    "Skin_Not_Sun_Exposed_Suprapubic": ["melanoma_skin", "nonmelanoma_skin"],
    "Skin_Sun_Exposed_Lower_leg": ["melanoma_skin", "nonmelanoma_skin"],
    "Stomach": ["stomach"],
    "Testis": ["testicular"],
    "Thyroid": ["thyroid"],
    "Uterus": ["corpus_uteri"],
    "Vagina": ["corpus_uteri"],
    "Whole_Blood": ["leukemia", "non_hodgkin_lymphoma", "multiple_myeloma_immunoproliferative"],
    "Cells_EBV-transformed_lymphocytes": [
        "leukemia",
        "non_hodgkin_lymphoma",
        "multiple_myeloma_immunoproliferative",
    ],
}


def split_list(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def split_float_list(value: object) -> list[float]:
    out = []
    for item in split_list(value):
        try:
            out.append(float(item))
        except ValueError:
            out.append(float("nan"))
    return out


def split_variant_id(value: object) -> tuple[object, object, object, object]:
    if pd.isna(value):
        return pd.NA, pd.NA, pd.NA, pd.NA
    parts = str(value).split(":")
    if len(parts) < 4:
        return pd.NA, pd.NA, pd.NA, pd.NA
    return parts[0], parts[1], parts[2], parts[3]


def normalize_direction(value: object) -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value in {"+", "plus", "positive"}:
        return "+"
    if value in {"-", "minus", "negative"}:
        return "-"
    return value


def cancer_site_matches(tissues: list[str]) -> str:
    sites = set()
    for tissue in tissues:
        sites.update(GTEX_TISSUE_TO_CANCER_SITE.get(tissue, []))
    return ";".join(sorted(sites))


def read_st26(path: Path) -> pd.DataFrame:
    st26 = pd.read_excel(path, sheet_name="ST26", header=3, dtype=object)
    st26 = st26.rename(
        columns={
            "UKBPPP ProteinID": "protein_id",
            "Top variant ID": "variant_id_hg37",
            "Direction (ALT)": "pqtl_alt_direction",
            "Gene name": "gene_symbol",
            "Coloc count": "gtex_coloc_tissue_count",
            "Coloc tissues": "gtex_coloc_tissues",
            "Direction (pQTL top variant, ALT)": "gtex_eqtl_alt_directions",
            "Top variant ID.1": "gtex_top_variant_ids_hg38",
            "Coloc PP.H4": "gtex_coloc_pph4_values",
        }
    )
    return st26[st26["protein_id"].notna() & st26["variant_id_hg37"].notna()].copy()


def expand_long(st26: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in st26.itertuples(index=False):
        tissues = split_list(row.gtex_coloc_tissues)
        directions = split_list(row.gtex_eqtl_alt_directions)
        top_variants = split_list(row.gtex_top_variant_ids_hg38)
        pph4_values = split_float_list(row.gtex_coloc_pph4_values)
        n = max(len(tissues), len(directions), len(top_variants), len(pph4_values))
        for i in range(n):
            rows.append(
                {
                    "protein_id": row.protein_id,
                    "variant_id_hg37": row.variant_id_hg37,
                    "gene_symbol": row.gene_symbol,
                    "pqtl_alt_direction": normalize_direction(row.pqtl_alt_direction),
                    "gtex_tissue": tissues[i] if i < len(tissues) else "",
                    "gtex_eqtl_alt_direction": normalize_direction(directions[i]) if i < len(directions) else "",
                    "gtex_top_variant_id_hg38": top_variants[i] if i < len(top_variants) else "",
                    "gtex_coloc_pph4": pph4_values[i] if i < len(pph4_values) else float("nan"),
                }
            )
    long = pd.DataFrame(rows)
    if long.empty:
        return long
    long["direction_concordant"] = (
        long["pqtl_alt_direction"].ne("")
        & long["gtex_eqtl_alt_direction"].ne("")
        & long["pqtl_alt_direction"].eq(long["gtex_eqtl_alt_direction"])
    ).astype(int)
    long["mapped_cancer_sites"] = long["gtex_tissue"].map(lambda x: cancer_site_matches([x]))
    return long


def summarize_st26(st26: pd.DataFrame, long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame()
    tissue_counts = long.groupby(["protein_id", "variant_id_hg37"])["gtex_tissue"].nunique()
    max_tissues = max(int(tissue_counts.max()), 1)
    grouped = (
        long.groupby(["protein_id", "variant_id_hg37"], dropna=False)
        .agg(
            gene_symbol=("gene_symbol", "first"),
            pqtl_alt_direction=("pqtl_alt_direction", "first"),
            gtex_coloc_tissue_count=("gtex_tissue", "nunique"),
            gtex_coloc_tissues=("gtex_tissue", lambda x: ";".join(sorted(set(filter(None, map(str, x)))))),
            gtex_coloc_max_pph4=("gtex_coloc_pph4", "max"),
            gtex_coloc_mean_pph4=("gtex_coloc_pph4", "mean"),
            gtex_direction_concordant_fraction=("direction_concordant", "mean"),
            gtex_cancer_site_matches=("mapped_cancer_sites", lambda x: ";".join(sorted(set(";".join(x).split(";")) - {""}))),
        )
        .reset_index()
    )
    grouped["gtex_tissue_specificity_score"] = 1.0 - (
        grouped["gtex_coloc_tissue_count"].clip(lower=1) - 1
    ) / max(max_tissues - 1, 1)
    grouped["snp_tissue_reg_prior"] = grouped["gtex_coloc_max_pph4"].fillna(0) * (
        0.65 + 0.35 * grouped["gtex_tissue_specificity_score"]
    )
    grouped["snp_tissue_reg_prior"] = grouped["snp_tissue_reg_prior"].clip(lower=0, upper=1)
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, default=Path("文献资料补充/Nature-PQTL-补充材料.xlsx"))
    parser.add_argument("--edges", type=Path, default=Path("model_data/st10_core/st10_pqtl_edges.csv"))
    parser.add_argument("--out", type=Path, default=Path("model_data/st10_core/st10_snp_tissue_regulatory_priors.csv"))
    parser.add_argument(
        "--out-long",
        type=Path,
        default=Path("model_data/st10_core/st10_snp_tissue_regulatory_priors.long.csv"),
    )
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    edges = pd.read_csv(args.edges, dtype={"rsid": str})
    st26 = read_st26(args.xlsx)
    long = expand_long(st26)
    summary = summarize_st26(st26, long)

    chr_pos_ref_alt = summary["variant_id_hg37"].map(split_variant_id)
    summary[["variant_chr_hg37", "variant_pos_hg37", "variant_ref", "variant_alt"]] = pd.DataFrame(
        chr_pos_ref_alt.tolist(), index=summary.index
    )

    keep_edge = [
        "protein_id",
        "variant_id_hg37",
        "rsid",
        "cis_trans",
        "gene_symbol",
        "regulatory_mode",
        "annotated_gene",
    ]
    matched = edges[keep_edge].merge(
        summary,
        on=["protein_id", "variant_id_hg37"],
        how="left",
        suffixes=("", "_st26"),
    )
    matched["has_gtex_eqtl_coloc"] = matched["gtex_coloc_tissue_count"].notna().astype(int)
    for col in [
        "gtex_coloc_tissue_count",
        "gtex_coloc_max_pph4",
        "gtex_coloc_mean_pph4",
        "gtex_direction_concordant_fraction",
        "gtex_tissue_specificity_score",
        "snp_tissue_reg_prior",
    ]:
        matched[col] = pd.to_numeric(matched[col], errors="coerce").fillna(0.0)
    for col in ["gtex_coloc_tissues", "gtex_cancer_site_matches"]:
        matched[col] = matched[col].fillna("")

    long = long.merge(
        edges[["protein_id", "variant_id_hg37", "rsid"]].drop_duplicates(),
        on=["protein_id", "variant_id_hg37"],
        how="inner",
    )
    matched.to_csv(args.out, index=False)
    long.to_csv(args.out_long, index=False)

    summary_path = args.out.with_suffix(".summary.txt")
    with summary_path.open("w") as out:
        out.write(f"st10_edges: {len(edges)}\n")
        out.write(f"matched_edges_with_gtex_eqtl_coloc: {int(matched['has_gtex_eqtl_coloc'].sum())}\n")
        out.write(f"matched_proteins: {matched.loc[matched['has_gtex_eqtl_coloc'].eq(1), 'protein_id'].nunique()}\n")
        out.write(f"matched_snps: {matched.loc[matched['has_gtex_eqtl_coloc'].eq(1), 'rsid'].nunique()}\n")
        out.write(f"long_tissue_rows: {len(long)}\n")
        out.write(f"mean_snp_tissue_reg_prior: {matched['snp_tissue_reg_prior'].mean():.4f}\n")
        out.write(f"mean_prior_among_matches: {matched.loc[matched['has_gtex_eqtl_coloc'].eq(1), 'snp_tissue_reg_prior'].mean():.4f}\n")
        out.write("\nTop GTEx tissues among matched ST10 edges:\n")
        out.write(long["gtex_tissue"].value_counts().head(25).to_string())
        out.write("\n\nMapped cancer sites among matched ST10 edges:\n")
        sites = long["mapped_cancer_sites"].str.split(";", expand=True).stack()
        sites = sites[sites.ne("")]
        out.write(sites.value_counts().head(25).to_string())
        out.write("\n")
    print(summary_path.read_text())


if __name__ == "__main__":
    main()
