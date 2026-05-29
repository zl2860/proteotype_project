#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


HEADER_ROWS = {
    "ST12": 2,
    "ST13": 2,
    "ST14": 2,
    "ST15": 3,
    "ST16": 3,
    "ST18": 3,
    "ST19": 2,
    "ST20": 2,
    "ST21": 4,
    "ST23": 3,
    "ST24": 3,
    "ST25": 2,
    "ST27": 4,
    "ST29": 3,
}

NUMERIC_FEATURE_PREFIXES = (
    "base_",
    "variant_",
    "edge_",
    "protein_",
    "hpa_",
    "gtex_",
    "st",
)


def clean_col(col: object) -> str:
    return re.sub(r"\s+", " ", str(col).strip())


def read_sheet(xlsx: Path, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx, sheet_name=sheet, header=HEADER_ROWS[sheet], dtype=object)
    df.columns = [clean_col(c) for c in df.columns]
    return df.dropna(how="all").reset_index(drop=True)


def to_num(series: pd.Series | object, default: float = 0.0) -> pd.Series | float:
    if isinstance(series, pd.Series):
        return pd.to_numeric(series, errors="coerce").fillna(default)
    try:
        return float(series)
    except (TypeError, ValueError):
        return default


def split_values(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in re.split(r"[,;]\s*", str(value)) if x.strip() and x.strip() not in {"-", "nan"}]


def split_float_values(value: object) -> list[float]:
    out = []
    for item in split_values(value):
        try:
            out.append(float(item))
        except ValueError:
            pass
    return out


def normalize_variant_no_suffix(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text in {"-", "nan", "NaN"}:
        return ""
    parts = text.split(":")
    return ":".join(parts[:4]) if len(parts) >= 4 else text


def add_indicator(df: pd.DataFrame, col: str, prefix: str, values: list[str] | None = None) -> list[str]:
    if values is None:
        values = sorted(v for v in df[col].dropna().astype(str).unique() if v and v != "nan")
    out_cols = []
    for value in values:
        safe = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()
        out_col = f"{prefix}_{safe}"
        df[out_col] = df[col].fillna("").astype(str).eq(value).astype(np.float32)
        out_cols.append(out_col)
    return out_cols


def aggregate_evidence(catalog: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    catalog = catalog.copy()
    catalog["variant_key_norm"] = catalog["variant_id"].map(normalize_variant_no_suffix)
    catalog.loc[catalog["variant_key_norm"].eq(""), "variant_key_norm"] = catalog["rsid"].fillna("")
    catalog["edge_key"] = catalog["ukbppp_protein_id"].fillna("") + "|" + catalog["variant_key_norm"]
    evidence_types = sorted(catalog["evidence_type"].dropna().astype(str).unique())

    variant = catalog[catalog["variant_key_norm"].ne("")].groupby("variant_key_norm").agg(
        variant_evidence_count=("evidence_type", "nunique"),
        variant_source_sheet_count=("source_sheet", "nunique"),
    )
    for evidence in evidence_types:
        hits = catalog.loc[catalog["evidence_type"].eq(evidence) & catalog["variant_key_norm"].ne(""), "variant_key_norm"]
        variant[f"variant_ev_{evidence}"] = variant.index.isin(set(hits)).astype(np.float32)
    variant = variant.reset_index()

    edge_catalog = catalog[catalog["ukbppp_protein_id"].notna() & catalog["variant_key_norm"].ne("")]
    edge = edge_catalog.groupby("edge_key").agg(
        edge_evidence_count=("evidence_type", "nunique"),
        edge_source_sheet_count=("source_sheet", "nunique"),
    )
    for evidence in evidence_types:
        hits = edge_catalog.loc[edge_catalog["evidence_type"].eq(evidence), "edge_key"]
        edge[f"edge_ev_{evidence}"] = edge.index.isin(set(hits)).astype(np.float32)
    edge = edge.reset_index()
    return variant, edge


def build_fine_mapping(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST16").iloc[1:].reset_index(drop=True)
    rows = []
    for _, row in df.iterrows():
        protein_id = row.get("UKBPPP ProteinID")
        top_variant = normalize_variant_no_suffix(row.get("Top variant (defined by highest PIP)"))
        top_pip = to_num(row.get("Unnamed: 3"), np.nan)
        top_log10bf = to_num(row.get("Unnamed: 4"), np.nan)
        top_beta = to_num(row.get("Unnamed: 5"), np.nan)
        top_cond_logp = to_num(row.get("Unnamed: 7"), np.nan)
        cs_size = to_num(row.get("Credible set (truncated to top 1000 variants)"), np.nan)
        cs_variants = split_values(row.get("Unnamed: 13"))
        cs_pips = split_float_values(row.get("Unnamed: 15"))
        for i, variant in enumerate(cs_variants):
            variant_norm = normalize_variant_no_suffix(variant)
            if not variant_norm:
                continue
            rows.append(
                {
                    "protein_id": protein_id,
                    "variant_key_norm": variant_norm,
                    "st16_is_finemap_top": float(variant_norm == top_variant),
                    "st16_is_credible_set": 1.0,
                    "st16_top_pip": top_pip if variant_norm == top_variant else 0.0,
                    "st16_top_log10bf": top_log10bf if variant_norm == top_variant else 0.0,
                    "st16_top_beta_cond": top_beta if variant_norm == top_variant else 0.0,
                    "st16_top_cond_log10p": top_cond_logp if variant_norm == top_variant else 0.0,
                    "st16_cs_size": cs_size,
                    "st16_cs_variant_pip": cs_pips[i] if i < len(cs_pips) else np.nan,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["edge_key"] = out["protein_id"].fillna("") + "|" + out["variant_key_norm"]
    return (
        out.groupby("edge_key")
        .agg(
            st16_is_finemap_top=("st16_is_finemap_top", "max"),
            st16_is_credible_set=("st16_is_credible_set", "max"),
            st16_top_pip=("st16_top_pip", "max"),
            st16_top_log10bf=("st16_top_log10bf", "max"),
            st16_top_beta_cond=("st16_top_beta_cond", "first"),
            st16_top_cond_log10p=("st16_top_cond_log10p", "max"),
            st16_min_cs_size=("st16_cs_size", "min"),
            st16_max_cs_variant_pip=("st16_cs_variant_pip", "max"),
        )
        .reset_index()
    )


def build_st12_features(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST12")
    df["variant_key_norm"] = df["Variant ID"].map(normalize_variant_no_suffix)
    df["st12_is_high_impact_coding"] = 1.0
    df["st12_sift_score"] = to_num(df["SIFT_score"], np.nan)
    df["st12_fitcons_rankscore"] = to_num(df["integrated_fitCons_rankscore"], np.nan)
    df["st12_fitcons_score"] = to_num(df["integrated_fitCons_score"], np.nan)
    return df[
        [
            "variant_key_norm",
            "Gene",
            "Consequence",
            "IMPACT",
            "st12_is_high_impact_coding",
            "st12_sift_score",
            "st12_fitcons_rankscore",
            "st12_fitcons_score",
        ]
    ].rename(columns={"Gene": "st12_gene", "Consequence": "st12_consequence", "IMPACT": "st12_impact"})


def build_st13_features(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST13")
    df["variant_key_norm"] = df["Variant ID"].map(normalize_variant_no_suffix)
    df["st13_is_regulatory_noncoding"] = 1.0
    df["st13_motif_score_change"] = to_num(df["MOTIF_SCORE_CHANGE"], 0.0)
    df["st13_abs_motif_score_change"] = df["st13_motif_score_change"].abs()
    df["st13_has_tf_motif"] = df["MOTIF_NAME"].fillna("").astype(str).ne("-").astype(np.float32)
    df["st13_tf_count"] = df["TRANSCRIPTION_FACTORS"].map(lambda x: len(split_values(x))).astype(np.float32)
    return (
        df.groupby("variant_key_norm")
        .agg(
            st13_is_regulatory_noncoding=("st13_is_regulatory_noncoding", "max"),
            st13_abs_motif_score_change=("st13_abs_motif_score_change", "max"),
            st13_has_tf_motif=("st13_has_tf_motif", "max"),
            st13_tf_count=("st13_tf_count", "max"),
            st13_biotype=("BIOTYPE", "first"),
        )
        .reset_index()
    )


def build_st14_features(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST14")
    rows = []
    for _, row in df.iterrows():
        protein_id = row.get("UKBPPP ProteinID")
        for col, flag in [("rsID (lead variant)", "st14_is_lead_in_ld_with_pav"), ("rsID (PAV)", "st14_is_pav")]:
            rsid = row.get(col)
            if pd.isna(rsid):
                continue
            rows.append(
                {
                    "protein_id": protein_id,
                    "rsid": str(rsid),
                    flag: 1.0,
                    "st14_lead_is_pav": to_num(row.get("Lead variant is PAV"), 0.0),
                    "st14_pav_ld_r2": to_num(row.get("R2 with lead variant"), 0.0),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["edge_key_rsid"] = out["protein_id"].fillna("") + "|" + out["rsid"].fillna("")
    return out.groupby("edge_key_rsid").max(numeric_only=True).reset_index()


def build_st15_features(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST15")
    df["variant_key_norm"] = df["Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"].map(normalize_variant_no_suffix)
    df["edge_key"] = df["UKBPPP ProteinID"].fillna("") + "|" + df["variant_key_norm"]
    df["st15_is_replicated"] = df["novel/replicated"].fillna("").astype(str).eq("replicated").astype(np.float32)
    df["st15_is_novel"] = df["novel/replicated"].fillna("").astype(str).eq("novel").astype(np.float32)
    df["st15_has_orthogonal_tech"] = df["technology if replicated"].fillna("").astype(str).ne("").astype(np.float32)
    return df.groupby("edge_key").agg(
        st15_is_replicated=("st15_is_replicated", "max"),
        st15_is_novel=("st15_is_novel", "max"),
        st15_has_orthogonal_tech=("st15_has_orthogonal_tech", "max"),
    ).reset_index()


def build_st18_features(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST18")
    df["variant_key_norm"] = df["Top variant ID"].map(normalize_variant_no_suffix)
    df["edge_key"] = df["UKBPPP ProteinID"].fillna("") + "|" + df["variant_key_norm"]
    df["st18_pqtl_coloc_count"] = to_num(df["Coloc count"], 0.0)
    df["st18_pqtl_coloc_max_pph4"] = df["Coloc PP.H4"].map(lambda x: max(split_float_values(x) or [0.0]))
    df["st18_pqtl_coloc_mean_pph4"] = df["Coloc PP.H4"].map(lambda x: float(np.mean(split_float_values(x))) if split_float_values(x) else 0.0)
    return df.groupby("edge_key").agg(
        st18_pqtl_coloc_count=("st18_pqtl_coloc_count", "max"),
        st18_pqtl_coloc_max_pph4=("st18_pqtl_coloc_max_pph4", "max"),
        st18_pqtl_coloc_mean_pph4=("st18_pqtl_coloc_mean_pph4", "max"),
    ).reset_index()


def build_network_features(xlsx: Path) -> pd.DataFrame:
    frames = []
    for sheet, variant_col, protein_col, prefix in [
        ("ST20", "ID", "UKBPPP ProteinID", "st20_ppi"),
        ("ST21", "ID", "UKBPPP ProteinID", "st21_reciprocal_trans"),
        ("ST23", "Variant ID", "Target protein (UKBPPP ProteinID)", "st23_receptor_ligand"),
    ]:
        df = read_sheet(xlsx, sheet)
        df["variant_key_norm"] = df[variant_col].map(normalize_variant_no_suffix)
        df["edge_key"] = df[protein_col].fillna("") + "|" + df["variant_key_norm"]
        df[f"{prefix}_flag"] = 1.0
        if sheet in {"ST21", "ST23"}:
            beta_col = "BETA"
            logp_col = "LOG10P" if "LOG10P" in df.columns else "log10(p)"
            df[f"{prefix}_abs_beta"] = to_num(df[beta_col], 0.0).abs()
            df[f"{prefix}_log10p"] = to_num(df[logp_col], 0.0)
        frames.append(df[["edge_key"] + [c for c in df.columns if c.startswith(prefix)]])
    out = pd.concat(frames, ignore_index=True, sort=False).fillna(0)
    return out.groupby("edge_key").max(numeric_only=True).reset_index()


def build_st24_features(xlsx: Path) -> pd.DataFrame:
    df = read_sheet(xlsx, "ST24").iloc[1:].reset_index(drop=True)
    df["variant_key_norm"] = df["Variant ID"].map(normalize_variant_no_suffix)
    df["edge_key"] = df["UKBPPP ProteinID"].fillna("") + "|" + df["variant_key_norm"]
    df["st24_discovery_log10p"] = to_num(df["Unnamed: 11"], 0.0)
    df["st24_blood_adj_log10p"] = to_num(df["Unnamed: 14"], 0.0)
    df["st24_bmi_adj_log10p"] = to_num(df["Unnamed: 17"], 0.0)
    df["st24_season_fasting_adj_log10p"] = to_num(df["Unnamed: 20"], 0.0)
    adj = df[["st24_blood_adj_log10p", "st24_bmi_adj_log10p", "st24_season_fasting_adj_log10p"]]
    df["st24_min_adj_log10p"] = adj.min(axis=1)
    df["st24_min_adj_over_discovery"] = df["st24_min_adj_log10p"] / df["st24_discovery_log10p"].replace(0, np.nan)
    df["st24_min_adj_over_discovery"] = df["st24_min_adj_over_discovery"].replace([np.inf, -np.inf], np.nan).fillna(0)
    return df[
        [
            "edge_key",
            "st24_discovery_log10p",
            "st24_blood_adj_log10p",
            "st24_bmi_adj_log10p",
            "st24_season_fasting_adj_log10p",
            "st24_min_adj_log10p",
            "st24_min_adj_over_discovery",
        ]
    ]


def build_protein_features(xlsx: Path, hpa_path: Path) -> pd.DataFrame:
    st19 = read_sheet(xlsx, "ST19").iloc[1:].reset_index(drop=True).rename(
        columns={
            "UKBPPP ProteinID": "protein_id",
            "pQTL component": "protein_st19_cis_h2",
            "Unnamed: 2": "protein_st19_trans_h2",
            "Unnamed: 3": "protein_st19_all_pqtl_h2",
            "Polygenic component": "protein_st19_polygenic_h2",
            "Total heritability (TH)": "protein_st19_total_h2",
            "Proportion of TH explained by primary cis pQTL component": "protein_st19_cis_h2_fraction",
            "Proportion of TH explained by primary trans pQTL component": "protein_st19_trans_h2_fraction",
        }
    )
    st19 = st19[st19["protein_id"].notna()].copy()
    for col in [c for c in st19.columns if c.startswith("protein_st19_")]:
        st19[col] = to_num(st19[col], 0.0)

    if hpa_path.exists():
        hpa = pd.read_csv(hpa_path).rename(
            columns={
                "hpa_has_match": "protein_hpa_has_match",
                "rna_tissue_prior": "protein_hpa_rna_tissue_prior",
                "protein_tissue_prior": "protein_hpa_protein_tissue_prior",
                "hpa_tissue_prior": "protein_hpa_tissue_prior",
            }
        )
        keep = [
            "protein_id",
            "protein_hpa_has_match",
            "protein_hpa_rna_tissue_prior",
            "protein_hpa_protein_tissue_prior",
            "protein_hpa_tissue_prior",
            "Secretome location",
            "Protein class",
        ]
        hpa = hpa[[c for c in keep if c in hpa.columns]]
        hpa["protein_hpa_is_secreted"] = hpa.get("Secretome location", "").fillna("").astype(str).str.contains(
            "secreted|blood|extracellular", case=False, regex=True
        ).astype(np.float32)
        hpa["protein_hpa_is_membrane"] = hpa.get("Protein class", "").fillna("").astype(str).str.contains(
            "membrane", case=False, regex=False
        ).astype(np.float32)
        st19 = st19.merge(hpa.drop(columns=[c for c in ["Secretome location", "Protein class"] if c in hpa]), on="protein_id", how="left")

    return st19


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, default=Path("文献资料补充/Nature-PQTL-补充材料.xlsx"))
    parser.add_argument("--edges", type=Path, default=Path("model_data/st10_core/st10_pqtl_edges.csv"))
    parser.add_argument("--catalog", type=Path, default=Path("文献资料补充/pqtl_variant_catalog/pqtl_variant_catalog.csv"))
    parser.add_argument("--hpa-priors", type=Path, default=Path("model_data/st10_core/st10_protein_tissue_priors.csv"))
    parser.add_argument("--snp-tissue-priors", type=Path, default=Path("model_data/st10_core/st10_snp_tissue_regulatory_priors.csv"))
    parser.add_argument("--outdir", type=Path, default=Path("model_data/st10_core/features"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    edges = pd.read_csv(args.edges, dtype={"rsid": str})
    edges["variant_key_norm"] = edges["variant_id_hg37"].map(normalize_variant_no_suffix)
    edges["edge_key"] = edges["protein_id"].fillna("") + "|" + edges["variant_key_norm"]
    edges["edge_key_rsid"] = edges["protein_id"].fillna("") + "|" + edges["rsid"].fillna("")

    edge_features = edges[
        [
            "protein_id",
            "variant_id_hg37",
            "variant_key_norm",
            "rsid",
            "assay_target",
            "gene_symbol",
            "cis_trans",
            "regulatory_mode",
            "consequence_class",
            "impact",
            "edge_key",
            "edge_key_rsid",
        ]
    ].copy()
    edge_features["base_beta"] = to_num(edges["beta"], 0.0)
    edge_features["base_abs_beta"] = edge_features["base_beta"].abs()
    edge_features["base_se"] = to_num(edges["se"], 0.0)
    edge_features["base_log10p"] = to_num(edges["log10p"], 0.0)
    edge_features["base_log1p_log10p"] = np.log1p(edge_features["base_log10p"].clip(lower=0))
    edge_features["base_a1freq"] = to_num(edges["a1freq"], 0.0)
    edge_features["base_is_cis"] = edges["cis_trans"].eq("cis").astype(np.float32)
    edge_features["base_is_trans"] = edges["cis_trans"].eq("trans").astype(np.float32)
    edge_features["base_annotated_gene_matches_target"] = to_num(edges["annotated_gene_matches_target"], 0.0)
    edge_features["base_snp_n_target_proteins_log1p"] = np.log1p(to_num(edges["snp_n_target_proteins"], 0.0))
    edge_features["base_snp_is_pleiotropic_pqtl"] = to_num(edges["snp_is_pleiotropic_pqtl"], 0.0)
    add_indicator(edge_features, "regulatory_mode", "base_regmode")
    add_indicator(edge_features, "consequence_class", "base_consequence")
    add_indicator(edge_features, "impact", "base_impact")

    catalog = pd.read_csv(args.catalog, dtype={"rsid": str}, low_memory=False)
    variant_ev, edge_ev = aggregate_evidence(catalog)
    edge_features = edge_features.merge(variant_ev, on="variant_key_norm", how="left")
    edge_features = edge_features.merge(edge_ev, on="edge_key", how="left")

    for feature_df in [
        build_fine_mapping(args.xlsx),
        build_st15_features(args.xlsx),
        build_st18_features(args.xlsx),
        build_network_features(args.xlsx),
        build_st24_features(args.xlsx),
    ]:
        if not feature_df.empty:
            edge_features = edge_features.merge(feature_df, on="edge_key", how="left")

    st14 = build_st14_features(args.xlsx)
    if not st14.empty:
        edge_features = edge_features.merge(st14, on="edge_key_rsid", how="left")

    st12 = build_st12_features(args.xlsx)
    st13 = build_st13_features(args.xlsx)
    snp_features = edges[["variant_key_norm", "rsid"]].drop_duplicates().copy()
    snp_features = snp_features.merge(variant_ev, on="variant_key_norm", how="left")
    snp_features = snp_features.merge(st12, on="variant_key_norm", how="left")
    snp_features = snp_features.merge(st13, on="variant_key_norm", how="left")

    if args.snp_tissue_priors.exists():
        tissue = pd.read_csv(args.snp_tissue_priors, dtype={"rsid": str})
        tissue_keep = [
            "protein_id",
            "variant_id_hg37",
            "has_gtex_eqtl_coloc",
            "gtex_coloc_tissue_count",
            "gtex_coloc_max_pph4",
            "gtex_coloc_mean_pph4",
            "gtex_direction_concordant_fraction",
            "gtex_tissue_specificity_score",
            "snp_tissue_reg_prior",
            "gtex_cancer_site_matches",
        ]
        tissue = tissue[[c for c in tissue_keep if c in tissue.columns]].drop_duplicates(["protein_id", "variant_id_hg37"])
        edge_features = edge_features.merge(tissue, on=["protein_id", "variant_id_hg37"], how="left")

    protein_features = build_protein_features(args.xlsx, args.hpa_priors)
    edge_features = edge_features.merge(protein_features, on="protein_id", how="left")

    numeric_cols = [
        col
        for col in edge_features.columns
        if col.startswith(NUMERIC_FEATURE_PREFIXES) and pd.api.types.is_numeric_dtype(edge_features[col])
    ]
    edge_features[numeric_cols] = edge_features[numeric_cols].fillna(0.0).astype(np.float32)
    snp_numeric_cols = [col for col in snp_features.columns if col.startswith(("variant_", "st12_", "st13_")) and pd.api.types.is_numeric_dtype(snp_features[col])]
    snp_features[snp_numeric_cols] = snp_features[snp_numeric_cols].fillna(0.0).astype(np.float32)

    edge_features.to_csv(args.outdir / "st10_edge_features.csv", index=False)
    snp_features.to_csv(args.outdir / "st10_snp_features.csv", index=False)
    protein_features.to_csv(args.outdir / "st10_protein_features.csv", index=False)
    with (args.outdir / "st10_edge_feature_names.txt").open("w") as out:
        out.write("\n".join(numeric_cols))
        out.write("\n")

    summary = args.outdir / "st10_feature_summary.txt"
    with summary.open("w") as out:
        out.write(f"edges: {len(edge_features)}\n")
        out.write(f"snps: {snp_features['variant_key_norm'].nunique()}\n")
        out.write(f"proteins: {protein_features['protein_id'].nunique()}\n")
        out.write(f"numeric_edge_features: {len(numeric_cols)}\n")
        out.write(f"edges_with_gtex_coloc: {int(edge_features.get('has_gtex_eqtl_coloc', pd.Series(0, index=edge_features.index)).fillna(0).sum())}\n")
        out.write(f"edges_with_finemap_top: {int(edge_features.get('st16_is_finemap_top', pd.Series(0, index=edge_features.index)).fillna(0).sum())}\n")
        out.write(f"edges_in_credible_set: {int(edge_features.get('st16_is_credible_set', pd.Series(0, index=edge_features.index)).fillna(0).sum())}\n")
        out.write(f"edges_with_regulatory_noncoding_variant: {int(edge_features.get('variant_ev_regulatory_non_coding_pqtl', pd.Series(0, index=edge_features.index)).fillna(0).sum())}\n")
        out.write(f"edges_with_high_impact_coding_variant: {int(edge_features.get('variant_ev_high_impact_coding_or_splicing_pqtl', pd.Series(0, index=edge_features.index)).fillna(0).sum())}\n")
        out.write(f"edges_with_network_evidence: {int((edge_features.get('st20_ppi_flag', 0) + edge_features.get('st21_reciprocal_trans_flag', 0) + edge_features.get('st23_receptor_ligand_flag', 0)).fillna(0).gt(0).sum())}\n")
    print(summary.read_text())


if __name__ == "__main__":
    main()
