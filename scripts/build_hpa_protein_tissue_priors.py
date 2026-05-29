#!/usr/bin/env python3

from __future__ import annotations

import argparse
import io
import math
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd


HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"

SPECIFICITY_WEIGHTS = {
    "Tissue enriched": 1.00,
    "Group enriched": 0.85,
    "Tissue enhanced": 0.70,
    "Low tissue specificity": 0.25,
    "Not detected": 0.00,
}

DISTRIBUTION_WEIGHTS = {
    "Detected in single": 1.00,
    "Detected in some": 0.75,
    "Detected in many": 0.45,
    "Detected in all": 0.20,
    "Not detected": 0.00,
}


def download_hpa_tsv(cache_path: Path, force: bool) -> Path:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists() and not force:
        return cache_path
    with urllib.request.urlopen(HPA_URL, timeout=120) as response:
        data = response.read()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = zf.namelist()
        if not names:
            raise RuntimeError("HPA archive did not contain a TSV file.")
        with zf.open(names[0]) as src, cache_path.open("wb") as dst:
            dst.write(src.read())
    return cache_path


def parse_max_expression(value: object) -> tuple[str, float]:
    if pd.isna(value) or not str(value).strip():
        return "", 0.0
    best_label = ""
    best_value = 0.0
    for part in str(value).split(";"):
        if ":" not in part:
            continue
        label, raw = part.rsplit(":", 1)
        try:
            numeric = float(raw.strip())
        except ValueError:
            continue
        if numeric > best_value:
            best_label = label.strip()
            best_value = numeric
    return best_label, best_value


def normalized_score(category: object, distribution: object, score: object, expression: object) -> float:
    category_weight = SPECIFICITY_WEIGHTS.get(str(category), 0.15)
    distribution_weight = DISTRIBUTION_WEIGHTS.get(str(distribution), 0.15)
    try:
        raw_score = float(score)
    except (TypeError, ValueError):
        score_weight = category_weight
    else:
        if raw_score <= 1.0:
            score_weight = raw_score
        else:
            score_weight = min(1.0, math.log1p(raw_score) / math.log1p(100.0))
    _, max_expr = parse_max_expression(expression)
    expr_weight = min(1.0, math.log1p(max_expr) / math.log1p(1000.0)) if max_expr > 0 else 0.0
    return max(
        0.0,
        min(1.0, 0.45 * category_weight + 0.35 * score_weight + 0.20 * max(distribution_weight, expr_weight)),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-manifest", type=Path, default=Path("model_data/st10_core/st10_protein_manifest.csv"))
    parser.add_argument("--hpa-cache", type=Path, default=Path("model_data/reference/hpa_proteinatlas.tsv"))
    parser.add_argument("--out", type=Path, default=Path("model_data/st10_core/st10_protein_tissue_priors.csv"))
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.protein_manifest)
    hpa_path = download_hpa_tsv(args.hpa_cache, args.force_download)
    hpa = pd.read_csv(
        hpa_path,
        sep="\t",
        usecols=[
            "Gene",
            "Uniprot",
            "RNA tissue specificity",
            "RNA tissue distribution",
            "RNA tissue specificity score",
            "RNA tissue specific nTPM",
            "Protein tissue specificity",
            "Protein tissue distribution",
            "Protein tissue specificity score",
            "Protein tissue specific Intensity",
            "Secretome location",
            "Protein class",
        ],
        dtype=str,
    )
    hpa["uniprot_key"] = hpa["Uniprot"].fillna("").str.split(",").str[0].str.strip()
    hpa_by_uniprot = hpa[hpa["uniprot_key"].ne("")].drop_duplicates("uniprot_key", keep="first")
    hpa_by_gene = hpa.drop_duplicates("Gene", keep="first")
    manifest["uniprot_key"] = manifest["uniprot"].fillna("").astype(str).str.split(",").str[0].str.strip()
    merged = manifest.merge(hpa_by_uniprot, on="uniprot_key", how="left", suffixes=("", "_hpa"))
    missing = merged["Gene"].isna()
    if missing.any():
        by_gene = manifest.loc[missing.to_numpy()].merge(
            hpa_by_gene,
            left_on="gene_symbol",
            right_on="Gene",
            how="left",
            suffixes=("", "_hpa"),
        )
        merged.loc[missing, hpa.columns] = by_gene[hpa.columns].to_numpy()

    rna_best = merged["RNA tissue specific nTPM"].map(parse_max_expression)
    protein_best = merged["Protein tissue specific Intensity"].map(parse_max_expression)
    merged["hpa_top_rna_tissue"] = [x[0] for x in rna_best]
    merged["hpa_top_rna_ntpm"] = [x[1] for x in rna_best]
    merged["hpa_top_protein_tissue"] = [x[0] for x in protein_best]
    merged["hpa_top_protein_intensity"] = [x[1] for x in protein_best]
    merged["rna_tissue_prior"] = [
        normalized_score(cat, dist, score, expr)
        for cat, dist, score, expr in zip(
            merged["RNA tissue specificity"],
            merged["RNA tissue distribution"],
            merged["RNA tissue specificity score"],
            merged["RNA tissue specific nTPM"],
        )
    ]
    merged["protein_tissue_prior"] = [
        normalized_score(cat, dist, score, expr)
        for cat, dist, score, expr in zip(
            merged["Protein tissue specificity"],
            merged["Protein tissue distribution"],
            merged["Protein tissue specificity score"],
            merged["Protein tissue specific Intensity"],
        )
    ]
    merged["hpa_tissue_prior"] = merged[["rna_tissue_prior", "protein_tissue_prior"]].max(axis=1)
    merged["hpa_has_match"] = merged["Gene"].notna().astype(int)

    keep = [
        "protein_id",
        "assay_target",
        "gene_symbol",
        "uniprot",
        "hpa_has_match",
        "RNA tissue specificity",
        "RNA tissue distribution",
        "RNA tissue specificity score",
        "hpa_top_rna_tissue",
        "hpa_top_rna_ntpm",
        "Protein tissue specificity",
        "Protein tissue distribution",
        "Protein tissue specificity score",
        "hpa_top_protein_tissue",
        "hpa_top_protein_intensity",
        "Secretome location",
        "Protein class",
        "rna_tissue_prior",
        "protein_tissue_prior",
        "hpa_tissue_prior",
    ]
    merged[keep].to_csv(args.out, index=False)

    summary = args.out.with_suffix(".summary.txt")
    with summary.open("w") as out:
        out.write(f"proteins: {len(merged)}\n")
        out.write(f"hpa_matches: {int(merged['hpa_has_match'].sum())}\n")
        out.write(f"mean_hpa_tissue_prior: {merged['hpa_tissue_prior'].mean():.4f}\n")
        out.write(f"median_hpa_tissue_prior: {merged['hpa_tissue_prior'].median():.4f}\n")
        out.write("\nRNA tissue specificity:\n")
        out.write(merged["RNA tissue specificity"].value_counts(dropna=False).to_string())
        out.write("\n\nProtein tissue specificity:\n")
        out.write(merged["Protein tissue specificity"].value_counts(dropna=False).to_string())
        out.write("\n")
    print(summary.read_text())


if __name__ == "__main__":
    main()
