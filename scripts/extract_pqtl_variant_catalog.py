#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


HEADER_ROWS = {
    "ST3": 2,
    "ST9": 4,
    "ST10": 4,
    "ST11": 5,
    "ST12": 2,
    "ST13": 2,
    "ST14": 2,
    "ST15": 3,
    "ST16": 3,
    "ST17": 3,
    "ST18": 3,
    "ST19": 2,
    "ST20": 2,
    "ST21": 4,
    "ST23": 3,
    "ST24": 3,
    "ST26": 3,
    "ST27": 4,
    "ST28": 3,
    "ST29": 3,
}


def clean_col(col: object) -> str:
    return re.sub(r"\s+", " ", str(col).strip())


def read_sheet(path: Path, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet, header=HEADER_ROWS[sheet], dtype=object)
    df.columns = [clean_col(c) for c in df.columns]
    return df.dropna(how="all").reset_index(drop=True)


def split_values(value: object) -> list[str]:
    if pd.isna(value):
        return []
    out = []
    for item in re.split(r"[,;]\s*", str(value)):
        item = item.strip()
        if item and item not in {"-", "nan"}:
            out.append(item)
    return out


def normalize_variant_id(value: object) -> str | None:
    if pd.isna(value):
        return None
    value = str(value).strip()
    if not value or value in {"-", "nan", "NaN"}:
        return None
    return value


def first_existing(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def add_rows(
    rows: list[dict],
    df: pd.DataFrame,
    sheet: str,
    evidence: str,
    variant_col: str,
    protein_col: str | None = None,
    rsid_col: str | None = None,
    cis_col: str | None = None,
    gene_col: str | None = None,
    consequence_col: str | None = None,
    impact_col: str | None = None,
    beta_col: str | None = None,
    logp_col: str | None = None,
    split_variant_col: bool = False,
) -> None:
    for _, r in df.iterrows():
        values = split_values(r.get(variant_col)) if split_variant_col else [normalize_variant_id(r.get(variant_col))]
        for variant in values:
            if not variant:
                continue
            rows.append(
                {
                    "variant_id": variant,
                    "rsid": normalize_variant_id(r.get(rsid_col)) if rsid_col else None,
                    "source_sheet": sheet,
                    "evidence_type": evidence,
                    "ukbppp_protein_id": normalize_variant_id(r.get(protein_col)) if protein_col else None,
                    "cis_trans": normalize_variant_id(r.get(cis_col)) if cis_col else None,
                    "gene_or_locus": normalize_variant_id(r.get(gene_col)) if gene_col else None,
                    "consequence": normalize_variant_id(r.get(consequence_col)) if consequence_col else None,
                    "impact": normalize_variant_id(r.get(impact_col)) if impact_col else None,
                    "beta": r.get(beta_col) if beta_col else None,
                    "log10p": r.get(logp_col) if logp_col else None,
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx", type=Path)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []

    for sheet in ["ST9", "ST10", "ST11", "ST15"]:
        df = read_sheet(args.xlsx, sheet)
        add_rows(
            rows,
            df,
            sheet,
            {
                "ST9": "discovery_significant_pqtl",
                "ST10": "combined_significant_pqtl",
                "ST11": "non_european_significant_pqtl",
                "ST15": "novel_or_replicated_discovery_pqtl",
            }[sheet],
            "Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)",
            protein_col="UKBPPP ProteinID",
            rsid_col="rsID",
            cis_col="cis/trans",
            gene_col=first_existing(df, ["Bioinfomatic annotated gene", "cis gene"]),
            consequence_col=first_existing(df, ["Annotated gene consequence"]),
            impact_col=first_existing(df, ["IMPACT"]),
            beta_col=first_existing(df, ["BETA", "BETA (discovery, wrt. A1)"]),
            logp_col=first_existing(df, ["log10(p)", "log10(p) (discovery)"]),
        )

    for sheet, evidence in [
        ("ST12", "high_impact_coding_or_splicing_pqtl"),
        ("ST13", "regulatory_non_coding_pqtl"),
    ]:
        df = read_sheet(args.xlsx, sheet)
        add_rows(
            rows,
            df,
            sheet,
            evidence,
            "Variant ID",
            gene_col=first_existing(df, ["Gene", "Ensembl ID"]),
            consequence_col=first_existing(df, ["Consequence", "BIOTYPE"]),
            impact_col=first_existing(df, ["IMPACT"]),
        )

    df16 = read_sheet(args.xlsx, "ST16")
    # First row is the second-level header for top variant / credible set columns.
    df16 = df16.iloc[1:].reset_index(drop=True)
    add_rows(
        rows,
        df16,
        "ST16",
        "fine_mapped_top_variant",
        "Top variant (defined by highest PIP)",
        protein_col="UKBPPP ProteinID",
        rsid_col="Unnamed: 2",
        cis_col="Unnamed: 9",
        beta_col="Unnamed: 5",
        logp_col="Unnamed: 7",
    )
    add_rows(
        rows,
        df16,
        "ST16",
        "fine_mapped_credible_set_variant",
        "Unnamed: 13",
        protein_col="UKBPPP ProteinID",
        cis_col="Unnamed: 9",
        split_variant_col=True,
    )

    for sheet, evidence, variant_col, protein_col in [
        ("ST17", "credible_set_contains_pav", "PAVs variant ID", "UKBPPP_ProteinID"),
        ("ST18", "pqtl_pqtl_colocalized_cis_top", "Top variant ID", "UKBPPP ProteinID"),
        ("ST18", "pqtl_pqtl_colocalized_other_top", "Top variant IDs", "UKBPPP ProteinID"),
        ("ST20", "trans_locus_interacting_protein", "ID", "UKBPPP ProteinID"),
        ("ST21", "reciprocal_trans_interaction", "ID", "UKBPPP ProteinID"),
        ("ST23", "receptor_ligand_trans_pqtl", "Variant ID", "Target protein (UKBPPP ProteinID)"),
        ("ST24", "pqtl_covariate_robustness", "Variant ID", "UKBPPP ProteinID"),
        ("ST26", "cis_pqtl_eqtl_colocalized_top", "Top variant ID", "UKBPPP ProteinID"),
        ("ST29", "inflammasome_trans_pqtl", "Sentinel Variant", "Protein"),
    ]:
        df = read_sheet(args.xlsx, sheet)
        add_rows(
            rows,
            df,
            sheet,
            evidence,
            variant_col,
            protein_col=protein_col,
            rsid_col=first_existing(df, ["rsID", "rsid"]),
            cis_col=first_existing(df, ["cis/trans", "Cis/Trans"]),
            gene_col=first_existing(df, ["Gene name", "Annotated/nearest gene", "Interacting protein at trans locus"]),
            beta_col=first_existing(df, ["BETA", "Beta"]),
            logp_col=first_existing(df, ["LOG10P", "log10(p)"]),
            split_variant_col=(sheet == "ST18" and variant_col == "Top variant IDs"),
        )

    df14 = read_sheet(args.xlsx, "ST14")
    for _, r in df14.iterrows():
        for col, evidence in [
            ("rsID (lead variant)", "lead_variant_in_ld_with_pav"),
            ("rsID (PAV)", "pav_in_ld_with_lead_variant"),
        ]:
            rsid = normalize_variant_id(r.get(col))
            if rsid:
                rows.append(
                    {
                        "variant_id": None,
                        "rsid": rsid,
                        "source_sheet": "ST14",
                        "evidence_type": evidence,
                        "ukbppp_protein_id": normalize_variant_id(r.get("UKBPPP ProteinID")),
                        "cis_trans": None,
                        "gene_or_locus": None,
                        "consequence": normalize_variant_id(r.get("PAV consequence")),
                        "impact": normalize_variant_id(r.get("PAV impact")),
                        "beta": None,
                        "log10p": None,
                    }
                )

    catalog = pd.DataFrame(rows)
    catalog["variant_key"] = catalog["variant_id"].fillna(catalog["rsid"])
    catalog = catalog.dropna(subset=["variant_key"]).drop_duplicates()
    catalog.to_csv(args.outdir / "pqtl_variant_catalog.csv", index=False)

    summary = (
        catalog.groupby(["source_sheet", "evidence_type"], dropna=False)
        .agg(n_rows=("variant_key", "size"), n_unique_variants=("variant_key", "nunique"))
        .reset_index()
        .sort_values(["source_sheet", "evidence_type"])
    )
    summary.to_csv(args.outdir / "pqtl_variant_catalog_summary.csv", index=False)

    sheet_catalog = []
    for sheet in pd.ExcelFile(args.xlsx).sheet_names:
        if sheet == "Contents":
            continue
        try:
            df = read_sheet(args.xlsx, sheet) if sheet in HEADER_ROWS else pd.read_excel(args.xlsx, sheet_name=sheet, dtype=object)
            sheet_catalog.append({"sheet": sheet, "rows": len(df), "cols": len(df.columns), "columns": "; ".join(map(clean_col, df.columns))})
        except Exception as exc:
            sheet_catalog.append({"sheet": sheet, "rows": None, "cols": None, "columns": f"ERROR: {exc}"})
    pd.DataFrame(sheet_catalog).to_csv(args.outdir / "supplement_sheet_catalog.csv", index=False)


if __name__ == "__main__":
    main()
