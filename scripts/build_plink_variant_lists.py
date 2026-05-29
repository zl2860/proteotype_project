#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


TIER_DEFINITIONS = {
    "tier00_core_st10_all": {
        "description": "Combined-cohort significant pQTL lead variants from ST10. Recommended first extraction.",
        "evidence": {"combined_significant_pqtl"},
    },
    "tier01_core_st10_cis": {
        "description": "ST10 cis pQTL lead variants only.",
        "evidence": {"combined_significant_pqtl"},
        "cis_trans": {"cis"},
    },
    "tier02_core_st10_trans": {
        "description": "ST10 trans pQTL lead variants only.",
        "evidence": {"combined_significant_pqtl"},
        "cis_trans": {"trans"},
    },
    "tier03_core_plus_finemap_top": {
        "description": "ST10 lead variants plus ST16 top-PIP fine-mapped variants.",
        "evidence": {"combined_significant_pqtl", "fine_mapped_top_variant"},
    },
    "tier04_functional_priority": {
        "description": "High-priority functional variants: coding/splicing, regulatory, PAV-linked, and PAV-containing credible sets.",
        "evidence": {
            "high_impact_coding_or_splicing_pqtl",
            "regulatory_non_coding_pqtl",
            "lead_variant_in_ld_with_pav",
            "pav_in_ld_with_lead_variant",
            "credible_set_contains_pav",
        },
    },
    "tier05_eqtl_coloc_cis": {
        "description": "GTEx eQTL-colocalized cis pQTL variants from ST26.",
        "evidence": {"cis_pqtl_eqtl_colocalized_top"},
    },
    "tier06_trans_network": {
        "description": "Trans-network variants: pQTL-pQTL coloc, interacting proteins, reciprocal trans, receptor-ligand, inflammasome.",
        "evidence": {
            "pqtl_pqtl_colocalized_cis_top",
            "pqtl_pqtl_colocalized_other_top",
            "trans_locus_interacting_protein",
            "reciprocal_trans_interaction",
            "receptor_ligand_trans_pqtl",
            "inflammasome_trans_pqtl",
        },
    },
    "tier07_discovery_and_non_eur": {
        "description": "Discovery significant pQTLs plus non-European ancestry pQTLs.",
        "evidence": {"discovery_significant_pqtl", "non_european_significant_pqtl"},
    },
    "tier08_full_finemap_credible_sets": {
        "description": "Full ST16 credible-set variants. Large sensitivity-analysis extraction.",
        "evidence": {"fine_mapped_credible_set_variant"},
    },
    "tier09_all_except_full_credible_sets": {
        "description": "All extracted catalog variants except full credible-set expansion.",
        "exclude_evidence": {"fine_mapped_credible_set_variant"},
    },
    "tier10_all_catalog": {
        "description": "All extracted variant keys, including full credible-set expansion.",
    },
}


def parse_variant_id(variant_id: object) -> dict[str, object]:
    if pd.isna(variant_id):
        return {}
    value = str(variant_id).strip()
    # Expected: CHROM:GENPOS(hg37):A0:A1:imp:v1
    m = re.match(r"^([^:]+):([0-9]+):([^:]+):([^:]+)(?::.*)?$", value)
    if not m:
        return {}
    chrom, pos, a0, a1 = m.groups()
    return {
        "chrom": chrom,
        "pos_hg37_from_variant_id": int(pos),
        "a0": a0,
        "a1": a1,
        "variant_id_no_suffix": f"{chrom}:{pos}:{a0}:{a1}",
    }


def variant_subset(catalog: pd.DataFrame, definition: dict) -> pd.DataFrame:
    out = catalog.copy()
    if "evidence" in definition:
        out = out[out["evidence_type"].isin(definition["evidence"])]
    if "exclude_evidence" in definition:
        out = out[~out["evidence_type"].isin(definition["exclude_evidence"])]
    if "cis_trans" in definition:
        out = out[out["cis_trans"].isin(definition["cis_trans"])]
    return out


def write_unique(series: pd.Series, path: Path) -> int:
    values = sorted(
        {
            str(x).strip()
            for x in series.dropna().tolist()
            if str(x).strip() and str(x).strip().lower() != "nan"
        }
    )
    path.write_text("\n".join(values) + ("\n" if values else ""))
    return len(values)


def write_coordinates(df: pd.DataFrame, path: Path, extra_cols: list[str] | None = None) -> None:
    cols = ["variant_id", "rsid", "chrom", "pos_hg37_from_variant_id", "a0", "a1", "variant_id_no_suffix"]
    if extra_cols:
        cols.extend(extra_cols)
    coord = df.dropna(subset=["chrom", "pos_hg37_from_variant_id", "a0", "a1"])[cols].drop_duplicates().copy()
    coord["pos_hg37_from_variant_id"] = coord["pos_hg37_from_variant_id"].astype("int64")
    coord.to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=Path("文献资料补充/pqtl_variant_catalog/pqtl_variant_catalog.csv"))
    parser.add_argument("--xlsx", type=Path, default=Path("文献资料补充/Nature-PQTL-补充材料.xlsx"))
    parser.add_argument("--outdir", type=Path, default=Path("plink_variant_lists"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(args.catalog, low_memory=False)

    parsed = catalog["variant_id"].apply(parse_variant_id).apply(pd.Series)
    catalog = pd.concat([catalog, parsed], axis=1)

    if args.xlsx.exists():
        st3 = pd.read_excel(args.xlsx, sheet_name="ST3", header=2, dtype=object)
        st3 = st3[["UKBPPP ProteinID", "Protein panel", "Gene symbol", "UniProt"]].drop_duplicates()
        catalog = catalog.merge(st3, left_on="ukbppp_protein_id", right_on="UKBPPP ProteinID", how="left")
    else:
        catalog["Protein panel"] = pd.NA

    catalog["is_oncology_panel"] = catalog["Protein panel"].isin(["Oncology", "Oncology_II"])

    # Per-variant manifest with aggregated evidence and tier flags.
    group_cols = ["variant_key"]
    agg = catalog.groupby(group_cols, dropna=False).agg(
        variant_id=("variant_id", lambda x: next((v for v in x.dropna().astype(str) if v), None)),
        rsid=("rsid", lambda x: next((v for v in x.dropna().astype(str) if v), None)),
        chrom=("chrom", lambda x: next((v for v in x.dropna().astype(str) if v), None)),
        pos_hg37_from_variant_id=("pos_hg37_from_variant_id", lambda x: next((v for v in x.dropna().tolist() if pd.notna(v)), None)),
        a0=("a0", lambda x: next((v for v in x.dropna().astype(str) if v), None)),
        a1=("a1", lambda x: next((v for v in x.dropna().astype(str) if v), None)),
        variant_id_no_suffix=("variant_id_no_suffix", lambda x: next((v for v in x.dropna().astype(str) if v), None)),
        evidence_types=("evidence_type", lambda x: ";".join(sorted(set(x.dropna().astype(str))))),
        source_sheets=("source_sheet", lambda x: ";".join(sorted(set(x.dropna().astype(str))))),
        cis_trans_values=("cis_trans", lambda x: ";".join(sorted(set(x.dropna().astype(str))))),
        proteins=("ukbppp_protein_id", lambda x: ";".join(sorted(set(x.dropna().astype(str)))[:20])),
        genes=("gene_or_locus", lambda x: ";".join(sorted(set(x.dropna().astype(str)))[:20])),
        protein_panels=("Protein panel", lambda x: ";".join(sorted(set(x.dropna().astype(str))))),
        is_oncology_panel=("is_oncology_panel", "max"),
    ).reset_index()

    tier_rows = []
    for tier, definition in TIER_DEFINITIONS.items():
        sub = variant_subset(catalog, definition)
        ids = set(sub["variant_key"].dropna().astype(str))
        agg[tier] = agg["variant_key"].astype(str).isin(ids)

        tier_dir = args.outdir / tier
        tier_dir.mkdir(exist_ok=True)
        n_variant_id = write_unique(sub["variant_id"], tier_dir / "extract.variant_id.txt")
        n_variant_no_suffix = write_unique(sub["variant_id_no_suffix"], tier_dir / "extract.variant_id_no_suffix.txt")
        n_rsid = write_unique(sub["rsid"], tier_dir / "extract.rsid.txt")
        n_any = write_unique(sub["variant_key"], tier_dir / "extract.variant_key.txt")

        write_coordinates(sub, tier_dir / "variants.coordinates_hg37.tsv")

        tier_rows.append(
            {
                "tier": tier,
                "description": definition["description"],
                "n_variant_key": n_any,
                "n_variant_id": n_variant_id,
                "n_variant_id_no_suffix": n_variant_no_suffix,
                "n_rsid": n_rsid,
                "directory": str(tier_dir),
            }
        )

    oncology = catalog[catalog["is_oncology_panel"] & catalog["evidence_type"].isin({"combined_significant_pqtl", "fine_mapped_top_variant"})]
    oncology_dir = args.outdir / "tier11_oncology_panel_st10_plus_finemap_top"
    oncology_dir.mkdir(exist_ok=True)
    n_variant_id = write_unique(oncology["variant_id"], oncology_dir / "extract.variant_id.txt")
    n_variant_no_suffix = write_unique(oncology["variant_id_no_suffix"], oncology_dir / "extract.variant_id_no_suffix.txt")
    n_rsid = write_unique(oncology["rsid"], oncology_dir / "extract.rsid.txt")
    n_any = write_unique(oncology["variant_key"], oncology_dir / "extract.variant_key.txt")
    write_coordinates(oncology, oncology_dir / "variants.coordinates_hg37.tsv", extra_cols=["ukbppp_protein_id", "Protein panel"])
    agg["tier11_oncology_panel_st10_plus_finemap_top"] = agg["variant_key"].astype(str).isin(set(oncology["variant_key"].dropna().astype(str)))
    tier_rows.append(
        {
            "tier": "tier11_oncology_panel_st10_plus_finemap_top",
            "description": "Oncology/Oncology_II Olink panel proteins, restricted to ST10 lead variants plus ST16 top-PIP variants.",
            "n_variant_key": n_any,
            "n_variant_id": n_variant_id,
            "n_variant_id_no_suffix": n_variant_no_suffix,
            "n_rsid": n_rsid,
            "directory": str(oncology_dir),
        }
    )

    agg.to_csv(args.outdir / "plink_variant_manifest.csv", index=False)
    pd.DataFrame(tier_rows).to_csv(args.outdir / "plink_variant_list_summary.csv", index=False)

    readme = args.outdir / "README.md"
    readme.write_text(
        "# PLINK2 pQTL variant extraction lists\n\n"
        "Each tier directory contains:\n\n"
        "- `extract.variant_id.txt`: IDs exactly as reported in the pQTL supplement, usually `CHROM:POS_hg37:A0:A1:imp:v1`.\n"
        "- `extract.variant_id_no_suffix.txt`: same as above but without `:imp:v1`, useful if the PLINK `.pvar` IDs omit the imputation suffix.\n"
        "- `extract.rsid.txt`: rsID list where available.\n"
        "- `extract.variant_key.txt`: variant ID when available, otherwise rsID. Use mainly as a generic reference list.\n"
        "- `variants.coordinates_hg37.tsv`: parsed chromosome, hg37 position, A0, A1, and rsID for server-side ID reconciliation.\n\n"
        "Recommended first extraction: `tier00_core_st10_all/extract.variant_id.txt`.\n"
        "If PLINK reports zero or low matches, retry with `extract.variant_id_no_suffix.txt` or `extract.rsid.txt`, after checking `.pvar`/`.bim` ID format.\n\n"
        "Example commands:\n\n"
        "```sh\n"
        "plink2 --pfile UKB_IMPUTED_PREFIX --extract plink_variant_lists/tier00_core_st10_all/extract.variant_id.txt --make-pgen --out ukb_pqtl_tier00\n"
        "plink2 --pfile UKB_IMPUTED_PREFIX --extract plink_variant_lists/tier00_core_st10_all/extract.rsid.txt --make-pgen --out ukb_pqtl_tier00_rsid\n"
        "plink2 --pfile UKB_IMPUTED_PREFIX --extract plink_variant_lists/tier00_core_st10_all/extract.variant_id.txt --export A --out ukb_pqtl_tier00_dosage\n"
        "```\n\n"
        "Coordinate caution: the position embedded in `variant_id` is hg37. The supplement also contains hg38 columns, but these extraction lists intentionally follow the hg37 ID because UKBB imputed genotype resources are commonly GRCh37/hg19 indexed.\n"
    )


if __name__ == "__main__":
    main()
