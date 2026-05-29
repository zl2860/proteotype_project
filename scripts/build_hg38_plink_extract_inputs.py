#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def clean_token(x: object) -> str | None:
    if pd.isna(x):
        return None
    x = str(x).strip()
    if not x or x in {"-", "nan", "NaN"}:
        return None
    return x


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, default=Path("文献资料补充/Nature-PQTL-补充材料.xlsx"))
    parser.add_argument("--sheet", default="ST10")
    parser.add_argument("--outdir", type=Path, default=Path("plink_variant_lists/hg38_st10_core"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(args.xlsx, sheet_name=args.sheet, header=4, dtype=object)
    needed = [
        "Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)",
        "CHROM",
        "GENPOS (hg38)",
        "rsID",
        "A1FREQ",
        "BETA",
        "SE",
        "log10(p)",
        "cis/trans",
        "UKBPPP ProteinID",
        "Assay Target",
        "Target UniProt",
    ]
    df = df[needed].copy()
    df.columns = [
        "variant_id_hg37",
        "chr",
        "pos_hg38",
        "rsid",
        "a1freq",
        "beta",
        "se",
        "log10p",
        "cis_trans",
        "ukbppp_protein_id",
        "assay_target",
        "target_uniprot",
    ]
    df = df.dropna(subset=["chr", "pos_hg38"]).copy()
    df["chr"] = df["chr"].astype(str).str.replace("^chr", "", regex=True)
    df["pos_hg38"] = df["pos_hg38"].astype("int64")
    df["rsid"] = df["rsid"].map(clean_token)
    df["coord_hg38"] = df["chr"].astype(str) + ":" + df["pos_hg38"].astype(str)

    coord = (
        df.sort_values(["chr", "pos_hg38", "rsid"], key=lambda s: s.astype(str))
        .drop_duplicates(["chr", "pos_hg38"])
        .copy()
    )
    coord.to_csv(args.outdir / "st10_core_hg38_coordinates_with_rsid.tsv", sep="\t", index=False)

    coord[["chr", "pos_hg38", "rsid"]].to_csv(
        args.outdir / "st10_core_hg38_chr_pos_rsid.tsv", sep="\t", index=False
    )
    coord[["coord_hg38"]].to_csv(
        args.outdir / "st10_core_hg38_chrpos.txt", sep="\t", index=False, header=False
    )
    coord["rsid"].dropna().drop_duplicates().sort_values().to_csv(
        args.outdir / "st10_core_rsid_only.txt", index=False, header=False
    )

    for chrom, sub in coord.groupby("chr", sort=False):
        chrom_dir = args.outdir / "by_chr"
        chrom_dir.mkdir(exist_ok=True)
        sub[["pos_hg38", "rsid"]].drop_duplicates().sort_values("pos_hg38").to_csv(
            chrom_dir / f"chr{chrom}.hg38_pos_rsid.tsv", sep="\t", index=False, header=False
        )

    readme = args.outdir / "README.md"
    readme.write_text(
        "# ST10 pQTL hg38 extraction inputs\n\n"
        "These files are built from `ST10. List of significant p<1.7e-11 pQTLs in combined cohort`.\n\n"
        "- `st10_core_hg38_coordinates_with_rsid.tsv`: full annotation table with hg38 position and rsID.\n"
        "- `st10_core_hg38_chr_pos_rsid.tsv`: compact `chr pos_hg38 rsid` table.\n"
        "- `st10_core_hg38_chrpos.txt`: one `chr:pos` hg38 coordinate per line.\n"
        "- `st10_core_rsid_only.txt`: rsID-only list.\n"
        "- `by_chr/chr*.hg38_pos_rsid.tsv`: per-chromosome `pos_hg38 rsid` files for matching against PLINK `.bim`.\n\n"
        "For your server `.bed/.bim/.fam` files, first match these hg38 positions to the actual variant IDs in `.bim`, then run PLINK/PLINK2 with `--extract` on those matched IDs.\n"
    )


if __name__ == "__main__":
    main()
