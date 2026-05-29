# ST10 pQTL PLINK extraction inputs

These files are built from `ST10. List of significant p<1.7e-11 pQTLs in
combined cohort` in the Nature pQTL supplementary workbook.

## Use rsID for the current UKB server

The server `.bim` examples show rsIDs in column 2:

```text
1  rs575272151  0  11008  G  C
```

So the recommended extraction key is:

- `st10_core_rsid_only.txt`: one rsID per line, 11,623 rsIDs.

Run on the UKB server after copying this whole `hg38_st10_core` folder:

```sh
cd /path/to/hg38_st10_core
bash extract_ukb_st10_core_by_rsid_on_server.sh
```

The script reads:

- UKB genotype files from `/data/public/UK_BIOBANK/processed_data/imputation_4_QC_UKBfull`
- per-chromosome files named `UKB.QC.chr*.bed/.bim/.fam`
- rsIDs from `st10_core_rsid_only.txt`

It writes to:

```text
/data/public/UK_BIOBANK/processed_data/proteomic_traj_snp/st10_core_rsid/
```

Main outputs:

- `matched_ids/chr*.matched_rsid.txt`: matched rsIDs per chromosome.
- `st10_core_rsid.matched_all.txt`: all matched rsIDs.
- `st10_core_rsid.matched_all_bim_map.tsv`: chromosome, rsID, BIM position,
  allele1, allele2.
- `bed_by_chr/st10_core_rsid.chr*.*`: extracted PLINK bed files.
- `raw_by_chr/st10_core_rsid.chr*.raw`: additive hard-call matrices.
- `st10_core_rsid.run_summary.txt`: input/matched counts.

## Override paths

```sh
export UKB_DIR=/data/public/UK_BIOBANK/processed_data/imputation_4_QC_UKBfull
export OUT_DIR=/data/public/UK_BIOBANK/processed_data/proteomic_traj_snp
export RSID_LIST=/path/to/hg38_st10_core/st10_core_rsid_only.txt
bash /path/to/hg38_st10_core/extract_ukb_st10_core_by_rsid_on_server.sh
```

## Coordinate files

The hg38 coordinate files are kept for annotation and troubleshooting, but they
are not the default extraction key for this UKB server:

- `st10_core_hg38_coordinates_with_rsid.tsv`
- `st10_core_hg38_chr_pos_rsid.tsv`
- `st10_core_hg38_chrpos.txt`
- `by_chr/chr*.hg38_pos_rsid.tsv`

If rsID matching is unexpectedly low, inspect the `.bim` column 2 IDs and compare
against `st10_core_rsid_only.txt`.
