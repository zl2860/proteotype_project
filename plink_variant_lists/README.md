# PLINK2 pQTL variant extraction lists

Each tier directory contains:

- `extract.variant_id.txt`: IDs exactly as reported in the pQTL supplement, usually `CHROM:POS_hg37:A0:A1:imp:v1`.
- `extract.variant_id_no_suffix.txt`: same as above but without `:imp:v1`, useful if the PLINK `.pvar` IDs omit the imputation suffix.
- `extract.rsid.txt`: rsID list where available.
- `extract.variant_key.txt`: variant ID when available, otherwise rsID. Use mainly as a generic reference list.
- `variants.coordinates_hg37.tsv`: parsed chromosome, hg37 position, A0, A1, and rsID for server-side ID reconciliation.

Recommended first extraction: `tier00_core_st10_all/extract.variant_id.txt`.
If PLINK reports zero or low matches, retry with `extract.variant_id_no_suffix.txt` or `extract.rsid.txt`, after checking `.pvar`/`.bim` ID format.

Example commands:

```sh
plink2 --pfile UKB_IMPUTED_PREFIX --extract plink_variant_lists/tier00_core_st10_all/extract.variant_id.txt --make-pgen --out ukb_pqtl_tier00
plink2 --pfile UKB_IMPUTED_PREFIX --extract plink_variant_lists/tier00_core_st10_all/extract.rsid.txt --make-pgen --out ukb_pqtl_tier00_rsid
plink2 --pfile UKB_IMPUTED_PREFIX --extract plink_variant_lists/tier00_core_st10_all/extract.variant_id.txt --export A --out ukb_pqtl_tier00_dosage
```

Coordinate caution: the position embedded in `variant_id` is hg37. The supplement also contains hg38 columns, but these extraction lists intentionally follow the hg37 ID because UKBB imputed genotype resources are commonly GRCh37/hg19 indexed.
