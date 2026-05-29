#!/usr/bin/env bash
set -euo pipefail

# Run this on the UKB server.
#
# Recommended:
#   cd /path/to/copied/hg38_st10_core
#   bash extract_ukb_st10_core_by_rsid_on_server.sh
#
# The UKB .bim files shown by the user have rsIDs in column 2, so this script
# matches by rsID instead of coordinate. The hg38 coordinate files in this folder
# are kept as annotation/reference files, not as the extraction key.

UKB_DIR="${UKB_DIR:-/data/public/UK_BIOBANK/processed_data/imputation_4_QC_UKBfull}"
OUT_DIR="${OUT_DIR:-/data/public/UK_BIOBANK/processed_data/proteomic_traj_snp}"
LIST_DIR="${LIST_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RSID_LIST="${RSID_LIST:-${LIST_DIR}/st10_core_rsid_only.txt}"

PLINK="${PLINK:-${UKB_DIR}/plink}"
PLINK2="${PLINK2:-${UKB_DIR}/plink2}"

BED_DIR="${OUT_DIR}/st10_core_rsid/bed_by_chr"
RAW_DIR="${OUT_DIR}/st10_core_rsid/raw_by_chr"
MATCH_DIR="${OUT_DIR}/st10_core_rsid/matched_ids"
LOG_DIR="${OUT_DIR}/st10_core_rsid/logs"

mkdir -p "${BED_DIR}" "${RAW_DIR}" "${MATCH_DIR}" "${LOG_DIR}"

if [[ ! -s "${RSID_LIST}" ]]; then
  echo "ERROR: rsID list not found or empty: ${RSID_LIST}" >&2
  echo "Set RSID_LIST=/path/to/st10_core_rsid_only.txt or run from the hg38_st10_core folder." >&2
  exit 2
fi

if [[ ! -x "${PLINK}" ]]; then
  echo "ERROR: PLINK executable not found: ${PLINK}" >&2
  exit 2
fi

if [[ ! -x "${PLINK2}" ]]; then
  PLINK2=""
fi

echo "UKB_DIR=${UKB_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "RSID_LIST=${RSID_LIST}"
echo "PLINK=${PLINK}"
if [[ -n "${PLINK2}" ]]; then
  echo "PLINK2=${PLINK2}"
else
  echo "PLINK2 not found; using PLINK --recode A for raw export."
fi
echo

echo "[1/5] Matching rsIDs against UKB.QC.chr*.bim column 2"
for chr in {1..22}; do
  bim="${UKB_DIR}/UKB.QC.chr${chr}.bim"
  matched="${MATCH_DIR}/chr${chr}.matched_rsid.txt"
  map="${MATCH_DIR}/chr${chr}.matched_rsid_to_bim.tsv"

  if [[ ! -s "${bim}" ]]; then
    echo "ERROR: missing BIM file: ${bim}" >&2
    exit 2
  fi

  awk 'BEGIN{OFS="\t"}
       NR==FNR{rs[$1]=1; next}
       ($2 in rs){print $2 > ids; print $1,$2,$4,$5,$6 > map}' \
       ids="${matched}" map="${map}" "${RSID_LIST}" "${bim}"

  sort -u "${matched}" -o "${matched}"
  echo "chr${chr}: matched $(wc -l < "${matched}") rsIDs"
done

cat "${MATCH_DIR}"/chr*.matched_rsid.txt | sort -u > "${OUT_DIR}/st10_core_rsid.matched_all.txt"
cat "${MATCH_DIR}"/chr*.matched_rsid_to_bim.tsv > "${OUT_DIR}/st10_core_rsid.matched_all_bim_map.tsv"

echo
echo "[2/5] Overall match count"
wc -l "${RSID_LIST}" "${OUT_DIR}/st10_core_rsid.matched_all.txt"

echo
echo "[3/5] Extracting matched rsIDs into per-chromosome PLINK bed files"
for chr in {1..22}; do
  matched="${MATCH_DIR}/chr${chr}.matched_rsid.txt"
  [[ -s "${matched}" ]] || continue

  "${PLINK}" \
    --bfile "${UKB_DIR}/UKB.QC.chr${chr}" \
    --extract "${matched}" \
    --make-bed \
    --out "${BED_DIR}/st10_core_rsid.chr${chr}" \
    > "${LOG_DIR}/st10_core_rsid.chr${chr}.makebed.log" 2>&1
done

echo
echo "[4/5] Exporting additive hard-call matrices per chromosome"
for chr in {1..22}; do
  prefix="${BED_DIR}/st10_core_rsid.chr${chr}"
  [[ -s "${prefix}.bed" ]] || continue

  if [[ -n "${PLINK2}" ]]; then
    "${PLINK2}" \
      --bfile "${prefix}" \
      --export A \
      --out "${RAW_DIR}/st10_core_rsid.chr${chr}" \
      > "${LOG_DIR}/st10_core_rsid.chr${chr}.exportA.log" 2>&1
  else
    "${PLINK}" \
      --bfile "${prefix}" \
      --recode A \
      --out "${RAW_DIR}/st10_core_rsid.chr${chr}" \
      > "${LOG_DIR}/st10_core_rsid.chr${chr}.recodeA.log" 2>&1
  fi
done

echo
echo "[5/5] Writing run summary"
{
  echo "RSID_LIST=${RSID_LIST}"
  echo "UKB_DIR=${UKB_DIR}"
  echo "OUT_DIR=${OUT_DIR}"
  echo
  echo "Input rsIDs and matched rsIDs:"
  wc -l "${RSID_LIST}" "${OUT_DIR}/st10_core_rsid.matched_all.txt"
  echo
  echo "Per-chromosome matched rsIDs:"
  for chr in {1..22}; do
    printf "chr%s\t%s\n" "${chr}" "$(wc -l < "${MATCH_DIR}/chr${chr}.matched_rsid.txt")"
  done
} > "${OUT_DIR}/st10_core_rsid.run_summary.txt"

echo "Done."
echo "Summary: ${OUT_DIR}/st10_core_rsid.run_summary.txt"
echo "Matched rsIDs: ${OUT_DIR}/st10_core_rsid.matched_all.txt"
echo "BIM map: ${OUT_DIR}/st10_core_rsid.matched_all_bim_map.tsv"
echo "Per-chromosome bed: ${BED_DIR}/"
echo "Per-chromosome additive matrices: ${RAW_DIR}/"
