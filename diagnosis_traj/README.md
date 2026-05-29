# UKBB diagnosis trajectory preprocessing

This folder contains cleaned diagnosis sequence data derived from UKBB inpatient
and cancer registry diagnosis records.

## Source file

- `diagnosis_long_for_traj.qs2`: long diagnosis table with `eid`, `source`,
  `array_index`, `diagnosis_icd`, and `diagnosis_date`.

## Cleaning command

Run from the project root:

```sh
Rscript scripts/clean_diagnosis_sequence.R diagnosis_traj/diagnosis_long_for_traj.qs2 diagnosis_traj/processed
```

The script requires `data.table` and `qs2`. If `arrow` is installed, it also
writes Parquet copies.

## Main outputs

- `processed/diagnosis_events_clean.qs2`: deduplicated event-level table.
- `processed/person_diagnosis_sequences.qs2`: one row per participant with
  whitespace-separated token IDs for transformer input.
- `processed/diagnosis_token_vocab.csv`: token ID lookup table.
- `processed/diagnosis_cleaning_summary.txt`: cohort and event summaries.

## Cancer trajectory labels

The script keeps all diagnoses, and adds oncology-oriented labels:

- `primary_cancer`: ICD-10 `C00-C76`
- `metastasis`: ICD-10 `C77-C79`
- `cancer_unspecified`: ICD-10 `C80`
- `haematological_cancer`: ICD-10 `C81-C96`
- `in_situ_neoplasm`: ICD-10 `D00-D09`
- `benign_neoplasm`: ICD-10 `D10-D36`
- `uncertain_neoplasm`: ICD-10 `D37-D48`
- `non_neoplasm`: all other diagnosis codes

`sequence_phase` is derived relative to each participant's first malignant
record and is intended for pre-cancer, index cancer, post-cancer, and metastasis
trajectory modeling.
