#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "diagnosis_traj/diagnosis_long_for_traj.qs2"
output_dir <- if (length(args) >= 2) args[[2]] else "diagnosis_traj/processed"

if (!requireNamespace("qs2", quietly = TRUE)) {
  stop("Package 'qs2' is required. Install with: install.packages('qs2')")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading: ", input_path)
dt <- qs2::qs_read(input_path)
setDT(dt)

required_cols <- c("eid", "source", "array_index", "diagnosis_icd", "diagnosis_date")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dt <- dt[, ..required_cols]
dt[, eid := as.integer(eid)]
dt[, source := fifelse(is.na(source) | source == "", "unknown", source)]
dt[, array_index := as.character(array_index)]
dt[, diagnosis_icd := toupper(gsub("[^A-Z0-9]", "", diagnosis_icd))]
dt[, diagnosis_date := as.IDate(diagnosis_date)]
dt <- dt[!is.na(eid) & diagnosis_icd != "" & !is.na(diagnosis_date)]

dt[, icd10_3 := substr(diagnosis_icd, 1L, pmin(3L, nchar(diagnosis_icd)))]
dt[, icd_letter := substr(diagnosis_icd, 1L, 1L)]
dt[, icd_num := suppressWarnings(as.integer(substr(diagnosis_icd, 2L, 3L)))]

chapter_label <- function(letter, num) {
  fifelse(letter %in% c("A", "B"), "infectious",
  fifelse(letter == "C" | (letter == "D" & num <= 49), "neoplasm",
  fifelse(letter == "D" & num >= 50, "blood_immune",
  fifelse(letter == "E", "endocrine_metabolic",
  fifelse(letter == "F", "mental_behavioral",
  fifelse(letter == "G", "nervous_system",
  fifelse(letter == "H", "eye_ear",
  fifelse(letter == "I", "circulatory",
  fifelse(letter == "J", "respiratory",
  fifelse(letter == "K", "digestive",
  fifelse(letter == "L", "skin",
  fifelse(letter == "M", "musculoskeletal",
  fifelse(letter == "N", "genitourinary",
  fifelse(letter == "O", "pregnancy",
  fifelse(letter == "P", "perinatal",
  fifelse(letter == "Q", "congenital",
  fifelse(letter == "R", "symptoms_signs",
  fifelse(letter %in% c("S", "T"), "injury_poisoning",
  fifelse(letter %in% c("V", "W", "X", "Y"), "external_causes",
  fifelse(letter == "Z", "health_service", "other"))))))))))))))))))))
}

site_label <- function(code3) {
  letter <- substr(code3, 1L, 1L)
  num <- suppressWarnings(as.integer(substr(code3, 2L, 3L)))
  out <- rep(NA_character_, length(code3))
  out[letter == "C" & num >= 0 & num <= 14] <- "lip_oral_pharynx"
  out[letter == "C" & num >= 15 & num <= 26] <- "digestive"
  out[letter == "C" & num >= 30 & num <= 39] <- "respiratory_intrathoracic"
  out[letter == "C" & num >= 40 & num <= 41] <- "bone_cartilage"
  out[letter == "C" & num >= 43 & num <= 44] <- "skin"
  out[letter == "C" & num >= 45 & num <= 49] <- "mesothelial_soft_tissue"
  out[letter == "C" & num == 50] <- "breast"
  out[letter == "C" & num >= 51 & num <= 58] <- "female_genital"
  out[letter == "C" & num >= 60 & num <= 63] <- "male_genital"
  out[letter == "C" & num >= 64 & num <= 68] <- "urinary"
  out[letter == "C" & num >= 69 & num <= 72] <- "eye_brain_cns"
  out[letter == "C" & num >= 73 & num <= 75] <- "thyroid_endocrine"
  out[letter == "C" & num >= 76 & num <= 80] <- "secondary_ill_defined"
  out[letter == "C" & num >= 81 & num <= 96] <- "haematological"
  out[letter == "D" & num >= 0 & num <= 9] <- "in_situ"
  out[letter == "D" & num >= 10 & num <= 36] <- "benign_neoplasm"
  out[letter == "D" & num >= 37 & num <= 48] <- "uncertain_neoplasm"
  out
}

dt[, icd_chapter := chapter_label(icd_letter, icd_num)]
dt[, cancer_site_broad := site_label(icd10_3)]

dt[, event_class := fcase(
  icd_letter == "C" & icd_num %between% c(77L, 79L), "metastasis",
  icd_letter == "C" & icd_num == 80L, "cancer_unspecified",
  icd_letter == "C" & icd_num %between% c(0L, 76L), "primary_cancer",
  icd_letter == "C" & icd_num %between% c(81L, 96L), "haematological_cancer",
  icd_letter == "D" & icd_num %between% c(0L, 9L), "in_situ_neoplasm",
  icd_letter == "D" & icd_num %between% c(37L, 48L), "uncertain_neoplasm",
  icd_letter == "D" & icd_num %between% c(10L, 36L), "benign_neoplasm",
  default = "non_neoplasm"
)]

dt[, is_precancer_flag := event_class %in% c("in_situ_neoplasm", "uncertain_neoplasm")]
dt[, is_malignant_flag := event_class %in% c("primary_cancer", "haematological_cancer", "metastasis", "cancer_unspecified")]
dt[, is_metastasis_flag := event_class == "metastasis"]
dt[, is_neoplasm_flag := event_class != "non_neoplasm"]
dt[, token_icd3 := paste("ICD3", icd10_3, sep = "_")]
dt[, token_stage := paste("DX", event_class, icd10_3, sep = "_")]

source_priority <- c(cancer_registry = 1L, inpatient = 2L, unknown = 9L)
dt[, source_rank := fifelse(source %in% names(source_priority), source_priority[source], 8L)]
dt[, array_rank := suppressWarnings(as.integer(array_index))]
dt[is.na(array_rank), array_rank := 999999L]

setorder(dt, eid, diagnosis_date, source_rank, array_rank, diagnosis_icd)
dt <- unique(dt, by = c("eid", "source", "diagnosis_date", "diagnosis_icd"))
dt[, event_index_all := seq_len(.N), by = eid]
dt[, event_index_day := seq_len(.N), by = .(eid, diagnosis_date)]
dt[, days_since_first_dx := as.integer(diagnosis_date - min(diagnosis_date)), by = eid]

first_cancer <- dt[is_malignant_flag == TRUE, .(first_malignant_date = min(diagnosis_date)), by = eid]
dt <- first_cancer[dt, on = "eid"]
dt[, days_to_first_malignant := as.integer(diagnosis_date - first_malignant_date)]
dt[, sequence_phase := fcase(
  is.na(first_malignant_date), "no_malignant_recorded",
  diagnosis_date < first_malignant_date, "pre_malignant",
  diagnosis_date == first_malignant_date & is_malignant_flag, "index_malignant",
  diagnosis_date == first_malignant_date, "same_day_index_context",
  diagnosis_date > first_malignant_date & is_metastasis_flag, "post_malignant_metastasis",
  diagnosis_date > first_malignant_date, "post_malignant",
  default = "unknown"
)]

vocab <- dt[, .(
  n_events = .N,
  n_eid = uniqueN(eid),
  first_date = min(diagnosis_date),
  last_date = max(diagnosis_date),
  icd10_3 = first(icd10_3),
  event_class = first(event_class),
  cancer_site_broad = first(cancer_site_broad)
), by = .(token = token_stage)]
setorder(vocab, -n_events, token)
vocab[, token_id := .I]
setcolorder(vocab, c("token_id", "token"))

dt <- vocab[, .(token_stage = token, token_id)][dt, on = "token_stage"]
setorder(dt, eid, event_index_all)

person_sequence <- dt[, .(
  n_events = .N,
  n_unique_icd3 = uniqueN(icd10_3),
  first_dx_date = min(diagnosis_date),
  last_dx_date = max(diagnosis_date),
  first_malignant_date = first(first_malignant_date),
  has_precancer = any(is_precancer_flag),
  has_malignant = any(is_malignant_flag),
  has_metastasis = any(is_metastasis_flag),
  tokens_all = paste(token_id, collapse = " "),
  tokens_neoplasm = paste(token_id[is_neoplasm_flag], collapse = " "),
  dates_all = paste(as.character(diagnosis_date), collapse = " ")
), by = eid]

summary_tables <- list(
  source = dt[, .N, by = source][order(-N)],
  event_class = dt[, .(n_events = .N, n_eid = uniqueN(eid)), by = event_class][order(-n_events)],
  phase = dt[, .(n_events = .N, n_eid = uniqueN(eid)), by = sequence_phase][order(-n_events)],
  person_flags = person_sequence[, .N, by = .(has_precancer, has_malignant, has_metastasis)][order(-N)]
)

events_out <- file.path(output_dir, "diagnosis_events_clean.qs2")
seq_out <- file.path(output_dir, "person_diagnosis_sequences.qs2")
vocab_out <- file.path(output_dir, "diagnosis_token_vocab.csv")
summary_out <- file.path(output_dir, "diagnosis_cleaning_summary.txt")

message("Writing outputs to: ", output_dir)
qs2::qs_save(dt, events_out)
qs2::qs_save(person_sequence, seq_out)
fwrite(vocab, vocab_out)

if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(dt, file.path(output_dir, "diagnosis_events_clean.parquet"))
  arrow::write_parquet(person_sequence, file.path(output_dir, "person_diagnosis_sequences.parquet"))
}

sink(summary_out)
cat("Input:", input_path, "\n")
cat("Output directory:", output_dir, "\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
cat("Events:", nrow(dt), "\n")
cat("Individuals:", uniqueN(dt$eid), "\n")
cat("Unique raw ICD codes:", uniqueN(dt$diagnosis_icd), "\n")
cat("Unique ICD3 codes:", uniqueN(dt$icd10_3), "\n")
cat("Date range:", as.character(min(dt$diagnosis_date)), "to", as.character(max(dt$diagnosis_date)), "\n\n")
for (nm in names(summary_tables)) {
  cat("##", nm, "\n")
  print(summary_tables[[nm]])
  cat("\n")
}
sink()

message("Done.")
message("Events: ", events_out)
message("Sequences: ", seq_out)
message("Vocabulary: ", vocab_out)
message("Summary: ", summary_out)
