#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
events_path <- if (length(args) >= 1) args[[1]] else "diagnosis_traj/processed/diagnosis_events_clean.qs2"
fam_dir <- if (length(args) >= 2) args[[2]] else "st10_core_rsid/bed_by_chr"
outdir <- if (length(args) >= 3) args[[3]] else "model_data/traj"

if (!requireNamespace("qs2", quietly = TRUE)) {
  stop("Package 'qs2' is required.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading diagnosis events: ", events_path)
dt <- qs2::qs_read(events_path)
setDT(dt)
dt[, diagnosis_icd := toupper(gsub("[^A-Z0-9]", "", diagnosis_icd))]
dt[, diagnosis_date := as.IDate(diagnosis_date)]
dt[, code3 := substr(diagnosis_icd, 1L, pmin(3L, nchar(diagnosis_icd)))]
dt[, letter := substr(diagnosis_icd, 1L, 1L)]
dt[, num := suppressWarnings(as.integer(substr(diagnosis_icd, 2L, 3L)))]

fam_files <- list.files(fam_dir, pattern = "\\.fam$", full.names = TRUE)
if (length(fam_files) == 0) stop("No .fam files found in ", fam_dir)
fam <- rbindlist(lapply(fam_files, fread, header = FALSE))
setnames(fam, c("fid", "iid", "pat", "mat", "sex", "phenotype"))
fam[, eid := as.integer(iid)]
fam <- unique(fam[, .(eid, sex)])
dt <- dt[eid %in% fam$eid]

site_from_cancer <- function(letter, num) {
  fifelse(letter == "C" & num %between% c(0L, 6L), "lip_oral_cavity",
  fifelse(letter == "C" & num %between% c(7L, 8L), "salivary_glands",
  fifelse(letter == "C" & num %between% c(9L, 10L), "oropharynx",
  fifelse(letter == "C" & num == 11L, "nasopharynx",
  fifelse(letter == "C" & num %between% c(12L, 13L), "hypopharynx",
  fifelse(letter == "C" & num == 15L, "esophagus",
  fifelse(letter == "C" & num == 16L, "stomach",
  fifelse(letter == "C" & num == 18L, "colon",
  fifelse(letter == "C" & num %between% c(19L, 20L), "rectum",
  fifelse(letter == "C" & num == 21L, "anus",
  fifelse(letter == "C" & num == 22L, "liver_intrahepatic_bile_ducts",
  fifelse(letter == "C" & num == 23L, "gallbladder",
  fifelse(letter == "C" & num == 25L, "pancreas",
  fifelse(letter == "C" & num == 32L, "larynx",
  fifelse(letter == "C" & num %between% c(33L, 34L), "lung_trachea_bronchus",
  fifelse(letter == "C" & num == 43L, "melanoma_skin",
  fifelse(letter == "C" & num == 44L, "nonmelanoma_skin",
  fifelse(letter == "C" & num == 45L, "mesothelioma",
  fifelse(letter == "C" & num == 46L, "kaposi_sarcoma",
  fifelse(letter == "C" & num == 50L, "female_breast",
  fifelse(letter == "C" & num == 51L, "vulva",
  fifelse(letter == "C" & num == 52L, "vagina",
  fifelse(letter == "C" & num == 53L, "cervix_uteri",
  fifelse(letter == "C" & num == 54L, "corpus_uteri",
  fifelse(letter == "C" & num == 56L, "ovary",
  fifelse(letter == "C" & num == 60L, "penis",
  fifelse(letter == "C" & num == 61L, "prostate",
  fifelse(letter == "C" & num == 62L, "testis",
  fifelse(letter == "C" & num %between% c(64L, 65L), "kidney_renal_pelvis",
  fifelse(letter == "C" & num == 67L, "bladder",
  fifelse(letter == "C" & num %between% c(70L, 72L), "brain_cns",
  fifelse(letter == "C" & num == 73L, "thyroid",
  fifelse(letter == "C" & num == 81L, "hodgkin_lymphoma",
  fifelse(letter == "C" & (num %between% c(82L, 86L) | num == 96L), "non_hodgkin_lymphoma",
  fifelse(letter == "C" & num %in% c(88L, 90L), "multiple_myeloma_immunoproliferative",
  fifelse(letter == "C" & num %between% c(91L, 95L), "leukemia",
  fifelse(letter == "C" & num %between% c(77L, 79L), "metastasis", NA_character_)))))))))))))))))))))))))))))))))))))
}

dt[, cancer_site := site_from_cancer(letter, num)]
dt[letter == "C" & num %between% c(18L, 21L), cancer_site := "colorectal"]

add_precancer <- function(code, code3, letter, num) {
  fifelse(code3 == "D000" | code %chin% c("K132", "K135", "K137"), "lip_oral_cavity",
  fifelse(code3 == "D000", "oropharynx",
  fifelse(code3 == "D000", "hypopharynx",
  fifelse(code3 == "D001" | code == "K227" | code == "K220", "esophagus",
  fifelse(code3 == "D002" | code == "K294" | code == "K318", "stomach",
  fifelse(code3 %chin% c("D010", "D011", "D012", "D013") | code3 %chin% paste0("D12", 0:8) | code %chin% c("K635", "K621"), "colorectal",
  fifelse(code3 == "D015" | substr(code, 1, 3) %chin% c("K74", "B18") | code == "K703", "liver_intrahepatic_bile_ducts",
  fifelse(code3 == "D015" | code == "K824" | substr(code, 1, 3) == "K80", "gallbladder",
  fifelse(code == "D136" | code == "K862" | code3 %chin% c("D017", "D019"), "pancreas",
  fifelse(code3 == "D020" | code == "J383", "larynx",
  fifelse(code3 %chin% c("D021", "D022"), "lung_trachea_bronchus",
  fifelse(substr(code, 1, 3) == "D03", "melanoma_skin",
  fifelse(substr(code, 1, 3) == "D04" | code == "L570", "nonmelanoma_skin",
  fifelse(substr(code, 1, 3) == "D05" | code == "N608", "female_breast",
  fifelse(code3 == "D071" | code3 %chin% c("N900", "N901", "N902", "N903") | code == "L900", "vulva",
  fifelse(code3 == "D072" | code3 %chin% c("N890", "N891", "N892", "N893"), "vagina",
  fifelse(substr(code, 1, 3) == "D06" | code3 %chin% c("N871", "N872", "N879"), "cervix_uteri",
  fifelse(code3 == "D070" | code %chin% c("N850", "N840"), "corpus_uteri",
  fifelse(code == "D391" | substr(code, 1, 3) == "N80", "ovary",
  fifelse(code3 == "D074" | code == "N480", "penis",
  fifelse(code3 == "D075" | code == "N423", "prostate",
  fifelse(code3 == "D076", "testis",
  fifelse(code3 == "D091", "kidney_renal_pelvis",
  fifelse(code3 == "D090", "bladder",
  fifelse(code == "D34" | code == "E312", "thyroid",
  fifelse(code == "D472" | code == "D477", "multiple_myeloma_immunoproliferative",
  fifelse(substr(code, 1, 3) == "D46" | code %chin% c("D471", "D473", "D474"), "leukemia", NA_character_)))))))))))))))))))))))))))
}

dt[, precancer_site := add_precancer(diagnosis_icd, code3, letter, num)]

strict_prefix <- function(code, code3) {
  substr(code, 1, 3) %chin% c("D00", "D01", "D02", "D03", "D04", "D05", "D06") |
    code3 %chin% c("D070", "D071", "D072", "D074", "D075", "D076", "D090", "D091") |
    code %chin% c("N871", "N872", "D391", "D472") |
    substr(code, 1, 3) == "D46" |
    code %chin% c("D471", "D473", "D474")
}

cancer_events <- dt[!is.na(cancer_site) | !is.na(precancer_site),
  .(eid, diagnosis_date, diagnosis_icd, code3, cancer_site, precancer_site, letter, num)]

cancer_events[, stage := fcase(
  !is.na(cancer_site) & cancer_site == "metastasis", "metastasis",
  !is.na(cancer_site), "index_cancer",
  !is.na(precancer_site) & strict_prefix(diagnosis_icd, code3), "precancer_strict",
  !is.na(precancer_site), "risk_or_precancer_broad",
  default = "other"
)]
cancer_events[is.na(cancer_site), cancer_site := precancer_site]
cancer_events <- cancer_events[!is.na(cancer_site)]
cancer_events[, token := paste(cancer_site, stage, sep = "__")]

setorder(cancer_events, eid, diagnosis_date, token)
cancer_events <- unique(cancer_events, by = c("eid", "diagnosis_date", "token"))
cancer_events[, event_index := seq_len(.N), by = eid]

vocab <- cancer_events[, .N, by = token][order(-N, token)]
special <- data.table(token = c("<PAD>", "<BOS>", "<EOS>", "<UNK>"), N = NA_integer_)
vocab <- rbind(special, vocab, fill = TRUE)
vocab[, token_id := .I - 1L]
setcolorder(vocab, c("token_id", "token", "N"))
cancer_events <- vocab[, .(token, token_id)][cancer_events, on = "token"]
setorder(cancer_events, eid, event_index)

seqs <- cancer_events[, .(
  n_events = .N,
  first_event_date = min(diagnosis_date),
  last_event_date = max(diagnosis_date),
  tokens = paste(token, collapse = " "),
  token_ids = paste(token_id, collapse = " "),
  has_index_cancer = any(stage == "index_cancer"),
  has_metastasis = any(stage == "metastasis"),
  has_precancer_strict = any(stage == "precancer_strict")
), by = eid]
seqs <- fam[, .(eid, sex)][seqs, on = "eid"]

fwrite(cancer_events, file.path(outdir, "cancer_site_events.csv"))
fwrite(seqs, file.path(outdir, "cancer_site_sequences.csv"))
fwrite(vocab, file.path(outdir, "cancer_site_token_vocab.csv"))
qs2::qs_save(cancer_events, file.path(outdir, "cancer_site_events.qs2"))
qs2::qs_save(seqs, file.path(outdir, "cancer_site_sequences.qs2"))

sink(file.path(outdir, "cancer_site_sequence_summary.txt"))
cat("Individuals in fam:", uniqueN(fam$eid), "\n")
cat("Individuals with cancer-site sequence:", nrow(seqs), "\n")
cat("Events:", nrow(cancer_events), "\n")
cat("Tokens:", nrow(vocab), "\n\n")
print(cancer_events[, .N, by = .(cancer_site, stage)][order(-N)][1:80])
sink()

message("Done: ", outdir)
