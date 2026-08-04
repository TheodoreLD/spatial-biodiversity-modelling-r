# =============================================================================
# 00_standardization.R — taxonomic name standardization.
# Resolves the spiny-trunk species names against the World Flora Online (WFO) and
# the World Checklist of Vascular Plants (WCVP), with synonym resolution via the
# Taxonomic Name Resolution Service (TNRS package; Maitner & Boyle 2024), for
# nomenclatural consistency across sources.
#
# NB TNRS queries a web service -> this step needs an internet connection.
# Input : data/spiny_trunk_species.xlsx  (column `Species`)   [paths$species_xlsx]
# Output: data/processed/species_standardized.csv             [paths$data_processed]
#
# SCOPE: this step reports and documents the nomenclature; it does not feed the spatial
# pipeline. 01_biogeography.R reads the per-species range files by filename, so the
# analysis uses the names as compiled in the xlsx. The CSV written here is the record
# that those names resolve cleanly against WCVP/WFO, for the manuscript and reviewers.
# =============================================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(readr); library(TNRS)
})
source("R/config.R")
create_project_dirs()

infile  <- require_file(paths$species_xlsx, "spiny-trunk species list")
outfile <- file.path(paths$data_processed, "species_standardized.csv")

# ---- submitted names ---------------------------------------------------------
raw_names <- read_excel(infile)$Species |>
  str_replace_all("_", " ") |> str_squish()
raw_names <- unique(raw_names[!is.na(raw_names) & raw_names != ""])
message(length(raw_names), " unique names submitted to TNRS (WCVP + WFO)")

# ---- resolve to WCVP + WFO via TNRS (best single match per name) --------------
res <- TNRS(taxonomic_names = raw_names, sources = c("wcvp", "wfo"), matches = "best")

standardized <- res |>
  transmute(
    submitted_name   = Name_submitted,
    accepted_species = Accepted_species,      # binomial used downstream
    accepted_name    = Accepted_name,         # full accepted name (with author)
    taxonomic_status = Taxonomic_status,      # Accepted / Synonym / No opinion ...
    match_source     = Source,                # wcvp / wfo
    overall_score    = Overall_score
  ) |>
  arrange(submitted_name)

write_csv(standardized, outfile)
message(sprintf("standardization complete -> %s  (%d names; %d not flagged 'Accepted')",
                outfile, nrow(standardized),
                sum(standardized$taxonomic_status != "Accepted", na.rm = TRUE)))
