# Compare the current output/tables/*.csv with the baseline fingerprints (SHA-256)
# recorded from the first full run of run_all.R (2026-11-02).
# Usage (after running run_all.R): source("baseline/check_baseline.R")
# Needs R >= 4.5 (tools::sha256sum). Only the CSVs are compared strictly; see
# baseline/README.md for why plot files are not.

manifest <- read.table(here::here("baseline", "output_sha256_run1.txt"),
                       col.names = c("sha256", "path"), stringsAsFactors = FALSE)
csv <- manifest[grepl("^output/tables/.*[.]csv$", manifest$path), ]
current <- unname(tools::sha256sum(here::here(csv$path)))
result <- data.frame(path = csv$path, identical = !is.na(current) & current == csv$sha256)
print(result, row.names = FALSE)
if (all(result$identical)) {
  message("All ", nrow(result), " baseline CSVs are byte-identical.")
} else {
  warning("Some baseline CSVs differ or are missing; outputs changed. Stop and report (policy section 3).")
}
