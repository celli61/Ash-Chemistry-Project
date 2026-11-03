############################################################################
# Run the whole analysis, top to bottom, in a fresh R session.
#
# The scripts share objects through the global environment, so they must run
# in this order (each one needs objects created by the ones before it):
#   1. Alignments                    reads raw data (data_path.R), aligns and detects peaks
#   2. Standards                     needs 1
#   3. Standard_Curves               needs 1, 2
#   4. Standard_Algorithmic_Detection needs 1, 2, 3; writes output/tables/
#   5. Plotting                      needs 1, 2, 4; writes output/plots/
#   6. Absorbance                    needs 1, 3
#
# NH_Ash_Chemistry_Alignments_2019.R is NOT run here. It reuses object names
# from Alignments (e.g. doe_pk_tab, doe_meta) and would overwrite them, and it
# is not yet decided whether it belongs in the pipeline.
#
# Usage: restart R (so the workspace is empty), then source("run_all.R").
############################################################################

if (length(ls(globalenv())) > 0) {
  stop("Workspace is not empty. Restart R (Session > Restart R) and run again.")
}

scripts <- c(
  "NH_Ash_Chemistry_Alignments.R",
  "NH_Ash_Chemistry_Standards.R",
  "NH_Ash_Chemistry_Standard_Curves.R",
  "NH_Ash_Chemistry_Standard_Algorithmic_Detection.R",
  "NH_Ash_Chemistry_Plotting.R",
  "NH_Ash_Chemistry_Absorbance.R"
)

for (s in scripts) {
  message("Running ", s)
  # local = FALSE: scripts hand objects to each other via the global environment
  source(here::here(s), local = FALSE)
}
