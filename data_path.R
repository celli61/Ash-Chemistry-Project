The # data_path(): where the raw data lives on THIS machine.
#
# The raw data is too large for git, so each person keeps their own copy
# (e.g. on an external drive) and tells R where it is, once, in a local
# `.Renviron` file in the project folder (gitignored; see `.Renviron.example`):
#
#   ASH_CHEM_DATA_DIR=E:/ash_chem/raw
#
# The folder must contain the year folders (2020/Seperate_Comma, ...).
# If the variable is not set, the default is the project's own data/raw/.
# The raw data is read-only: nothing in this project writes to this folder.

data_path <- function(...) {
  root <- Sys.getenv("ASH_CHEM_DATA_DIR", unset = "")
  if (!nzchar(root)) root <- here::here("data", "raw")
  if (!dir.exists(root)) {
    stop("Raw data folder not found: ", root,
         "\nIf the data is on an external drive, check that it is plugged in and ",
         "that ASH_CHEM_DATA_DIR in .Renviron points to the right folder ",
         "(restart R after editing .Renviron). See data_path.R.", call. = FALSE)
  }
  file.path(root, ...)
}
