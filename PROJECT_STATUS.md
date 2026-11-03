# Project status

> **AI-drafted (Claude Sonnet 5.5); last updated 2026-11-02, log entry 2026-11-02-16.** A short, factual snapshot so a
> new session or a new device can see where the project stands. The detailed record is `AI_LOG.md`; what the code
> does is in `PROJECT_SUMMARY.md`. Keep this file current whenever a milestone, decision or open item changes.
> It contains no interpretation of results.

## Where things stand

- The R pipeline (`run_all.R`, six scripts) has been run end to end once, on the real data, on the researcher's
  machine (log 2026-11-02-13). It exited with status 0 and produced 8 CSVs (`output/tables/`) and 85 plot files
  (`output/plots/`). `output/` is gitignored.
- Raw data is outside the repo (~20 GB), found through `ASH_CHEM_DATA_DIR` in a local `.Renviron`
  (see `.Renviron.example`, `data_path.R`). It is listed as a cleared path in `AI_USE_POLICY.md` §12.
- Repo is set up with relative paths, an `renv` lockfile, a driver script, dependency headers in each script, and
  the AI-use policy and log.
- **Baseline set (log 2026-11-02-16):** `baseline/` holds SHA-256 fingerprints of run-1 outputs and
  `check_baseline.R`. A second identical run reproduced all 8 CSVs byte for byte (PDFs match once dates are
  ignored; 8 label-placement PNGs differ). Use `source("baseline/check_baseline.R")` after every refactor.
- Statistical steps named in the original plan (NMDS, Random Forest) are not in this repo.

## Pending checks (not done)

- [ ] **Comparison with earlier historical results (policy §3):** none done; no earlier outputs are on file. Needs
      the researcher's older outputs if they exist.
- [ ] **34 warnings** from the runs were not captured or reviewed.
- [ ] **Decision for the researcher:** 8 `labelled_peaks/*.png` differ between runs (probably `ggrepel` label
      placement with no seed set). Whether to set a seed in that plot call is not decided.
- [ ] **Clean-machine test:** `renv::restore()` and the Python environment that `chromConverter` installed on first
      use are untested/unrecorded.

## Open questions for the researcher

- Role of `NH_Ash_Chemistry_Alignments_2019.R` (excluded from `run_all.R`).
- Whether to restore the missing log entry 2026-11-02-01 (setup, same day; it was lost when the log was replaced).
- Whether the pinned GitHub versions of `chromatographR` and `chromConverter` are the intended ones.
- Whether calibration regressions are meant to be fitted (comments describe them; no `lm()` found).
- Advisor authorization is still blank in the `AI_LOG.md` header.

## Decisions made so far

- Model pinned to Claude Sonnet 5.5; `data/raw/` and `output/` gitignored.
- Raw-data location set per machine by environment variable (option A, 2026-11-02).
- Scripts run through a driver (`run_all.R`, option A over saved checkpoints).
- No analysis decisions (tests, models, thresholds, exclusions) have been made with AI help.

## Possible next steps (researcher chooses)

See the list given in chat on 2026-11-02: baseline check, warnings review, settings block for hardcoded values,
peak-identification report, quantification options, and the NMDS / Random Forest stage. None is started.
