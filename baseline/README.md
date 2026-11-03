# Output baseline (2026-11-02)

Reference fingerprints for checking that refactoring does not change results (policy §3). *AI-drafted; log entry
2026-11-02-16.*

- `output_sha256_run1.txt`: SHA-256 of all 93 files in `output/` after the first full `run_all.R` run
  (8 CSVs in `output/tables/`, 85 files in `output/plots/`).
- `check_baseline.R`: compares the current CSVs with those fingerprints.
- The files themselves are **not** in git (`output/` is gitignored). A copy of run 1 is kept outside the repository on
  the researcher's machine only (`baseline_outputs_2026-11-02_run1/`).

## Repeatability check (run 1 vs. a second identical run, same machine and settings)

- **CSVs: 8 of 8 byte-identical.** The numeric results were reproduced exactly.
- **PDFs: 21 of 21 identical** once embedded creation/modification dates are ignored (raw bytes differ by ~10 bytes
  each because of those dates). Use the date-stripped comparison, not raw checksums, for PDFs.
- **PNGs: 56 of 64 identical; 8 differ**, all `output/plots/labelled_peaks/<year>/<site>.png`. **Reasoned, not
  confirmed:** these come from `ggrepel::geom_text_repel`, which places labels with a random start, and no seed is
  set anywhere in the scripts (`set.seed` appears 0 times). This would change label positions only, not data.
  Setting a `seed` in that call would make them repeatable; not done, as it is a plotting change for the researcher
  to approve.

Not covered: different machines, R/package versions other than those in `renv.lock`, or runs with changed code.
