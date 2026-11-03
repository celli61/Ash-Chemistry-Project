# What the code in this project does (plain-language summary)

> **AI-drafted (Claude Sonnet 5.5), 2026-11-02, log entry 2026-11-02-12.** Written from reading the R scripts and
> their comments only; no data was opened and the pipeline's results are not described. The researcher must check it
> against the code before relying on it or using any of it in a thesis or disclosure.

## The research question the scripts belong to

The header of `NH_Ash_Chemistry_Alignments.R` states the project goal: to test the ontogenetic resistance
hypothesis in green and white ash trees across four size classes in New Hampshire, by looking at how the mix of
defensive chemicals in the phloem (inner bark) changes after emerald ash borer attack or after treatment with the
plant hormone methyl jasmonate. The scripts here cover the **chemistry data processing** part. According to a note
in the header, the later steps (NMDS ordination and Random Forest) were cut out of these files, so no statistical
testing happens in this code.

## The big picture

Plant samples were run through an HPLC instrument (Shimadzu, diode-array detector), which produces a
**chromatogram** for each sample: a curve showing how much light is absorbed (here at 280 nm) as different
chemicals come off the column over time. Each bump ("peak") is roughly one chemical. The time at which a peak
appears is its **retention time**. The code:

1. reads those raw files,
2. lines up the chromatograms from different samples so the same chemical appears at the same time in each,
3. finds the peaks,
4. compares peaks to known chemical **standards** run on the same instrument to see which known compounds seem to
   be present,
5. measures how large the peaks are, and
6. makes tables and figures.

## The scripts, in the order they run (`run_all.R`)

**1. `NH_Ash_Chemistry_Alignments.R` (the big one)**
- Reads the 2020 and 2021 raw chromatogram files from the data folder (found through `data_path.R`).
- Separates the runs into standards and samples, and samples by field site. Sites named in the code: 2020 – Doe
  Farm, French Conservation Area, East Foss Farm, Jennings Forest; 2021 – Tuttle Swamp, Powder Major's forest,
  Town of Lee forest, Little River Conservation Area.
- Cleans each chromatogram (removes baseline noise and trims to the useful time window).
- Aligns retention times between samples (`correct_rt`) so peaks line up.
- Finds peaks (`get_peaks`, "EGH" peak shape) using settings chosen by hand for each group of samples and each
  standard, builds a table of peaks by sample (`get_peaktable`), and attaches each sample's information (site,
  treatment, tree size, etc.) from the metadata.

**2. `NH_Ash_Chemistry_Standards.R`**
- Records, for each known standard compound, the retention time and a tolerance window, in tables the later steps
  can loop over. In 2020 these are 12 standards (gallic, chlorogenic and vanillic acids, eleutheroside B, coumaric
  and ferulic acids, verbascoside, oleuropein, pinoresinol, ligustrocide, luteolin, apagenin). The 2021 section
  also refers to syringin.

**3. `NH_Ash_Chemistry_Standard_Curves.R`**
- Processes the "standard mix" runs (several standards together, at several known concentrations such as 1, 0.5,
  0.25 and 0.125 mg/ml): aligns them, finds their peaks, and builds standard-mix peak tables.
- Collects known concentrations alongside measured peak values for each standard. The comments say this is to
  build calibration models for estimating unknown concentrations; I did not find any `lm()` call in the scripts,
  so whether a regression is actually fitted is **unconfirmed**.

**4. `NH_Ash_Chemistry_Standard_Algorithmic_Detection.R`**
- Defines helper functions: a search for a retention time within a tolerance (`linear_search`), a table of each
  peak's strongest absorbance wavelength, and several versions of a **standard presence table**.
- A presence table says, for each sample, whether a peak was found near each standard's retention time (and, in
  one version, whether its spectrum also matches). Then it builds presence-rate tables per site.
- Quantifies peaks (`quantify_peaks`: the comments describe area-under-the-curve contributions) and writes one
  CSV per site, `<site>_peak_and_meta.csv`, to `output/tables/`. These are the only table files the scripts
  write.

**5. `NH_Ash_Chemistry_Plotting.R`**
- Defines plotting functions and then saves about 88 figures to `output/plots/` (colour-blind-safe Okabe-Ito
  palette): a chromatogram per sample marked where a standard was matched, presence of standards by tree size and
  by treatment, detection rates for different detection methods, peak quantity plots, and labelled chromatograms.

**6. `NH_Ash_Chemistry_Absorbance.R`**
- Defines a small function that turns each peak's intensity into an absorbance value (log10 of a baseline
  intensity divided by the peak intensity) and runs it on the first standard mix and the Doe Farm samples. It
  prints results; nothing is saved.

## Scripts and files outside the main flow

- `NH_Ash_Chemistry_Alignments_2019.R` repeats the alignment for 2019 Doe Farm samples. It reuses object names
  from script 1, so it is not part of `run_all.R`; its role is undecided.
- `run_all.R` runs scripts 1–6 in order in one R session (they pass objects to each other through R's memory,
  not through saved files).
- `data_path.R` finds the raw data on whatever computer you are using (set in a local `.Renviron`).
- `renv.lock` fixes the R package versions. `AI_USE_POLICY.md` and `AI_LOG.md` hold the AI-use rules and record.

## Things to check or know

- Many numbers (peak-finding thresholds, smoothing, tolerance windows) were tuned by hand and differ by group.
  The reasons are only in scattered code comments.
- Comments in the code note that some chemicals may be split into two peaks by the alignment.
- Sample numbers used inside the plotting calls (for example run index 43 for Doe Farm) are typed in by hand.
- This summary describes what the code is set up to do. Whether each step worked as intended on the real data has
  not been verified by me.
