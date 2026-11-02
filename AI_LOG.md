# AI Use Log

Append-only record of AI assistance in this project. Claude Code adds an entry at the end of each qualifying task (see `AI_USE_POLICY.md` §6); the researcher completes **Researcher review**. Entries are never edited or deleted — corrections are new entries referencing the original ID.

Entry IDs are `YYYY-MM-DD-NN` (NN = order within that day). Commits containing AI work reference the ID in their message.

**Project model:** Claude Code with Claude Sonnet 5.5 (`claude-sonnet-5-5`), pinned in `.claude/settings.json`. Any change of model is recorded as its own entry.

**Advisor authorization:** _(researcher fills in: who authorized AI use, date, agreed scope, and where it is documented — e.g. email, prospectus, meeting notes)_

**Data cleared for AI use:** see `AI_USE_POLICY.md` §12

## Template

```
### [YYYY-MM-DD-NN] Short title
- **Tool / model:** Claude Code — Claude Sonnet 5.5 (`claude-sonnet-5-5`)
- **Request:** what the researcher asked, briefly
- **Purpose:** code writing / debugging / refactoring / code review / method explanation / writing feedback / disclosure compilation / other
- **What was done:** specific actions — files created/changed, functions, commands run
- **Files changed:** `path/to/file.R` (functions/sections), …  — or "None"
- **Affects results:** Yes / No — one-line reason
- **Verification:** what was run and what was observed — or "Not verified" and why
- **Methodological decisions:** options presented and which the researcher chose — or "None"
- **AI-drafted prose:** Yes (where) / No
- **Uncertainties / [UNVERIFIED]:** anything not confirmed — or "None"
- **Rule overrides:** any policy rule the researcher waived for this task — or "None"
- **Commit:** `[AI] … (log YYYY-MM-DD-NN)` — or "Not committed"
- **Researcher review:** _accepted / modified / rejected — notes_
```

---

## Entries

### [2026-11-02-00] AI use policy and log created
- **Tool / model:** Claude (claude.ai chat) — Claude Opus 5.5
- **Request:** Draft a ruleset for Claude Code covering proper AI use and record-keeping for this research project, compliant with LSU requirements
- **Purpose:** other — policy and documentation drafting
- **What was done:** Drafted `AI_USE_POLICY.md` (rules Claude Code loads via `CLAUDE.md`) and this log with its template; revised in the same session to pin Claude Sonnet 5.5 as the project model and to implement LSU requirements (Code of Student Conduct §10.1, Faculty Senate GAI Guidelines Dec 2025, Graduate School AI disclosure guidance), compiled from public LSU web sources on 2026-11-02
- **Files changed:** `AI_USE_POLICY.md`, `AI_LOG.md` (new)
- **Affects results:** No — documentation only
- **Verification:** Not applicable
- **Methodological decisions:** None
- **AI-drafted prose:** Yes — both files are AI-drafted; researcher to review and edit
- **Uncertainties / [UNVERIFIED]:** LSU has no finalized research-specific AI policy; requirements are drawn from the closest applicable sources and should be confirmed with the advisor. Advisor authorization and cleared/restricted data paths not yet recorded (policy §12).
- **Rule overrides:** None
- **Commit:** Not committed
- **Researcher review:** _

### [2026-11-02-02] Updated policy review, baseline tag, and audit reports
- **Tool / model:** Claude Code — Claude Sonnet 5.5 (`claude-sonnet-5-5`)
- **Request:** Re-read the updated `AI_USE_POLICY.md` and `AI_LOG.md`; record `data/raw/` as the confirmed raw-data location and gitignore it; tag and commit the baseline (step 3a) and push; report on manual steps (3e), randomness (3f) and authorship (4); propose a plan for path fixes (3b) without editing code
- **Purpose:** other — project setup and git baseline; code review (read-only audit of the R scripts)
- **What was done:**
  - Re-read both replaced files in full.
  - Edited `AI_USE_POLICY.md` §2: replaced the "Researcher: adjust this path" note with a statement that `data/raw/` (inside the project, gitignored) is the confirmed raw-data location (researcher-approved edit).
  - Appended `data/raw/` to `.gitignore`.
  - Created annotated tag `pre-ai-baseline` on commit `651718d` (HEAD before any AI work, per the researcher).
  - Committed `AI_USE_POLICY.md`, `AI_LOG.md`, `CLAUDE.md`, `.claude/settings.json`, `.gitignore` as one commit, then pushed the commit and the tag to `origin`.
  - Read-only audit of the `*.R` scripts and `git log` for the 3e/3f/4 reports. No data files were opened (policy §2/§12: nothing cleared yet); no R script was run; no R code was edited.
- **Files changed:** `AI_USE_POLICY.md` (§2 one bullet), `.gitignore` (added `data/raw/`), `AI_LOG.md` (this entry); `CLAUDE.md` and `.claude/settings.json` committed unchanged from entry 2026-11-02-01
- **Affects results:** No — no analysis code or data touched
- **Verification:** `git tag -n1` and `git rev-parse` confirmed the tag points at `651718d`. Audit findings are from reading script text and grep output; nothing was executed.
- **Methodological decisions:** None
- **AI-drafted prose:** No
- **Uncertainties / [UNVERIFIED]:**
  - Entry `2026-11-02-01` (setup and conformance check, same day) is absent from the replaced `AI_LOG.md`. The replacement file contained only entry `-00`. Per the append-only rule I have not re-created it; its content is: created `CLAUDE.md` and `.claude/settings.json`, and ran a read-only audit of the scripts (no analysis affected). The researcher should decide whether to restore it.
  - "Before any AI involvement" for `pre-ai-baseline` is the researcher's statement; git history cannot confirm it (all commits are authored by `celli61`).
  - Advisor authorization (log header) and `AI_USE_POLICY.md` §12 cleared/restricted paths are still blank; the researcher was reminded.
- **Rule overrides:** None
- **Commit:** `[AI] Add AI use policy, log, Claude Code config; gitignore data/raw (log 2026-11-02-02)`
- **Researcher review:** _
