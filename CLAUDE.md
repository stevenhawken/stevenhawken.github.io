# Website & Dynamic CV — Maintenance & Refresh Workflow

Personal site + dynamic CV. **Live:** https://www.stevenhawken.ca
**Repo:** stevenhawken.github.io (this folder)

## Architecture (read before changing anything)
- Quarto **website** project. Renders to `docs/` (`output-dir: docs` in `_quarto.yml`).
- Deploys via **GitHub Pages from the `main` branch `/docs` folder. No CI / no Action.**
  → The rendered `docs/` MUST be committed for the live site to change.
- Custom domain via `CNAME` (www.stevenhawken.ca); post-render hooks copy `CNAME`
  into `docs/` and touch `docs/.nojekyll`.
- `execute: freeze: auto`, and `_freeze/` is gitignored (cache is local-only).

## How the CV is built (`cv.qmd`)
- Manual sections come from **`cv_data.yml`** (education, employment, grants,
  trainee_supervision, invited_talks_seminars, academic_honours, teaching,
  thesis_* examiner roles, editorial_review, preprints, etc.).
- Publications/citations are pulled **live at render time** from Google Scholar
  (`scholar`, ID `aTyysxgAAAAJ`) and PubMed (`rentrez`).
- Manual publication control: `extra_pmids:` forces in a paper the auto-pull
  misses; `exclude_pmids:` drops a wrong-attribution (a different S. Hawken).

## ⚠️ The freeze gotcha
`freeze: auto` only re-executes when the `.qmd` *source* changes — NOT when
`cv_data.yml` or the live Scholar/PubMed data changes. So a plain `quarto render`
silently reuses stale CV output. **`refresh.sh` deletes `_freeze/cv` every run**,
which is what makes new publications AND cv_data.yml edits take effect.

## Refreshing (the normal workflow)
Open Claude Code in this repo and run **`/refresh-cv`**. It:
1. Asks what's new and edits the right `cv_data.yml` sections (or `extra_pmids`),
   plus any prose pages.
2. Runs `./refresh.sh` (clears `_freeze/cv` → re-pulls pubs → renders to `docs/`).
3. Shows the diff + a summary, commits with a descriptive message.
4. Asks for confirmation, then pushes to `main`. Live site updates in ~1–2 min.

Publication-only refresh: run `/refresh-cv` and say "nothing to add, just pull pubs."

## Key files
- `_quarto.yml` — site config, nav, output-dir, post-render hooks
- `cv.qmd` / `cv_data.yml` — CV render logic / CV data + PMID overrides
- `index.qmd`, `software_and_web_tools.qmd`, `resources.qmd` — prose pages
- `refresh.sh` — deterministic build (clear freeze → render → status)
- `.claude/commands/refresh-cv.md` — the guided refresh+publish command

## Troubleshooting
- **Scholar pull errors** — Google Scholar rate-limits scrapers; just rerun.
  PubMed (`rentrez`) is more reliable.
- **`git push` prompts for credentials** — auth isn't set up for non-interactive
  use; fix with `gh auth login` or an SSH key before relying on the command.
- **Live site didn't change** — confirm `docs/` was actually committed and pushed
  to `main` (Pages serves the committed `docs/`, not your local render).
- **Editing `cv_data.yml`** — YAML is whitespace-sensitive; keep the existing
  2-space indentation and quoting; edit only the relevant block.