---
description: Guided CV/site refresh — update cv_data.yml, rebuild, commit, and publish
allowed-tools: Bash(./refresh.sh), Bash(git *), Read, Edit, Grep, Glob
---

You are refreshing the live CV/website (stevenhawken.ca). Follow these steps in order.

## 1. Ask what's new

Ask the user what has changed since the last refresh. Possibilities:

- New grants, invited talks, trainees, teaching, honours, editorial/review or
  thesis examiner roles → edit the matching section of `cv_data.yml`.
  YAML is whitespace-sensitive: keep the existing 2-space indentation and
  quoting style; edit only the relevant block.
- A publication the auto-pull missed → add its PMID to `extra_pmids:` in
  `cv_data.yml`. A wrongly attributed paper (different S. Hawken) → add its
  PMID to `exclude_pmids:`.
- Prose changes → edit `index.qmd`, `software_and_web_tools.qmd`, or
  `resources.qmd`.
- "Nothing to add, just pull pubs" → skip straight to step 2.

## 2. Build

Run `./refresh.sh`. It clears `_freeze/cv` (required — `freeze: auto` won't
re-execute the CV otherwise), re-pulls publications from Google Scholar and
PubMed, and renders to `docs/`.

If the Google Scholar pull errors (rate limiting), rerun `./refresh.sh` once
or twice before troubleshooting further.

## 3. Review

Show the user:
- `git status --short` and a concise summary of the diff (source files AND
  what changed in `docs/` — e.g., how many new publications appeared).
- Sanity checks: new publications present? citation counts updated? no
  sections accidentally emptied?

## 4. Commit and publish

- Stage everything (`git add -A`) — the rendered `docs/` MUST be committed;
  GitHub Pages serves the committed `docs/` on `main`, there is no CI.
- Commit with a descriptive message summarizing what changed.
- Ask the user to confirm, then `git push origin main`.
- Remind them the live site updates in ~1–2 minutes at
  https://www.stevenhawken.ca.

If `git push` prompts for credentials, stop and tell the user to set up
`gh auth login` or an SSH key.
