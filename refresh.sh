#!/usr/bin/env bash
# Deterministic site/CV rebuild.
# Clears the CV freeze cache (freeze: auto won't re-execute cv.qmd when only
# cv_data.yml or live Scholar/PubMed data changed), then renders to docs/.
set -euo pipefail
cd "$(dirname "$0")"

echo "==> Clearing _freeze/cv (forces CV re-execution)"
rm -rf _freeze/cv

echo "==> Rendering site (pulls publications live from Scholar + PubMed)"
quarto render

echo "==> Done. Changed files:"
git status --short

echo ""
echo "Next: review changes, then commit docs/ and push to main:"
echo "  git add -A && git commit -m 'Update CV' && git push origin main"
