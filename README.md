# Ovary Macrophage Library

日本語・English対応の、卵巣マクロファージ候補に焦点を置く計算済み解析ライブラリ。Static, read-only catalogue for revisiting single-cell results.

## Included now

- One study: Morris et al. (2022), *A single-cell atlas of the cycling murine ovary*. [eLife](https://doi.org/10.7554/eLife.77239), [source data](https://osf.io/924fz/).
- Default viewer: QC-filtered **1,377 macrophage-candidate cells / 7 groups**, not a definitive cell-type annotation.
- Whole-ovary **54,748 cells / 27 clusters** retained separately for context; **25,289 genes**, **14 saved analyses** total.
- No completed cross-study integration. Historical 915-cell, differential-expression and GSEA reports are not bundled.
- No analysis server, AI model calls, external fonts or tracking. Browser draws stored UMAPs and computes the documented light log-normalization formula from saved counts; it does not recluster or integrate.

## View and publish

Serve this directory as static files. There is no npm dependency or website build service requirement. GitHub Pages can publish the root of a reviewed public repository from its main branch. `.nojekyll` preserves the static files as supplied. All app assets use repository-relative paths; no custom domain is required. Do not upload private analysis source directories, account files, credentials, or unrelated notes.

Optional local preview from this directory: `python3 -m http.server 8041` then open `http://localhost:8041/`.

Use `index.html?lang=en` or `?lang=ja`. Viewer direct links accept `?run=<saved-run-id>&lang=en`. Unknown run IDs fall back to QC1377. For whole ovary use `studies/morris-2022/viewer/?run=395b054678cad32d3e`. Language is remembered where browser storage is available.

`robots.txt` requests no indexing, but is not access control: published files are public.

## Add a study

1. Compute/QC/annotate in the analysis environment, keeping source data and provenance there. Export only a reviewed viewer snapshot under `studies/<study-id>/viewer/`.
2. Add a record to `data/studies.json`. Preserve the existing bilingual field shapes and include a stable ID, paper title/DOI, source URL, species, biological conditions, candidate-selection caveat, cell/gene counts, methods, version, snapshot date, and saved analysis records with IDs/parent IDs/paths.
3. Adapt the new study's viewer entry/deep links; do not relabel old coordinates as new analyses. Supply study-specific bilingual `page_title`, `included_scope` and `interpretation` fields; do not reuse Morris-specific biological claims. The JSON is an extensible schema, not a scientific analysis engine.
4. Run `python3 scripts/build.py` to regenerate catalogue pages. No scientific assets are changed by this command.
5. Run `python3 scripts/verify.py`, `node scripts/check-viewer.cjs`, `node scripts/check-stage-filter.cjs` and JavaScript syntax checks. Then manually verify the deployed page, language switch, candidate default, gene expression and paths in a browser. Hash verification is not visual verification.
6. Publish updated static files. Keep prior snapshot data and the source archive in the research backup, and preserve stable URL names when possible.

## Add integrated results

`data/integrations.json` starts with an empty `integrations` array. Add an entry only when the result exists and is validated. Each entry must contain the documented required fields, a bilingual title, contributing study IDs, species, biological conditions, counts, methods/batch-correction choices, donor metadata, version, date and a working page path. Put its reviewed viewer/result page under `integrations/<integration-id>/` before running the renderer. Do not present a side-by-side comparison as batch-corrected joint analysis. Retain uncorrected expression for biological interpretation and donor-aware inference.

## Third-party code

Plotly is bundled only inside the existing study viewer. Its license is included as `studies/morris-2022/viewer/PLOTLY_LICENSE.txt`. Scientific data and original paper rights remain with their respective sources; this repository does not relicense them.

## Stage-filter semantics

Each saved UMAP panel has a stage filter (all/proestrus/estrus/metestrus/diestrus/unknown), a label-source selector, and an evidence selector. Existing expression panels and expression CSV export follow the selected cells. UMAP coordinates and original normalization denominators remain unchanged. Selecting a saved child inherits these filters; the parent remains visible. Filtering resets lasso selections to avoid exporting hidden cells. Explicit palette/range settings remain intact.

QC1377 defaults to **analysis stages (including inferred assignments)**, explicitly labelled: 317 proestrus, 486 estrus, 193 metestrus, 381 diestrus. There are **927 direct author-cell matches and 450 library-supported inferred matches**. Switch to author-mapped stage to retain those 450 as **unknown**, or select direct author matches only. Direct-only stage counts are 217 / 266 / 161 / 283, respectively. These are cell counts, not independent biological replicate counts.

Missing per-cell stage fields in a child are inherited from its saved parent using source row IDs, never by positional coincidence. Whole-ovary metadata has no explicit `stage_evidence` field; its evidence is shown as unrecorded and direct-only selection is unavailable. No claim of confirmed evidence is inferred from a stage label alone. Each gene panel also has **「4周期を並べて比較」 / “Compare four stages”**. It opens proestrus/estrus/metestrus/diestrus together, using the same saved coordinates and fixed full-analysis axis ranges. All four panels share the gene’s full-analysis automatic color limits, or the explicit minimum/maximum entered for that gene. Palette, range, language, stage-label source and evidence changes refresh the comparison. The single-stage selector still controls the original expression plot and CSV only; it does not remove stages from the comparison. Author-unknown cells are not assigned to any of the four named stages. Empty stage panels are explicit. Closing the comparison, gene or analysis, or resetting the view releases its Plotly plots. No clustering, coordinates or scientific files are recalculated.

The added comparison is covered by non-browser contracts; actual browser rendering and mobile interaction remain unverified.
