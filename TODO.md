# TODO — Challenges of Incorporating Rotation Information into Crop Insurance Rates

Open issues only. Rebuilt 2026-09-02 (completed P0/sample-size items dropped;
bibliography audit findings added as P0-bib and P4). Ordered by priority:
blockers first, then table/figure regeneration, then inference, then prose
verification, then bibliography.

---

## P0 — Blockers (build-breaking or committed-artifact contradictions)

### Undefined citations in `rotations_losses.tex`

BibTeX keys are case-sensitive; the five below do not resolve and print as `(?)`
with a `Citation undefined` warning. `just_pope_rotations.tex` and
`rotations_insurance.tex` use the correct keys — only `rotations_losses.tex` is
broken.

- [ ] line 79  — `\citet{justPope1979}` → `\citet{Just1979}` (line 95 same file already uses `Just1979`)
- [ ] line 95  — `\citep[QDANN;][]{ma2024}` → `{Ma2024}` (line 131 same file already uses `Ma2024`)
- [ ] line 490 — `\citet{maQDANN2024}` → `\citet{Ma2024}` (figure source caption)
- [ ] line 173 — `\citep{abatzoglou2013}` → `{Abatzoglou2013}` (GRIDMET)
- [ ] line 173 — `\citep{abatzoglou2018}` → `{Abatzoglou2018}` (TerraClimate)

### `rotations_losses.tex` is stale vs the current draft

- [ ] Abstract (line 79) says panel covers **2009–2022**; intro (line 95) says
      **2003–2022**. Reconcile to the actual sample window.
- [ ] Abstract + line 97 say **"four" structural sequence features**, but line 97
      lists three (recent soy, no consecutive soy, moderate total soy harvests),
      and `TODO.md` / `just_pope_rotations.tex` describe the **three-feature**
      Z-vector (nsoy dropped). Fix the count and the feature list.
- [ ] Decide which file is the manuscript of record. README names
      `rotations_losses.tex`; content suggests `just_pope_rotations.tex` is
      current. Align the README and archive the other(s) under `old/`.

---

## P1 — Regenerate tables/figures on the AWC-free control vector, reconcile prose

AWC (`rootznaws_mean`) was dropped from the R control vector. Re-run the
table/figure exports in `corn_analysis_full.R`, then tick each and reconcile the
prose numbers listed.

- [ ] `tables/corn_rot.tex` — regenerated Sep 1; 3/6 spot-checks pass
      (S-S-C-S-S-C=5.94, C-S-C-S-S-C=4.82, C-S-C-C-C-C=−0.77). Still verify
      S-C-C-C-C-C=−1.08, S-C-S-C-S-C=3.63, C-C-C-C-S-C=4.41 against prose.
- [ ] `tables/corn_rci.tex` — table regen fixed (cols (1)/(2) no longer
      byte-identical; headers "No controls" / "Weather and soil controls").
      Sec 5.5 prose STILL STALE: it cites RCI 1.73 negative, 2.00–2.65 ≈
      reference, positives at 2.74/4.24/4.74. Regenerated levels are 1.41, 1.73,
      2, 2.24, 2.45, 2.65, 3, 3.24, 3.46 — no 2.74/4.24/4.74, and 2.24/2.45/2.65
      are all significantly positive. Rewrite lines 525 and 837; clear inline
      TODOs at 492/502/525.
- [ ] `tables/corn_rot_vpd.tex` — verify Sec 5.8 (dry-July −1.29/−1.88; vpd terms).
- [ ] `tables/corn_rci_vpd.tex` — verify Sec 5.9 (positive interactions at RCI 2.00/2.24/2.45).
- [ ] `tables/corn_jp_moments.tex` (Appendix) — verify Sec 5.2 C-C-C-C-S-C
      variance; resolve the "11.3 units" vs −44,885.8 scale mismatch and restate units.
- [ ] `tables/corn_jp_skew.tex` (Appendix) — verify Sec 5.2 skewness values
      (C-S-C-S-C-C −0.31, S-S-C-C-C-C −0.33, C-C-C-S-S-C −0.37).
- [ ] `tables/corn_rci_jp.tex` (Appendix) — verify; confirm omitted RCI reference level.
- [ ] `figures/corn_rot_plot.png`, `corn_var_plot.png`, `corn_jp_plot.png` —
      regenerate from updated coefficients. Resolve `corn_jp_plot` x-axis
      unit/scale (labeled bu/acre but spans ≈ −200/500).
- [ ] `figures/corn_rci_plot.png` — regenerate alongside `corn_rci.tex`.

---

## P2 — Inference and specification consistency

- [ ] **Bootstrap the three-feature variance/skew stages.** Inference of record
      is the three-stage bootstrap (Sec 4.3), not analytic SEs. Re-confirm the
      marginal `soy_gap` variance effect (p<0.10 analytic — may weaken) and the
      skewness significance stars (Sec 5.2 note flags this) against bootstrap draws.
- [ ] **FE inconsistency in VPD tables.** `corn_rot_vpd` / `corn_rci_vpd` use
      field (`tile_field_ID`) FE; the main mean model (eq. 1) uses county (FIPS)
      FE. Re-estimate under a common FE structure or state the difference so
      sequence coefficients are comparable across tables.
- [ ] **VPD double-counting check.** VPD models include both continuous July VPD
      (`vpd_7`) AND the July dry-season bins. Confirm intended; if the bins
      capture the nonlinearity, drop continuous `vpd_7` or state it holds the
      within-bin slope.
- [ ] **VPD coefficient scale/units.** `vpd_6/7/8` (≈ 54.9 / −18.9 / 11.7) are
      large vs a bu/acre outcome. If VPD is in kPa, sanity-check per-kPa marginal
      effects against the observed within-field VPD range before quoting any.
- [ ] **RCI×Dry interactions.** Non-monotone and individually noisy
      (RCI3×Dry=+5.36 vs RCI3.24×Dry=−5.32). Use a joint test or a linear-in-RCI
      interaction rather than reading individual factor-level coefficients before
      any drought-buffering claim.

---

## P3 — Prose verification (independent of nsoy/AWC)

- [ ] **RCI observed range on the 49-sequence wheat-inclusive sample** (Sec 2.4).
      Prior corn-soy-only range was 1.41–4.74; recompute and state from the
      current sample.
- [ ] **"18 distinct RCI values"** (Sec 2.4 and Sec 4.1) — verify the count on
      the current wheat-inclusive factor-RCI sample; may have changed.
- [ ] **RCI reference level** — several tables show RCI=1.41 with its own
      coefficient while the notes call it the omitted reference; Sec 2.4 says
      monoculture has RCI=0 but the sample range excludes 0. Clarify which level
      is actually omitted, consistently across `corn_rci`, `corn_rci_vpd`,
      `corn_rci_jp`.
- [ ] **`late_soy` sign convention** — the variable runs −6 to 0 (closer to 0 =
      more recent soybean). Add one clause at its definition (Sec 4.1) so the
      positive coefficient reads correctly (recency raises yield).
- [ ] **`tables/corn_lasso.tex`** — confirm whether AWC was in the LASSO control
      set (if so, re-run). Verify 29 selected / 28.4% / cv.glmnet agreement
      (26/89.7%). Convert sequence labels from raw CDL numeric codes to the
      C/S/W letter convention.
- [ ] **`corn_summary_stats.tex`** — if regenerated on the current sample,
      verify the bucket numbers (perfect 198.0/28.0; monoculture 187.2/33.2).
- [ ] **Soy-side tables** — earlier soy figures (e.g. `soy_lasso` N=419,886)
      predate the current re-run; re-pull against the current sample tiers.

---

## P4 — Bibliography

Findings from the 2026-09-02 reference audit (`bibliography.bib` vs all
`\cite*` calls, verified against published sources).

### Factual errors in `bibliography.bib`

- [ ] **`wu2025`** — `pages = {111--114}` and `number = {20}` are copy-paste
      artifacts from the `shi2013` entry above it. Correct record: Wu, Davis &
      Sohngen (2025), *Carbon Balance and Management* **20**, article **6**,
      doi:10.1186/s13021-025-00293-5. Drop `number = {20}`, keep `volume = {20}`
      / `issue = {1}`, set article number to `6`.
- [ ] **`Benami2026`** — `year = {2026}` is wrong. arXiv:2510.05108 was posted
      October **2025**. Change year to 2025.

### Mislabeled keys (render OK, but the key name misstates the year)

- [ ] **`Gentry2015`** → paper is **2013** (*Agronomy Journal* 105(2):295–303).
      `year` field is correct; rename key to `Gentry2013` and update the three
      `\citep{Gentry2015}` calls (`just_pope_rotations.tex` / `rotations_insurance.tex`
      / `rotations_losses.tex` line ~108/130/131).
- [ ] **`Tiemann2016`** → paper is **2015** (*Ecology Letters* 18(8):761–771).
      Rename key to `Tiemann2015` and update the `\citet{Tiemann2016}` call
      (line ~160/161 in both JP drafts).
- [ ] **`hennesy2006`** → author is Henne**ss**y; key drops an `s`. Uncited, so
      low priority — fix if the entry is ever used.

### Coverage gaps

- [ ] **`USDA2015`** (CropScape / CDL) is defined but never `\cite`d. The CDL is
      central to the data section (referenced only in prose). Add a `\cite` at
      first mention of the Cropland Data Layer, or drop the entry.
- [ ] Add the **Elhorst (2014)** BibTeX entry to `bibliography.bib`
      (spatial-panel reference).
- [ ] Unused entries carried only by `old.tex` (not in the build):
      `won2024`, `aglasan2024`, `yu2024`, `ortiz_bobea2021`. Also unused:
      `gammans2025`. Decide keep-for-later vs prune.

### Style nits (won't change output under `apalike`)

- [ ] `Antle2010` — author `{Antle, john M.}` (lowercase "john").
- [ ] `Cameron2015` — `month = {3u}` typo.
- [ ] `doi = {https://doi.org/10.1111/...}` in `won2024`, `aglasan2024`,
      `yu2024`, `Farmaha2016`, etc. — strip the URL prefix from the `doi` field.
- [ ] All three `.tex` set `\bibliographystyle{apalike}`, but the repo ships
      `mplainnat.bst` and the README cites it. Pick one.
- [ ] `ortiz_bobea2021` is `@inbook` but uses `journal =` instead of
      `booktitle =` for "Handbook of Agricultural Economics".
