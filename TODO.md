# TODO — Challenges of Incorporating Rotation Information into Crop Insurance Rates

Open issues only. Rebuilt 2026-09-02 after the reference audit and the
prose/table reconciliation pass on `just_pope_rotations.tex`.

**Manuscript of record: `just_pope_rotations.tex`.** `rotations_insurance.tex` and
`rotations_losses.tex` were moved to `old/` (commit `157fae0`); README and Makefile
updated to match (P0, 2026-09-02).

Ordered: build/spec blockers, then table/figure verification, then inference,
then prose verification, then bibliography.

---

## P0 — Blockers

- [x] **README pointed at the wrong manuscript.** Rewrote the stale sections of
      `README.md` for the live corn-only, wheat-in-history, B=999, field-FE paper
      and the `corn_analysis_full.R` pipeline: Overview, Repository structure,
      Data, Running the analysis, Tables/figures inventory, Building the paper,
      Econometric framework. `Makefile` target fixed (`main` → `just_pope_rotations`,
      + extra `pdflatex` pass). Remaining README nits are cosmetic (exact runtime,
      package list completeness). NOTE: `just_pope_rotations.tex` line ~665 has
      `\bibliography{bibliography.bib}`; the `.bib` extension can trip a local
      `bibtex` (looks for `bibliography.bib.bib`) — drop to `{bibliography}` if
      `make` fails at the bibtex step (fine on Overleaf). Ties into the P4
      apalike-vs-`mplainnat.bst` decision.
- [x] **`corn_jp_skew.tex` scale mismatch — deleted.** The ~2000× factor is s³
      (s ≈ 12.4 bu/acre): `corn_jp_skew.tex` reported the raw `resid_cube`
      coefficients under a "standardized" caption. Standardized stage-3 skewness
      with B=999 bootstrap SEs is already col (2) of `corn_jp_moments.tex` (the
      table the manuscript uses), so `corn_jp_skew.tex` was pure redundancy.
      `git rm`'d the file; removed `write_jp_skew_tex()` + its call from
      `code/write_jp_tables.R`; updated the appendix comment in
      `just_pope_rotations.tex`. (`corn_jp_var.tex` is likewise orphaned by the
      live manuscript but left in place — out of P0 scope.)

---

## P1 — Table/figure verification (regenerated artifacts vs prose)

Resolved in the reconciliation pass and not repeated here: `corn_rot`,
`corn_rci` (regenerated + §5.5/§7 prose rewritten), `corn_jp_moments` /
`corn_jp_var` / `corn_jp_skew` (bootstrap SEs, sane scale; §5.2 magnitudes and
skewness values updated), `zvector`, `corn_lasso` (29 sequences, C/S/W labels),
FE/clustering notes across all tables, RCI reference level (= 0), RCI value
count/range.

Verified 2026-09-02 (P1 pass):

- [x] `tables/corn_rot.tex` last three vs §5.1: S-C-C-C-C-C = −1.07 ✓;
      S-C-S-C-S-C = 3.637 (prose said "3.63" → fixed to 3.64); C-C-C-C-S-C =
      4.419 (prose said "4.41" → fixed to 4.42).
- [x] `tables/corn_rot_vpd.tex` vs §5.8: orderings correct (W-C-C-S-W-C,
      C-S-C-W-S-C top; C-S-C-C-C-C / S-C-C-C-C-C / S-S-C-C-C-C negative). VPD
      signs correct (vpd_6 +54.9, vpd_7 −18.9, vpd_8 +11.7). Drought bins
      −1.29 / −1.88, both n.s. ✓.
- [x] `tables/corn_rci_vpd.tex` vs §5.9: MISMATCH fixed. Significant "somewhat
      dry" interactions are RCI = 1.41 (1.42*) and 2.00 (1.61**) only; prose
      named "2.00, 2.24, 2.45" — 2.24 (1.58, n.s.) and 2.45 (1.87, n.s.) are
      not significant. Prose reworded; joint-test route still open (P2).
- [x] `figures/corn_rot_plot.png`, `corn_var_plot.png`: consistent with current
      coefficients — same commit (344d810) as `corn_rot.tex`, which is unchanged
      since; both show the 28 sequences (23 CS + 5 CSW) with ordering/magnitudes
      matching the table. Cosmetic: both have x/y axis labels swapped
      ("Crop sequence" on the numeric axis).
- [x] `figures/corn_jp_plot.png`: x-axis issue already resolved in the Sep-2
      regen — axis spans ≈ −1 to 9 bu/acre (ticks 0/4/8), matches `corn_rot.tex`.
      Stale inline TODO at the figure removed.
- [x] `figures/score_yield.png`: y-axis shows 14 distinct rotation-score values
      from −2.68 to +1.88 exactly as expected. ✓

Sequence-count reconciliation (28, not 48/49): the regenerated artifacts use 28
non-monoculture sequences (23 corn-soy + 5 corn-soy-wheat; 29 incl. monoculture,
= `corn_soy_patterns` in `rotation_setup_wa.R`). Stale "48" / "31 + 17" /
"49-sequence" prose corrected across `just_pope_rotations.tex` (lines 198, 200,
202, 157, 297, 336, 391–392, 415, 466, 565) and the `rotation_setup_wa.R` header
comments. (§7 "corn-soybean-only sample" wording is still P3.)

---

## P2 — Inference and specification

- [x] **Sequence-dummy variance + skewness bootstrap** — already done. Verified
      2026-09-02: `corn_jp_moments.tex` col (1) SEs = sqrt(diag(`jp_vcov_var.txt`))
      to 4 d.p. for all 28 sequences; col (2) SEs = sqrt(diag(`jp_vcov_skew.txt`))
      rescaled by one 1/s³ factor (s ≈ 12.4 bu/acre). Re-derived stars match the
      printed table exactly. Nothing weakens vs analytic: C-C-C-C-S-C variance
      z = −1.67 (p = 0.095, still *); skewness C-C-C-S-S-C p = 0.062,
      S-S-C-C-C-C p = 0.038, C-S-C-S-C-C p < 0.001; none significantly positive.
      §5.2 inline NOTE rewritten to record this (was wrongly claiming analytic SEs).
- [x] **Z-vector (three-feature) variance bootstrap** — DONE 2026-09-03.
      `code/zvector_bootstrap_var.R` run → `tables/zvector_boot_var.txt`
      (B = 999, field-pairs design). Bootstrap SEs: late_soy 0.286 (p=0.040, **),
      soy_gap 0.375 (p≈0, ***), soy_cons 1.354 (p≈0, ***) — all significant,
      vs analytic county+year 0.90 / 1.96 / 3.94 (only soy_gap marginal). Gap
      persists under field clustering too (script's analytic_se 0.86 / 1.68 /
      3.72), so it is the bootstrap's generated-regressor correction, not the
      clustering dimension. `zvector.tex` col (2) updated; `tables_combined.r`
      now re-applies the bootstrap SEs automatically after `etable()` and no
      longer writes duplicate FE rows. Prose updated: abstract (§ line 91),
      §Z-vector variance para, Discussion finding 2, Conclusion, Implications.
      CAVEAT for co-authors: the bootstrap SEs are ~4–5× tighter than either
      analytic clustering; prose flags this and leans hardest on soy_gap.
- [x] **VPD units** — confirmed kPa (from `corn_summary_stats.tex`: July `vpd_7`
      mean ≈ 1.1, SD ≈ 0.65 kPa). §5.8 prose + inline NOTE updated: `vpd_7`
      = −18.9/kPa ⇒ ≈ −12 bu/acre per raw SD (less after FE); `vpd_6` = +54.9/kPa
      flagged as not structurally identified (June weather-block collinearity),
      `vpd_8` likewise. Text now interprets only the July sign; no per-kPa
      marginals quoted for June/August.
- [x] **VPD double-counting (spec decision)** — RESOLVED 2026-09-03: continuous
      monthly VPD (`vpd_6/vpd_7/vpd_8`) dropped from `corn_rot_vpd` /
      `corn_rci_vpd`; the July drought effect is now identified through the bins
      alone. Comment added at `corn_analysis_full.R` ~line 419; both tables
      regenerated (vpd_6/7/8 rows gone, coefficients + within-R² shifted).
      §5.8 prose rewritten (removed the per-kPa continuous-VPD discussion; drought
      bins now 0.92 / 1.64 bu/acre, both n.s.); §5.9 penalty magnitude 2.7 → 2.5,
      RCI3×Dry now +6.67* noted; table notes blocks updated.
- [x] **RCI × Dry interactions** — §5.9 prose (fixed in P1) no longer makes a
      drought-buffering claim; it states only that the pattern is sign-flipping
      and imprecise. Inline TODO downgraded to a NOTE: a formal joint Wald test
      on the 9 RCI×Dry terms would upgrade "absence of evidence" to "evidence of
      absence" but needs a `corn_rci_vpd` re-run to save the sub-vcov. Not a
      blocker.

---

## P3 — Prose verification

Resolved 2026-09-02 (P3 pass):

- [x] **`late_soy` sign convention** (§4.2, Z-vector definition) — `late_soy` is
      `-min(soy_pos)` in `corn_analysis_full.R` (~line 576): a non-positive integer
      from −6 to 0, closer to 0 = more recent soybean, 0 for monoculture. Added a
      clause at the $Z_{it}$ definition stating the coding and range and noting a
      positive coefficient means greater soybean recency raises yield. The §5.7 /
      §7 interpretations were already correct in sign.
- [x] **`tables/corn_lasso.tex`** — AWC (`rootznaws_mean`) is NOT in the LASSO
      control set. `corn_selection.r` sources `rotation_setup_wa.R`, whose
      `all_controls` (→ `all_controls_cols`, passed as `controls=` to
      `lasso_select_sequences`) is weather + GRIDMET soil-moisture only
      (`pr/GDD/EDD/vpd/soil` for Jun–Aug); no gSSURGO soil vars. (`rootznaws_mean`
      appears in `all_controls` only in `just_pope.r` line 117, a different
      pipeline that does not write `corn_lasso.tex`.) No re-run needed. Counts
      already reconciled: table has 29 sequence rows; stale "28 of 102" comment in
      `corn_selection.r` line 187 is cosmetic.
- [x] **`corn_summary_stats.tex`** — committed table (`eedf565`): monoculture
      187.2 (33.2) ✓; perfect rotation 194.2 (31.1), not the 198.0 / 28.0 the
      prose claimed; transitioning 33.5 is the highest SD of the three, not
      monoculture's 33.2. §2.6 prose fixed to 194.2 / 31.1 and the "highest
      variability" clause reworded to cover both monoculture and transitioning.
      (Did not re-run the R pipeline; verified prose against the committed table.)
- [x] **"corn-soybean-only sample"** — §7 (Conclusion) reworded: "a sample
      dominated by corn-soybean rotations, with only five winter-wheat sequences
      entering on thin support" (5 of 28, matching §3.2 / §5.1), and the follow-on
      clause changed to "holds more broadly for rotations built around a third
      crop" so it no longer implies wheat is absent.

---

## P4 — Bibliography

From the 2026-09-02 reference audit (`bibliography.bib` vs all `\cite*` calls,
verified against published sources). `Elhorst2014` has since been added and is
resolved.

Applied 2026-09-02 unless noted.

### Factual errors

- [x] **`wu2025`** — dropped the copy-pasted `pages = {111--114}`; now
      `volume = {20}`, `number = {1}`, `pages = {6}` (article number), added bare
      `doi = {10.1186/s13021-025-00293-5}`.
- [x] **`Benami2026`** — `year` → `2025`, added `month = {10}` (arXiv:2510.05108,
      Oct 2025). Key left as `Benami2026`; `\cite` calls render "(2025)" now.

### Mislabeled keys

- [x] **`Gentry2015` → `Gentry2013`** — key renamed; both `\citep` calls in
      `just_pope_rotations.tex` updated.
- [x] **`Tiemann2016` → `Tiemann2015`** — key renamed; `\citet` call updated.
- [x] **`hennesy2006` → `hennessy2006`** — key fixed (uncited, no call to update).

### Coverage / style

- [x] **`USDA2015`** — `\citep[CDL;][]{USDA2015}` added at first CDL mention (§3.1).
- [x] `Antle2010` — `john` → `John`.
- [x] `Cameron2015` — `month = {3u}` → `{3}`.
- [x] URL-prefixed `doi` fields — stripped `https://doi.org/` and
      `https://dx.doi.org/` from all 22 entries.
- [x] `ortiz_bobea2021` — `journal =` → `booktitle =`, added `publisher = {Elsevier}`.
- [x] **Unused entries** — decision: keep all (`aglasan2024`, `won2024`,
      `yu2024`, `ortiz_bobea2021`, `gammans2025`, `lobell2015`, `hennessy2006`).
      Harmless; several are used by `old.tex`.
- [x] **`apalike` vs `mplainnat.bst`** — decision: author-year, keep `apalike`.
      Fixed the preamble inconsistency: natbib options
      `[square,sort,comma,numbers]` → `[round,sort,comma,authoryear]`. `mplainnat.bst`
      now unused (left in repo). In-text cites render as "Author (year)".
