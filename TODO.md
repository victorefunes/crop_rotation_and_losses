# TODO — Challenges of Incorporating Rotation Information into Crop Insurance Rates

Consolidated from the paper's inline TODOs plus the changes made this session
(nsoy drop, AWC drop, three-feature Z-vector, score update). Ordered by priority:
blockers first, then table/figure regeneration, then remaining prose/verification.

---

## P0 — Blockers (prose now contradicts the committed artifacts)

- [x] **Regenerate `tables/zvector.tex` as the THREE-feature model.**
      DONE. `tables/zvector.tex` now holds the three-feature model (late_soy 0.4470***,
      soy_gap 0.6928***, soy_cons 1.210**; variance soy_gap −3.748*, p<0.10), N = 1,785,861,
      field + year FE, clustered at county + year. Matches prose Sec 5.4 / 6.1 exactly.
      Redundant `tables/corn_zvector.tex` removed. `\input` unchanged (line 468). Notes
      block under the table corrected: FE is field (tile_field_ID) + year, clustering is
      county + year (was mis-stated as county FE / field clustering).
- [x] **Regenerate `figures/score_yield.png` from the three-feature `score_plot_df`.**
      DONE. Regenerated figure's y-axis (rotation score) runs −2.68 → +1.88 with exactly
      14 distinct values (1.88, 1.18, 1.01, 0.74, 0.56, 0.49, 0.04, −0.33, −0.4, −0.89,
      −1.34, −1.79, −2.23, −2.68). Old −144/+229 span is gone; matches prose (lines 120, 730).

---

## P1 — Regenerate all tables/figures on the AWC-free control vector

AWC (`rootznaws_mean`) was dropped from the R control vector this session. Every
table whose spec included it must be re-exported, and any bushel figures quoted in
prose re-checked against the new values. Re-run the table/figure exports in
`corn_analysis_full.R`, then tick each and reconcile the prose numbers listed.

- [~] `tables/corn_rot.tex` (sequence-level mean). Regenerated (Sep 1). Verified 3/6:
      S-S-C-S-S-C=5.944, C-S-C-S-S-C=4.825, C-S-C-C-C-C=−0.768 match prose (5.94/4.82/−0.77).
      STILL VERIFY: S-C-C-C-C-C=−1.08, S-C-S-C-S-C=3.63, C-C-C-C-S-C=4.41.
- [~] `tables/corn_rci.tex` (factor RCI mean/variance).
      [x] Standing issue FIXED in the regen: cols (1)/(2) no longer byte-identical;
          headers now read "No controls" / "Weather and soil controls".
      [ ] Sec 5.5 prose STILL STALE: cites RCI 1.73 negative, 2.00–2.65 ≈ reference,
          positives at 2.74/4.24/4.74. Regenerated table has levels 1.41, 1.73, 2, 2.24,
          2.45, 2.65, 3, 3.24, 3.46 — no 2.74/4.24/4.74, and 2.24/2.45/2.65 are all
          significantly positive, not indistinguishable from reference. Rewrite lines
          525 and 837, clear inline TODOs at 492/502/525.
- [ ] `tables/corn_rot_vpd.tex` — verify Sec 5.8 (dry-July −1.29/−1.88; vpd terms).
- [ ] `tables/corn_rci_vpd.tex` — verify Sec 5.9 (positive interactions at RCI 2.00/2.24/2.45).
- [ ] `tables/corn_jp_moments.tex` (Appendix) — verify Sec 5.2 C-C-C-C-S-C variance;
      resolve the "11.3 units" vs current −44,885.8 scale mismatch and restate units.
- [ ] `tables/corn_jp_skew.tex` (Appendix) — verify Sec 5.2 skewness values
      (C-S-C-S-C-C −0.31, S-S-C-C-C-C −0.33, C-C-C-S-S-C −0.37).
- [ ] `tables/corn_rci_jp.tex` (Appendix) — verify; confirm omitted RCI reference level.
- [ ] `figures/corn_rot_plot.png`, `corn_var_plot.png`, `corn_jp_plot.png` — regenerate
      from the updated coefficients. Resolve corn_jp_plot x-axis unit/scale (labeled
      bu/acre but spans ≈ −200/500).
- [ ] `figures/corn_rci_plot.png` — regenerate alongside corn_rci.tex.

---

## P2 — Inference and specification consistency

- [ ] **Bootstrap the three-feature variance/skew stages.** Paper's inference of record
      is the three-stage bootstrap (Sec 4.3), not analytic SEs. Re-confirm the
      marginal soy_gap variance effect (p<0.10 analytic — may weaken) and the skewness
      significance stars (Sec 5.2 note flags this) against bootstrap draws.
- [ ] **FE inconsistency in VPD tables.** `corn_rot_vpd` / `corn_rci_vpd` use field
      (tile_field_ID) FE; the main mean model (eq. 1) uses county (FIPS) FE. Either
      re-estimate under a common FE structure or state the difference explicitly so
      sequence coefficients are comparable across tables.
- [ ] **VPD double-counting check.** VPD models include both continuous July VPD (vpd_7)
      AND the July dry-season bins. Confirm intended; if bins capture nonlinearity,
      consider dropping continuous vpd_7 or state it holds the within-bin slope.
- [ ] **VPD coefficient scale/units.** vpd_6/7/8 (≈54.9/−18.9/11.7) are large vs a
      bu/acre outcome; if VPD is in kPa, sanity-check per-kPa marginal effects against
      the observed within-field VPD range before quoting any in text.
- [ ] **RCI×Dry interactions.** Non-monotone and individually noisy
      (RCI3×Dry=+5.36 vs RCI3.24×Dry=−5.32). Use a joint test or a linear-in-RCI
      interaction rather than reading individual factor-level coefficients before
      any drought-buffering claim.

---

## P3 — Remaining prose verification (independent of nsoy/AWC)

- [ ] **RCI observed range on the 49-sequence wheat-inclusive sample** (Sec 2.4).
      Prior corn-soy-only range was 1.41–4.74; recompute and state from current sample.
- [ ] **"18 distinct RCI values"** (Sec 2.4 and Sec 4.1) — verify the count on the
      current wheat-inclusive factor-RCI sample; may have changed.
- [ ] **RCI reference level** — several tables show RCI=1.41 with its own coefficient
      while the notes call it the omitted reference; and Sec 2.4 says monoculture has
      RCI=0 but the sample range excludes 0. Clarify which level is actually omitted,
      consistently across corn_rci, corn_rci_vpd, corn_rci_jp.
- [ ] **`late_soy` sign convention** — the variable runs −6 to 0 (higher/closer-to-0 =
      more recent soybean). Add one clause at its definition (Sec 4.1) so the positive
      coefficient reads correctly (recency raises yield).
- [ ] **`tables/corn_lasso.tex`** — confirm whether AWC was in the LASSO control set
      (if so, re-run). Verify 29 selected / 28.4% / cv.glmnet agreement (26/89.7%).
      Convert sequence labels from raw CDL numeric codes to the C/S/W letter convention.
- [ ] **`corn_summary_stats.tex`** — if regenerated on the current sample, verify the
      bucket numbers (perfect 198.0/28.0; monoculture 187.2/33.2).
- [ ] **Soy-side tables** — earlier soy figures (e.g. soy_lasso N=419,886) predate the
      current re-run; re-pull against the current sample tiers.

---

## P4 — Bibliography / misc

- [ ] Add the Elhorst (2014) BibTeX entry to `bibliography.bib` (spatial-panel reference).
- [x] Sample-size prose reconciled. The current re-run puts EVERY field-level corn table
      (corn_rot, corn_rci, corn_jp_moments, corn_jp_var, corn_jp_skew, zvector,
      corn_rot_vpd, corn_rci_vpd, corn_rci_jp) on one common sample of 1,785,861;
      only corn_lasso differs (1,799,074). Sec 2.4 / Sec 2.7 prose (lines ~209, 231) and
      the comment block rewritten: 1,796,806 is now described as the pre-VPD-filter
      qualifying count only, not an estimation-sample tier. The old claim that JP-moment
      and Z-vector models use 1,796,806 was removed (it was false against the regen).

---

## Done this session (for reference)

- [x] P0: `tables/zvector.tex` rewritten as three-feature model; `corn_zvector.tex` removed;
      table Notes FE/clustering corrected; `score_yield.png` verified against prose.
- [x] P4: sample-size prose collapsed to one common corn sample (1,785,861); stale
      three-tier story and false "JP/Z-vector use 1,796,806" claim removed.
- [x] Dropped `nsoy` from Z-vector (collinear with soy_gap, r=0.65; sign unstable).
- [x] Rewrote Z-vector / score as three features (Sec 4.1, 4.4, 5.4, 5.7, 6.1, 6.2, conclusion).
- [x] Softened variance claims to "marginal/directional" (Option A) throughout.
- [x] Updated score numbers everywhere (−2.68/+1.88, 14 distinct, mean 0.16, median 0.56, perfect +0.49).
- [x] Dropped AWC from Sec 3.4 prose and the Z-vector table note (unidentified: infinite SE).
- [x] Removed the nsoy sign-flip and stale score TODO blocks; replaced with action notes.
