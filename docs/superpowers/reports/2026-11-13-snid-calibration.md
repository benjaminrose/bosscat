# SNID-SAGE calibration on the 105 BOSS ToO spectra (2026-11-13)

SNID-SAGE 1.2.5, templates v1.2.5 (optical grid), all 105 spectra in
`data/spectra/`. Reference: NGSF (`data/manifest.csv`: `z_fit`, `sn_family`)
and TNS (`tns_type`, 19 spectra). Outputs in `runs/` (not committed).

## Headline

1. **The high-z calibration criterion was not met: 0 of 5 recovered.**
   `zmax = 0.7` changes nothing on these data. Three of the five fail in SNID
   (S/N 1.2 to 1.7), one is fitted at the wrong (low) redshift, and one is
   probably an NGSF error (see below). Thresholds therefore could not be
   calibrated on the high-z set.
2. **`zmax` 0.7 vs 0.3: no new disagreements at rebin 1.** Every broad type is
   identical and no consensus redshift moves by more than 0.01. `zmax = 0.7`
   is kept (harmless, needed for future high-z data) but is unproven.
3. **Rebinning dropped (`rebin_factor = 1`).** Factors 4 and 8 raise type
   agreement, but not significantly once three factors are tried, and they do
   not improve redshifts (more catastrophic outliers). Details below; this is
   the closest call in the report.
4. **Flag thresholds kept at SNID-SAGE defaults**, which makes
   `snid_flag_low_quality` equivalent to `snid_success == False` (SNID-SAGE
   already discards matches below `lapmin = 0.3` and HσLAP-CCC 1.5).
5. **`snid_consensus_age`, `snid_consensus_age_err` and `snid_r` are always
   0.0** (every success in every run). They are library placeholders.

## Runs

Command: `MPLBACKEND=Agg uv run bosscat data/spectra/*.fits --zmax Z
--rebin-factor K -o runs/zZ_rebinK.csv`, then `python -m bosscat.evaluate`.
Three runs in parallel on 12 cores; wall time per run 638 to 666 s
(~6.3 s per spectrum per process).

| config | wall (s) | success | NGSF agree | TNS agree | med Δz | MAD Δz | \|Δz\|<0.02 | \|Δz\|>0.05 | low_snr | nonsn |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| z0.3 rebin1 | 638 | 92 | 61/88 (0.69) | 13/18 (0.72) | 0.0020 | 0.0030 | 81/92 | 7 | 34 | 8 |
| z0.7 rebin1 | 666 | 92 | 61/88 (0.69) | 13/18 (0.72) | 0.0020 | 0.0030 | 82/92 | 7 | 34 | 8 |
| z0.7 rebin2 | 661 | 90 | 59/86 (0.69) | 13/18 (0.72) | 0.0019 | 0.0030 | 80/90 | 7 | 20 | 9 |
| z0.7 rebin4 | 666 | 95 | 66/90 (0.73) | 15/19 (0.79) | 0.0019 | 0.0029 | 82/95 | 10 | 11 | 6 |
| z0.7 rebin8 | 666 | 95 | 69/90 (0.77) | 16/19 (0.84) | 0.0019 | 0.0032 | 84/95 | 8 | 5 | 5 |

Δz = `snid_consensus_redshift - z_fit`. Agreement is broad type
(`bosscat.evaluate.broad_type`); non-SN SNID types count as disagreement.
`snid_flag_low_quality` and `snid_flag_zedge` were 0 among successes in
every run. Maximum consensus redshift in any rebin-1 run: 0.251.

## Rebinning

Grid spacing: BOSS pixels are Δlog10 λ = 1.0e-4 (69 km/s). SNID-SAGE
resamples to 1024 bins over 2500 to 10000 Å, Δln λ = 1.354e-3 (406 km/s),
i.e. 5.9 BOSS pixels per SNID bin (checked on the processed spectrum's
`log_wave`). Factors 2 and 4 stay finer than the SNID grid; factor 8
(552 km/s) is coarser. SNID-SAGE is given flux only, no `ivar`, so the
gain from rebinning is the inverse-variance weighting it otherwise lacks.

Paired comparison against rebin 1 (z0.7, all 105 spectra; SNID failure
counts as disagreement; two-sided sign test on gained vs lost):

| factor | type changed vs rebin 1 | NGSF gained/lost (p) | TNS gained/lost (p) | \|Δz\|>0.05 |
| --- | --- | --- | --- | --- |
| 2 | 4/105 | 0/2 (0.50) | 0/0 | 7 |
| 4 | 8/105 | 6/1 (0.13) | 2/0 (0.50) | 10 |
| 8 | 15/105 | 10/2 (0.04) | 3/0 (0.25) | 8 |

Decision: keep `rebin_factor = 1`. The rule was to ship a factor only if
agreement is more stable than at factor 1. Type agreement rises at 4 and 8,
but the best p (0.04 at factor 8) is not significant after trying three
factors (Bonferroni 0.12), redshift residuals do not improve, catastrophic
redshift errors rise (7 to 10 at factor 4), and results move for up to
15 objects. Two of the new high-z redshifts are wrong: AT2026fxx (NGSF Ic,
z 0.161) goes to Ic at z 0.36 (factor 4) and 0.34 (factor 8), which only
`zmax > 0.3` allows. Gains come mostly from low-S/N spectra that fail at
factor 1, plus two Galaxy fits that become II (smoothing of narrow host
lines). **Factor 4 is a credible alternative for the user to weigh**: it
loses no information on SNID's grid and recovers 3 more successful fits.

## Flag calibration (5 objects with NGSF `z_fit > 0.3`)

Redshift tolerance: |Δz| < 0.02. Among all rebin-1 successes the median
|Δz| is 0.004, the 80th percentile 0.009, and the median SNID consensus
error 0.003, so 0.02 is ~4.5σ of the core; beyond it the residuals form a
tail (0.032 to 0.49). The choice does not matter here: every high-z miss is
off by > 0.4.

| name | NGSF type, z_fit | TNS | SNID (z0.7, rebin 1) | lap | HσLAP-CCC | S/N | flags |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AT2026fvw | SLSN-I, 0.644 | | fail | | | 1.2 | low_quality, low_snr |
| AT2025cee | Ia-02cx, 0.667 | | fail | | | 1.2 | low_quality, low_snr |
| AT2025muz | II, 0.438 | | fail | | | 1.7 | low_quality, low_snr |
| SN2025ewp | II-flash, 0.528 | SN II | II, 0.041 | 0.69 | 1.57 | 10.1 | none |
| SN2025kpa (MJD 60822) | SLSN-I, 0.493 | SN Ia-91T | Ia, 0.076 | 0.66 | 1.56 | 12.0 | none |

- SN2025kpa was observed the night before (MJD 60821), where NGSF gives
  super-Chandra Ia at z 0.075 and SNID gives Ia at 0.074. With TNS also Ia,
  the NGSF z 0.493 SLSN-I fit at 60822 is probably wrong, so SNID's answer is
  probably right. Restricting SNID to z 0.3 to 0.7 finds the SLSN-I solution
  at z 0.494 (HσLAP-CCC 2.37 at factor 4, higher than the low-z Ia at 1.96):
  a real type/redshift degeneracy that SNID's clustering resolves to low z.
- SN2025ewp: right type, wrong redshift. Restricting to z ≥ 0.3 gives no
  match at all (factors 1 and 4).
- The three failures stay failures at factor 4, except AT2025muz (II at
  z 0.264, still wrong).

Separation on the full sample (rebin 1, 84 SN-type successes, "wrong" =
|Δz| > 0.02, n = 10):

| flag rule | flags wrong | flags right |
| --- | --- | --- |
| HσLAP-CCC < 1.5 (default) | 0/10 | 0/74 |
| HσLAP-CCC < 2.0 | 8/10 | 45/74 |
| HσLAP-CCC < 2.5 | 9/10 | 61/74 |
| S/N < 3 | 5/10 | 29/74 |

Type-agreement failures are no better separated: good and bad fits have
indistinguishable lap and HσLAP-CCC distributions (bad fits include
high-S/N Galaxy matches with HσLAP-CCC up to 11).

Decisions (values unchanged from provisional; now justified):

- `min_lap = 0.3`, `min_hsigma_lap_ccc = 1.5`: SNID-SAGE defaults, since the
  data cannot discriminate. Consequence: `snid_flag_low_quality` marks only
  failed fits. Option for the user: 2.0 catches both high-z misses and 8/10
  wrong redshifts but flags 61% of correct ones (wrong-z rate 15% when
  flagged vs 6% when not). SNID-SAGE's own cluster-quality scale calls < 2.5
  "Very Low", but that is a penalized top-5 cluster score, not this metric.
- `min_snr = 3.0`: all 13 SNID failures at rebin 1 have S/N 0.96 to 2.82.
  `snid_snr_median` is measured after rebinning, so this threshold must be
  revisited if a rebin factor is adopted (S/N scales roughly as √factor).
- `zedge_margin = 0.02`: untestable (no rebin-1 consensus z above 0.252); kept at
  ~6× the median consensus redshift error.

Uncertainty: n = 5, of which 2 produced fits and 1 of those is likely an
NGSF error. Nothing above is a calibration in the statistical sense.

## consensus_age and r

In all five runs (464 successful fits) `snid_consensus_age`,
`snid_consensus_age_err` and `snid_r` are exactly 0.0. In SNID-SAGE 1.2.5,
`result.r` comes from the vectorized peak finder, which hard-codes
`r_value = 0.0`, and `consensus_age` reads a `getzt` key (`age_enhanced`)
that only the no-match fallback dicts define. Recommendation: drop the three
columns (or keep them only behind a check that they are non-zero).
`result.template_age` (age of the best template) is populated and could
replace them as `snid_template_age`; not done, needs sign-off.

## Air vs vacuum (measured, not enabled)

BOSS wavelengths are vacuum. Template convention: **not determinable**. The
template HDF5 attributes, `template_index.json` and `templates_meta.json`
state nothing; the package docs are silent; the line-overlay code defaults
to air wavelengths (indirect only). Centroiding narrow lines in the two
emission-line galaxy templates (kcSB1, kcSB2) sits 0.0 ± 0.15 SNID bins from
vacuum and 0.2 bins from air, which is below the measurement scatter and says
nothing about the SN templates.

Effect of converting BOSS vacuum to air (IAU formula) before SNID:

| spectrum | Δ consensus z | type | HσLAP-CCC |
| --- | --- | --- | --- |
| AT2026jea (Ia, z 0.079) | -0.00030 | unchanged | 1.96 → 2.03 |
| AT2026fav (Ia, z 0.079) | -0.00028 | unchanged | 5.62 → 5.61 |

This is the expected 2.8e-4 (1+z) shift, ~10% of the consensus redshift
error. Negligible for classification.

## Galactic dereddening (measured, not enabled)

Cardelli, Clayton and Mathis (1989) law, R_V = 3.1, with the BOSS `ebv`
column (the shared column is `ebv`, not `boss_ebv`; sample median 0.029,
max 0.150). Two highest-`ebv` successful fits:

| spectrum | E(B-V) | Δ consensus z | type | best template, HσLAP-CCC |
| --- | --- | --- | --- | --- |
| AT2026jea | 0.150 | +0.00001 | Ia-03fg, unchanged | sn2005M, 1.96 → 2.11 |
| AT2026fav | 0.125 | 0.00000 | Ia-norm, unchanged | sn2007if 5.62 → sn2019so 2.05 |

Redshift and type are unchanged, but on AT2026fav the best template changed
and HσLAP-CCC dropped from 5.6 to 2.0. SNID removes a continuum, so the
residual sensitivity is surprising; with n = 2 it is a reason to test on more
spectra before enabling, not a result.

## Caveats

- NGSF is a reference, not truth; it models host light and SNID does not, so
  type disagreements (especially Galaxy fits, 8 at rebin 1) are expected.
- Only 18 to 19 spectra have a TNS type; TNS agreement differences of 2 to 3
  objects are within noise.
- The calibration set is 5 spectra, 3 of them at S/N < 2.
- Spectra are low S/N (median ~3.6); 13 of 105 fail in SNID at rebin 1.
- SNID clips the BOSS spectrum to its 2500 to 10000 Å grid (the last ~400 Å
  of BOSS coverage is lost).
