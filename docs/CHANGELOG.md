# Changelog

## February 2026 — Phase II Complete

### New Notebooks

- **Notebook 05** (`05_phase2_density_coupling.ipynb`): Population split analysis, density-dependent coupling regression, mass-to-light sensitivity test, and M33 linear vs. tapered re-analysis.
- **Notebook 06** (`06_model_gallery.ipynb`): Head-to-head model gallery across 7 diverse galaxies with surface brightness predictor.

### New Results

- **Two-population discovery**: 84% of galaxies prefer the tapered model (interior solutions), 16% prefer the pure linear model (boundary solutions). The populations differ significantly in luminosity, surface brightness, and flat velocity ($p < 0.001$).
- **$k$–$\Sigma_0$ independence**: No significant correlation ($R^2 = 0.009$, $p = 0.34$) between the coupling constant and surface brightness — supports universality.
- **$\Upsilon$ sensitivity quantified**: $\omega$ varies 28% (LSB) to 134% (HSB) across plausible mass-to-light ratios. Model is most robust for gas-dominated systems.
- **Model gallery (80% accuracy)**: Surface brightness predictor correctly identifies the BIC-preferred model for 4/5 testable galaxies. The NGC 2841 failure suggests the tapered model may be more broadly applicable than the population split implies.

### New Figures

- `split_populations.png` — Population boxplots
- `k_vs_surface_brightness.png` — Density coupling test
- `upsilon_sensitivity.png` — Mass-to-light robustness
- `M33_linear_reanalysis.png` — Clean head-to-head M33 comparison
- `gallery_*.png` — Per-galaxy model comparisons (DDO154, DDO161, NGC0300, NGC2841, NGC3198, NGC7331, M33)

### New Tables

- `SPARC_unified_coupling_results.csv` — Re-fit of 19 boundary galaxies
- `upsilon_sensitivity.csv` — Three-galaxy $\Upsilon$ sensitivity grid
- `M33_linear_reanalysis.csv` — M33 linear vs. tapered comparison
- `model_gallery_validation.csv` — Gallery prediction vs. outcome

### Documentation Updates

- `README.md` — Complete rewrite with key results, gallery table, and updated figures
- `METHODOLOGY.md` — Added Phase II methods (sensitivity, population analysis, gallery validation, method versioning)
- `docs/RESULTS.md` — Full Phase I & II report with all figures and discussion questions
- `docs/CHANGELOG.md` — This file

---

## February 2026 — Phase I Complete

### Notebooks 01–04

- M33 calibration and pipeline validation (Notebook 01)
- Linear vs. quadrature mechanism test (Notebook 02)
- Tapered model development (Notebook 03)
- 118-galaxy SPARC batch analysis (Notebook 04)

### Key Results

- Linear mechanism confirmed ($\Delta$BIC = $-4,559$)
- Rational taper reduces M33 RMSE by 69%
- Candidate universal coupling constant $k \approx 2.8$
- 100% convergence rate across 118 galaxies
- Median RMSE = 5.9 km/s, median $\chi^2_\nu$ = 1.47
