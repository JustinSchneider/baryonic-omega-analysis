# Phase I Analysis: The "Tapered Linear" Scaling Law

Empirical validation of the Omega model and discovery of a Baryonic Coupling Constant ($k \approx 2.5$).

**Authors:** Justin Schneider, David C. Flynn, Jim Cannaliato
**Date:** February 2026

---

## 1. Executive Summary

This analysis rigorously tested the empirical "Linear Omega" model ($V = V_{bary} + \omega R$) against the SPARC galaxy database (Lelli et al. 2016).

**Key Findings:**

1.  **Pipeline Verification:** The codebase produces internally consistent baryonic velocity profiles from Corbelli et al. (2014) surface density data, confirming reproducibility of the ingestion-to-analysis pipeline (Section 2).
2.  **Mechanism Identification:** Statistical testing (BIC) on M33 strongly favors a **Linear (Kinematic)** addition over a **Quadrature (Dynamic/Halo)** addition ($\Delta\text{BIC} = 4559$). This suggests the anomaly is a velocity-dependent effect (e.g., frame dragging/vorticity), not a mass-dependent potential well.
3.  **The "Tapered" Correction:** The pure linear model diverges at large radii. We introduced a **Rational Taper** ($V_{boost} = \frac{\omega R}{1 + R/R_t}$) which reduced M33 residual error by **69%** (RMSE 31.1 km/s $\to$ 9.5 km/s).
4.  **Universal Scaling Law:** Batch analysis of 118 SPARC galaxies reveals a candidate coupling constant. Among the 89 well-fit galaxies ($\chi^2_{red} < 5$), the transition radius scales with the disk scale length with a median coupling factor of $k \approx 2.5$.

---

## 2. Methodology & Validation

**Objective:** Verify pipeline correctness before testing new physics.

- **Target:** M33, using surface density data extracted from Corbelli et al. (2014), Table 1 (58 radial bins, $R = 0.24 \text{–} 22.72$ kpc).
- **Process:** Implemented a Casertano (1983) thin-disk solver to convert HI gas and stellar surface densities ($\Sigma_{HI}$, $\Sigma_\star$) into velocity contributions ($V_{gas}$, $V_{disk}$), then combined them into $V_{bary}$.
- **Round-trip check:** The $V_{bary}$ computed during notebook analysis matches the $V_{bary}$ computed during database ingestion with 0.00% deviation across all 48 radial bins above $R > 2$ kpc. This confirms that the pipeline is deterministic and internally self-consistent.

**Important caveat:** This validation is a **self-consistency test** — both the "test" and "reference" $V_{bary}$ values are computed by our own Casertano solver from the same Corbelli (2014) surface density inputs. It verifies that no numerical errors, rounding artifacts, or parameter mismatches were introduced between the ingestion and analysis stages. It does **not** constitute an independent comparison against Corbelli's published velocity decomposition, which would require their original velocity component data (not available in their Table 1, which provides surface densities only). An independent cross-check against Corbelli's Fig. 7 mass models remains a desirable future validation step.

- **Implication:** Any deviations in the model fits below are attributable to the physical models tested, not to software bugs or data handling errors.

---

## 3. Model Comparison: Force vs. Flow

**Objective:** Determine if $\omega$ acts like a "Dark Matter Halo" (Force) or a "Background Current" (Flow).

We fit two functional forms to the M33 data (58 data points):

- **Model A (Linear):** $V_{obs} = V_{bary} + \omega R$ — Kinematic/Coriolis-like.
- **Model B (Quadrature):** $V_{obs} = \sqrt{V_{bary}^2 + (\omega R)^2}$ — Dynamic/Centrifugal-like, analogous to standard DM halo modeling.

**Results (M33):**

| Metric | Linear (A) | Quadrature (B) |
|--------|-----------|----------------|
| $\omega$ (km/s/kpc) | 6.97 | 9.84 |
| RMSE (km/s) | 31.1 | 45.3 |
| BIC | 4160 | 8718 |
| $\chi^2_{red}$ | 72.9 | 152.9 |

$\Delta\text{BIC} = 4559$ in favor of the Linear model — "very strong" evidence per Kass & Raftery (1995).

**Conclusion:** The data strongly rejects the standard "Quadrature addition" used in Dark Matter halo modeling. The effect adds linearly to velocity, supporting a **kinematic origin** (e.g., rotation/shear). Both models yield poor absolute fits ($\chi^2_{red} \gg 1$), indicating the pure linear form is incomplete — motivating the tapered model in Section 4.

---

## 4. The "Tapered Linear" Solution

**Problem:** While Linear is the correct _mechanism_, the $V \propto R$ term grows unboundedly, conflicting with observed flat rotation curves at large radii.

**Hypothesis:** The kinematic coupling "saturates" or decouples from the disk at a specific characteristic radius ($R_t$).

**New Model Tested:**
$$V_{total} = V_{bary} + \frac{\omega R}{1 + \frac{R}{R_t}}$$

This rational function interpolates between pure linear growth ($\omega R$ for $R \ll R_t$) and a constant asymptotic velocity ($\omega R_t$ for $R \gg R_t$).

**Results (M33):**

| Metric | A: Pure Linear | B: Rational Taper | C: Tanh Taper |
|--------|---------------|-------------------|---------------|
| $\omega$ (km/s/kpc) | 6.97 | 42.97 | 22.69 (effective) |
| $R_t$ (kpc) | — | 1.98 | 3.09 |
| RMSE (km/s) | 31.1 | **9.5** | 13.5 |
| $\chi^2_{red}$ | 72.9 | **4.6** | 9.7 |
| BIC | 4160 | **264** | 551 |

- **Fit Quality:** RMSE dropped from **31.1 km/s** (Pure Linear) to **9.5 km/s** (Rational Taper) — a 69% reduction.
- **Physical Alignment:** The rational taper model tracks the flat outer rotation curve, eliminating the "check mark" divergence of the pure linear model.
- **Model selection:** BIC strongly favors the rational taper over both the pure linear ($\Delta\text{BIC} = 3896$) and tanh taper ($\Delta\text{BIC} = 287$) forms.

![M33 Tapered Linear Models — Comparison of pure linear, rational taper, and tanh taper fits to the M33 rotation curve (Corbelli 2014 data). Upper panel shows circular velocity vs. radius with model predictions extrapolated beyond the data limit (22.7 kpc). Lower panel shows residuals with error bars. The rational taper (green) provides the best fit with RMSE = 9.5 km/s.](../results/figures/M33_tapered_models.png)

*Figure 1: M33 rotation curve fits. The pure linear model (red) diverges beyond ~10 kpc. The rational taper (green, $R_t = 1.98$ kpc) saturates to a constant velocity boost, closely tracking the observed flat rotation. The tanh taper (blue) performs intermediately. Data points are from Corbelli et al. (2014).*

---

## 5. Discovery: The Universal Coupling Constant ($k$)

**Objective:** Determine if the transition radius ($R_t$) is random or tied to galaxy structure.

**Method:** We ran the $R_t = k \cdot R_d$ parameterization (see Methodology, Section 4) on 118 SPARC galaxies, fitting $\omega$ and $k$ simultaneously while fixing $R_d$ from the SPARC photometric catalog.

**The Scaling Law:**
$$R_t = k \cdot R_d$$

**Findings:**

- **Sample:** 118 galaxies successfully converged (all fits achieved numerical convergence).
- **Quality filter:** 89 galaxies (75%) achieved $\chi^2_{red} < 5$, our threshold for acceptable fit quality.
- **The Constant (quality-filtered):** The coupling factor $k$ has a median of **2.55** among the 89 well-fit galaxies. The full-sample median is 2.81 (mean = 6.14), with the distribution skewed by 22 galaxies (19%) whose $k$ values hit the parameter bounds ($k = 0.1$ or $k = 20$), indicating the taper is unconstrained for those systems.

![Coupling Constant Universality Check — Three-panel plot showing k vs. flat velocity, Hubble morphological type, and disk scale length for 118 SPARC galaxies, color-coded by quality flag. The horizontal dashed line marks the median k = 2.81 (full sample).](../results/figures/k_universality.png)

*Figure 2: Coupling constant $k = R_t / R_d$ as a function of galaxy properties. The constant shows no strong systematic trend with flat velocity, morphological type, or disk scale length, consistent with a universal value. Points are color-coded by SPARC quality flag ($Q$: 1=high, 2=medium, 3=low). The full-sample median is $k = 2.81$ (dashed line); restricting to well-fit galaxies ($\chi^2_{red} < 5$) yields $k = 2.55$.*

**Caveats:**

- 22 of 118 galaxies (19%) have $k$ values at the parameter bounds (3 at $k = 0.1$, 19 at $k = 20$), meaning the transition scale is effectively unconstrained for those systems. These are predominantly galaxies where the linear rise extends to the outermost data point (no observed turnover), or where $V_{bary}$ already overpredicts $V_{obs}$.
- The $k$ distribution is right-skewed (mean 6.14 vs. median 2.81 for the full sample), driven by the boundary pile-up at $k = 20$.
- The "universality" of $k$ is stronger when restricted to the well-fit subsample, which may introduce selection effects favoring galaxies whose rotation curves happen to exhibit a clear transition.

**Physical Interpretation (tentative):**
Among galaxies where the taper model is well-constrained, the transition from linear velocity growth to flat rotation occurs at approximately **$2\text{–}3$ disk scale lengths**. If confirmed, this would imply that the "dark" velocity influence is geometrically linked to the baryonic disk extent.

---

## 6. Scaling Relations

The batch analysis reveals an apparent inverse trend between $\omega$ and galaxy luminosity: lower-luminosity galaxies tend to have higher $\omega$ values. This is qualitatively consistent with a picture where $\omega$ compensates for the mass deficit more aggressively in less massive systems. However, the scatter is substantial, and a rigorous correlation analysis (controlling for distance, inclination, and quality flag) is deferred to Phase II.

---

## 7. Recommendations

1.  **Strengthen the $k$ measurement:** The 19% boundary-hit rate suggests the parameter space should be expanded or alternative functional forms explored for galaxies where the rational taper is a poor fit. A robust estimate of the $k$ distribution (e.g., using censored-data statistics to account for boundary hits) would strengthen any universality claim.
2.  **Independent M33 validation:** Obtain Corbelli's published velocity decomposition (e.g., from their Fig. 7 or by direct communication) to provide a true independent benchmark for our Casertano solver, beyond the current self-consistency check.
3.  **Investigate the $\omega$–Luminosity relation:** The apparent inverse correlation warrants formal statistical testing (Spearman rank, partial correlations controlling for distance).
4.  **Physical Interpretation:** The success of the **Linear + Tapered** form over the quadrature model is a robust result. Whether this supports a geometric/kinematic origin (e.g., frame dragging, cosmic vorticity) or is merely a phenomenological preference requires theoretical work beyond the scope of this empirical analysis.
