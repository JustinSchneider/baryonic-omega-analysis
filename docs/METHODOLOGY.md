# Methodology: Baryonic Omega Analysis

**Authors:** Justin Schneider, David C. Flynn, Jim Cannaliato
**Last Updated:** February 2026
**Phase:** I — Infrastructure & Calibration (complete through batch analysis)

---

## 1. Scientific Background

The empirical "Omega" ($\omega$) velocity correction model (Flynn & Cannaliato 2025) proposes a linear correction term to galaxy rotation curves. This project upgrades the original point-mass approximation by fitting $\omega$ to the **full baryonic mass decomposition** (Gas + Disk + Bulge), following the methods established by Corbelli & Salucci (2000) and the SPARC catalog (Lelli et al. 2016).

**Hypothesis:** Fitting $\omega$ to the baryonic potential should yield a tighter, more universal correlation than fitting to the Keplerian decline alone.

---

## 2. Data Sources

### 2.1 M33 Calibration Target — Corbelli et al. (2014)

- **Paper:** "Dynamical signatures of a $\Lambda$CDM-halo and the distribution of the baryons in M33" (A&A 572, A23)
- **Data:** Table 1 — 58 radial bins from $R = 0.24$ to $22.72$ kpc
- **Columns:** $V_r(R)$, $\sigma_V(R)$, $\Sigma_{HI}(R)$, $\Sigma_*(R)$
- **Physical Parameters:** $D = 0.84$ Mpc, $i = 52°$
- **Local file:** `data/extracted/corbelli2014_table1.csv`

M33 is not in the SPARC catalog. Since Corbelli provides surface densities rather than velocity components, we convert $\Sigma \to V_{circ}$ using the Casertano thin-disk method (Section 3.2).

### 2.2 SPARC Database — Lelli, McGaugh, & Schombert (2016)

- **Paper:** "SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves" (AJ, 152, 157)
- **Source:** http://astroweb.cwru.edu/SPARC/
- **Format:** Per-galaxy `.dat` files with columns: `Rad`, `Vobs`, `errV`, `Vgas`, `Vdisk`, `Vbul`, `SBdisk`, `SBbul`
- **Note:** Velocity components are provided at $\Upsilon = 1$. Mass-to-light ratios are applied during analysis.

---

## 3. Baryonic Velocity Computation

### 3.1 SPARC Galaxies (Pre-computed Components)

For SPARC galaxies, velocity components are provided directly. We compute the total baryonic velocity following Lelli et al. (2016) Eq. 2:

$$V_{bary}(R) = \sqrt{|V_{gas}| \cdot V_{gas} + \Upsilon_d \cdot |V_{disk}| \cdot V_{disk} + \Upsilon_b \cdot |V_{bulge}| \cdot V_{bulge}}$$

The $|V| \cdot V$ pattern preserves the sign so that negative $V^2$ contributions (e.g., from central gas depressions) reduce $V_{bary}$ rather than producing imaginary numbers.

**Implementation:** `src/physics.py :: compute_v_bary()`

### 3.2 Corbelli Data (Surface Density → Circular Velocity)

For M33, surface densities ($\Sigma$ in $M_\odot/\text{pc}^2$) are converted to circular velocity contributions using the **Casertano (1983) thin-disk method**:

1. The disk surface density profile is interpolated onto a fine, evenly-spaced radial grid (5000+ points)
2. The gravitational potential is computed by summing ring contributions using complete elliptic integrals (Binney & Tremaine 2008, Eq. 2.188)
3. Circular velocity: $V^2_{circ}(R) = R \cdot d\Phi/dR$, obtained by numerical differentiation
4. A softening length of $0.1 \times$ grid spacing avoids the $K(k^2)$ divergence at $R = R'$

For HI gas, a helium correction factor of 1.33 is applied: $\Sigma_{gas} = 1.33 \times \Sigma_{HI}$.

**Implementation:** `src/physics.py :: circular_velocity_thin_disk()`

### 3.3 Mass-to-Light Ratio Conventions

| Parameter | Symbol | Default | Allowed Range | Notes |
|-----------|--------|---------|---------------|-------|
| Disk M/L | $\Upsilon_{disk}$ | 0.5 | 0.3 – 0.8 | SPARC default (3.6 $\mu$m band) |
| Bulge M/L | $\Upsilon_{bulge}$ | 0.7 | 0.3 – 0.8 | SPARC default |

**Phase I (current):** $\Upsilon$ values are held constant at SPARC defaults.
**Future work:** $\Upsilon$ may be treated as free parameters in the fit (bounded 0.3–0.8).

---

## 4. Omega Fitting Method

### 4.1 Models

Two competing functional forms are tested (Notebook 02):

**Model A — Linear (Flynn & Cannaliato 2025):**

$$V_{model}(R) = V_{bary}(R) + \omega \cdot R$$

This treats $\omega$ as a kinematic velocity boost — velocities add directly, analogous to frame-dragging or Coriolis-type effects.

**Model B — Quadrature (Standard force addition):**

$$V_{model}(R) = \sqrt{V_{bary}^2(R) + (\omega \cdot R)^2}$$

This treats $\omega$ as a dynamical force contribution whose potential adds in quadrature with the baryonic potential, analogous to how dark matter halo models combine with baryonic components (Lelli et al. 2016). The force balance is:

$$\frac{V^2}{R} = F_{grav} + F_\omega \implies V^2 = V_{bary}^2 + (\omega R)^2$$

Both models have a single free parameter ($\omega$, in km/s/kpc). Model selection uses BIC (Section 4.4).

**Implementation:** `src/physics.py :: fit_omega()` (linear), `fit_omega_quadrature()` (quadrature)

### 4.2 Fitting Procedure

**Objective:** Minimize weighted chi-squared:

$$\chi^2 = \sum_i \left[ \frac{V_{obs,i} - V_{model,i}}{\sigma_{V,i}} \right]^2$$

**Implementation details:**

- **Optimizer:** `scipy.optimize.curve_fit` with `absolute_sigma=True`
- **Free parameter:** $\omega$ only (single-parameter fit)
- **Initial guess:** $\omega_0 = 0$
- **Bounds:** $[-50, 50]$ km/s/kpc
- **Error handling:** Zero or negative errors are replaced with the minimum nonzero error (or 1.0 km/s if all are zero)

### 4.3 Diagnostics

For each fit, we compute and store:

| Metric | Description |
|--------|-------------|
| $\omega \pm \sigma_\omega$ | Best-fit value and 1-$\sigma$ uncertainty from covariance matrix |
| $\chi^2_\nu$ | Reduced chi-squared ($\text{dof} = n - 1$) |
| RMSE | Root mean square of residuals (km/s) |
| Converged | Whether `curve_fit` converged successfully |
| Flag: $V_{obs} < V_{bary}$ | True if observed velocity falls below baryonic at any radius |

**Implementation:** `src/physics.py :: fit_omega()`, `fit_omega_quadrature()` → `OmegaFitResult` dataclass

### 4.4 Model Selection: BIC

Both models have $k = 1$ free parameter. We use the Bayesian Information Criterion (Schwarz 1978) for model comparison:

$$\text{BIC} = \chi^2 + k \ln(n)$$

where $n$ is the number of data points. Lower BIC is preferred. The difference $\Delta\text{BIC} = \text{BIC}_A - \text{BIC}_B$ follows the Kass & Raftery (1995) interpretation:

| $|\Delta\text{BIC}|$ | Evidence |
|---|---|
| < 2 | Not worth mentioning |
| 2–6 | Positive |
| 6–10 | Strong |
| > 10 | Very strong |

We also compute RMSE restricted to the outer disk ($R > 10$ kpc) where the two models diverge most, as a secondary discriminant.

**Implementation:** `src/physics.py :: compute_bic()`
**Results:** `results/tables/M33_model_comparison.csv`

### 4.5 Tapered Models (Notebook 03)

The pure linear model ($\omega R$) grows without bound, conflicting with observed flat rotation curves. We test two saturation forms:

**Model C — Rational Taper:**

$$V_{model}(R) = V_{bary}(R) + \frac{\omega \cdot R}{1 + R/R_t}$$

At small $R$ this reduces to $\omega R$ (linear); at large $R$ it saturates to $\omega R_t$ (constant). Two free parameters: $\omega$ (km/s/kpc) and $R_t$ (kpc).

**Model D — Tanh Taper:**

$$V_{model}(R) = V_{bary}(R) + V_{max} \cdot \tanh(R / R_t)$$

Two free parameters: $V_{max}$ (km/s) and $R_t$ (kpc). The effective inner-disk $\omega = V_{max} / R_t$.

**Fitting details:**
- **Optimizer:** `scipy.optimize.curve_fit` with `absolute_sigma=True`
- **Initial guesses:** $\omega_0 = 5$, $R_{t,0} = 5$ (rational); $V_{max,0} = 80$, $R_{t,0} = 5$ (tanh)
- **Bounds:** $\omega \in [0, 50]$, $R_t \in [0.1, 50]$ (rational); $V_{max} \in [0, 500]$, $R_t \in [0.1, 50]$ (tanh)
- **Degrees of freedom:** $\text{dof} = n - 2$ (two free parameters)

**Implementation:** `src/physics.py :: fit_omega_tapered()` (rational), `fit_omega_tanh()` (tanh)
**Results:** `results/tables/M33_tapered_results.csv`

### 4.6 Universal Coupling Parameterization (Notebook 04)

To test whether the transition radius $R_t$ is related to galaxy structure, we reparameterize:

$$R_t = k \cdot R_d$$

where $R_d$ is the disk scale length from the SPARC photometric catalog (fixed per galaxy) and $k$ is a dimensionless coupling constant (fit parameter). The model becomes:

$$V_{model}(R) = V_{bary}(R) + \frac{\omega \cdot R}{1 + R / (k \cdot R_d)}$$

If $k$ is approximately constant across galaxies, this implies the saturation scale is set by the baryonic disk.

**Fitting details:**
- **Free parameters:** $\omega$ and $k$ (two parameters per galaxy)
- **Bounds:** $\omega \in [0, 200]$, $k \in [0.1, 20]$
- **Initial guesses:** $\omega_0 = 5$, $k_0 = 2$
- **Quality filter:** Fits with $\chi^2_{red} < 5$ are considered well-constrained

**Implementation:** `src/physics.py :: fit_omega_tapered_kRd()`
**Results:** `results/tables/SPARC_tapered_batch_results.csv`

---

## 5. Validation Criteria

### 5.1 M33 Pipeline Self-Consistency

- **Test:** Round-trip check — compare $V_{bary}$ recomputed during analysis against $V_{bary}$ computed during database ingestion, both from the same Corbelli (2014) surface density inputs
- **Success criterion:** 0% deviation (deterministic pipeline)
- **Note:** This verifies internal consistency (no bugs, rounding errors, or parameter mismatches between pipeline stages). It does **not** constitute an independent comparison against Corbelli's published velocity decomposition, which would require their original $V_{gas}$ and $V_{stars}$ values (not available in their Table 1). An independent cross-check against Corbelli's Fig. 7 mass models remains a future validation step.
- **Implementation:** `src/physics.py :: compute_validation_metrics()`
- **Output:** `results/tables/M33_validation.csv`

### 5.2 Fit Convergence

- **Target:** `scipy.curve_fit` converges for >90% of SPARC galaxies
- **Phase I result:** 118/118 galaxies achieved numerical convergence (100%); 89/118 (75%) met the $\chi^2_{red} < 5$ quality threshold
- **Edge cases:** Galaxies with $V_{obs} < V_{bary}$ anywhere are flagged in the database

### 5.3 Database Integrity

- **Round-trip check:** Ingest SPARC data → SQLite → Pandas DataFrame, verify floats match original text files

---

## 6. Software Architecture

### Pipeline

```
Raw Data (.dat/.csv) → Ingestion (src/ingest.py) → SQLite DB
     ↓
Query profiles → Compute V_bary → Fit omega → Store results
     ↓
Export tables (CSV) + Generate figures (PNG)
```

### Database Schema (SQLite)

- **`galaxies`** — Metadata (distance, inclination, luminosity, quality flag, data source)
- **`radial_profiles`** — Per-galaxy radial curves ($R$, $V_{obs}$, $V_{err}$, $V_{gas}$, $V_{disk}$, $V_{bulge}$)
- **`omega_fits`** — Fit results ($\omega$, $\chi^2$, RMSE, method version, timestamp)

### Versioning

Each fit stores a `method_version` string so results from different fitting algorithms can be compared. Current versions:

| Version | Description |
|---------|-------------|
| `v1_fixed_ML` | Linear model, fixed $\Upsilon$ at SPARC defaults |
| `v1_quadrature` | Quadrature model, fixed $\Upsilon$ |
| `v1_rational_taper` | Rational taper, free $\omega$ and $R_t$ |
| `v1_tanh_taper` | Tanh taper, free $V_{max}$ and $R_t$ |
| `v2_kRd_taper` | Rational taper with $R_t = k \cdot R_d$, free $\omega$ and $k$ |

---

## 7. Sensitivity Analysis

For the M33 calibration, we test the sensitivity of $\omega$ to the disk mass-to-light ratio:

- $\Upsilon_{disk}$ varied from 0.3 to 0.8 in steps of 0.05
- $\Upsilon_{bulge}$ held fixed at 0.7 (M33 has negligible bulge)
- Results stored in `results/tables/M33_sensitivity.csv`

This establishes how strongly the omega measurement depends on the assumed stellar mass normalization.

---

## References

1. Binney, J. & Tremaine, S. (2008). *Galactic Dynamics* (2nd ed.). Princeton University Press.
2. Casertano, S. (1983). MNRAS, 203, 735. "Rotation curve of the edge-on spiral galaxy NGC 5907."
3. Corbelli, E. & Salucci, P. (2000). MNRAS, 311, 441. "The Extended Rotation Curve and the Dark Matter Halo of M33."
4. Corbelli, E., et al. (2014). A&A, 572, A23. "Dynamical signatures of a $\Lambda$CDM-halo and the distribution of the baryons in M33."
5. Flynn, D. C. & Cannaliato, J. (2025). "A New Empirical Fit to Galaxy Rotation Curves."
6. Kass, R. E. & Raftery, A. E. (1995). JASA, 90, 773. "Bayes Factors."
7. Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). AJ, 152, 157. "SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves."
8. Schwarz, G. (1978). Annals of Statistics, 6, 461. "Estimating the Dimension of a Model."
