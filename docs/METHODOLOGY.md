# Methodology: Baryonic Omega Analysis

**Authors:** Justin Schneider, David C. Flynn, Jim Cannaliato
**Last Updated:** February 2026
**Phase:** I — Infrastructure & Calibration

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

**Implementation:** `src/physics.py :: surface_density_to_v_circ()`

### 3.3 Mass-to-Light Ratio Conventions

| Parameter | Symbol | Default | Allowed Range | Notes |
|-----------|--------|---------|---------------|-------|
| Disk M/L | $\Upsilon_{disk}$ | 0.5 | 0.3 – 0.8 | SPARC default (3.6 $\mu$m band) |
| Bulge M/L | $\Upsilon_{bulge}$ | 0.7 | 0.3 – 0.8 | SPARC default |

**Phase I (current):** $\Upsilon$ values are held constant at SPARC defaults.
**Future work:** $\Upsilon$ may be treated as free parameters in the fit (bounded 0.3–0.8).

---

## 4. Omega Fitting Method

### 4.1 Model

$$V_{model}(R) = V_{bary}(R) + \omega \cdot R$$

The parameter $\omega$ (units: km/s/kpc) represents a linear velocity correction that scales with galactocentric radius. Physically, this may correspond to a correction to Newtonian dynamics or a coupling to a background field (Flynn & Cannaliato 2025).

### 4.2 Fitting Procedure

**Objective:** Minimize weighted chi-squared:

$$\chi^2 = \sum_i \left[ \frac{V_{obs,i} - V_{bary,i} - \omega \cdot R_i}{\sigma_{V,i}} \right]^2$$

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

**Implementation:** `src/physics.py :: fit_omega()` → `OmegaFitResult` dataclass

---

## 5. Validation Criteria

### 5.1 M33 Benchmark

- **Test:** Compare our $V_{bary}$ against Corbelli's published $V_{stars} + V_{gas}$
- **Success criterion:** Agreement within 5% at all radii $R > 2$ kpc
- **Implementation:** `src/physics.py :: compute_validation_metrics()`
- **Output:** `results/tables/M33_validation.csv`

### 5.2 Fit Convergence (Phase II)

- **Target:** `scipy.curve_fit` converges for >90% of SPARC galaxies
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

Each fit stores a `method_version` string (e.g., `v1_fixed_ML`) so results from different fitting algorithms can be compared. Current version: `v1_fixed_ML` (fixed mass-to-light ratios).

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
6. Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). AJ, 152, 157. "SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves."
