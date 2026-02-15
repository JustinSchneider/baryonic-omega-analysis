"""Core physics equations for baryonic velocity computation and omega fitting.

Implements the Corbelli method: decompose observed rotation curves into
baryonic components (gas, disk, bulge) and fit the omega correction parameter.

Key references:
  - Lelli, McGaugh, & Schombert (2016), Eq. 2 for the sign convention on
    velocity components.
  - Casertano (1983) for the thin-disk gravitational potential from surface
    density profiles via ring summation with elliptic integrals.
  - Binney & Tremaine (2008), Eq. 2.188 for the axisymmetric disk potential.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import ellipe, ellipk

from src.utils import setup_logger

logger = setup_logger(__name__)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

# Gravitational constant in units: pc * (km/s)^2 / M_sun
G_PC = 4.302e-3

# Multiplicative factor to convert HI mass to total gas mass (HI + He)
HELIUM_FACTOR = 1.33

# 1 kpc = 1000 pc
KPC_TO_PC = 1000.0


# ---------------------------------------------------------------------------
# Thin-disk gravitational potential (Casertano 1983)
# ---------------------------------------------------------------------------


def circular_velocity_thin_disk(
    r_eval: np.ndarray,
    r_profile: np.ndarray,
    sigma_profile: np.ndarray,
    helium_factor: float = 1.0,
) -> np.ndarray:
    """Compute circular velocity from a tabulated surface density profile.

    Uses the Casertano (1983) method: the disk surface density is interpolated
    onto a fine, evenly-spaced grid and the gravitational potential is summed
    from ring contributions using complete elliptic integrals. The circular
    velocity is obtained by numerical differentiation of the potential.

    A small softening length (equivalent to finite disk thickness) avoids the
    logarithmic singularity that occurs when evaluation and ring radii coincide.

    Based on Binney & Tremaine (2008), Eq. 2.188 for the ring potential:
        Phi(R) = -G * dM / (pi * (R + a)) * K(k^2)
    where k^2 = 4*R*a / (R + a)^2.

    Args:
        r_eval: Radii at which to evaluate V_circ (kpc).
        r_profile: Radii of the surface density samples (kpc).
        sigma_profile: Surface density at each r_profile (M_sun/pc^2).
        helium_factor: Multiplicative factor for sigma (use 1.33 for HI gas
            to account for helium; use 1.0 for stellar mass).

    Returns:
        Circular velocity (km/s) at each r_eval. Always non-negative.
    """
    r_eval = np.asarray(r_eval, dtype=np.float64)
    r_profile = np.asarray(r_profile, dtype=np.float64)
    sigma_profile = np.asarray(sigma_profile, dtype=np.float64)

    R_eval_pc = r_eval * KPC_TO_PC
    Rp = r_profile * KPC_TO_PC
    sigma = sigma_profile * helium_factor

    # Interpolate the surface density onto a fine, evenly-spaced grid.
    # This avoids discrete ring artifacts and ensures smooth potential.
    r_min_pc = max(Rp.min() * 0.5, 5.0)
    r_max_pc = max(Rp.max(), R_eval_pc.max()) * 1.3
    n_fine = max(5000, len(Rp) * 20)
    R_fine = np.linspace(r_min_pc, r_max_pc, n_fine)
    sigma_fine = np.interp(R_fine, Rp, sigma, left=0.0, right=0.0)

    # Ring widths on the fine grid (uniform)
    dr_fine = R_fine[1] - R_fine[0]

    # Ring masses
    dM_fine = 2.0 * np.pi * R_fine * sigma_fine * dr_fine

    # Softening: model the disk with a small effective thickness.
    # Use a fraction of the fine grid spacing — small enough for accuracy,
    # large enough to remove the K(k^2) divergence at k=1.
    eps = dr_fine * 0.1
    eps2 = eps ** 2

    # Compute potential on an evaluation grid that covers all requested radii.
    # Use a grid offset from the ring positions by half a step.
    R_grid = R_fine + dr_fine * 0.5

    Phi_grid = np.zeros_like(R_grid)
    for j in range(n_fine):
        if dM_fine[j] <= 0:
            continue
        # Softened denominator: sqrt((R + R')^2 + eps^2)
        sum_R_sq = (R_grid + R_fine[j]) ** 2 + eps2
        k2 = 4.0 * R_grid * R_fine[j] / sum_R_sq
        k2 = np.clip(k2, 0.0, 1.0 - 1e-10)
        Phi_grid += -2.0 * G_PC * dM_fine[j] / (np.pi * np.sqrt(sum_R_sq)) * ellipk(k2)

    # Numerical derivative: dPhi/dR
    dPhi_dR = np.gradient(Phi_grid, R_grid)

    # V^2 = R * dPhi/dR (dPhi/dR > 0 for net inward force)
    v2_grid = R_grid * dPhi_dR
    v2_grid = np.maximum(v2_grid, 0.0)
    v_grid = np.sqrt(v2_grid)

    # Interpolate to the requested evaluation radii
    v_circ = np.interp(R_eval_pc, R_grid, v_grid, left=0.0, right=0.0)

    return v_circ


@dataclass
class OmegaFitResult:
    """Container for all results from an omega fit."""

    galaxy_id: str
    omega_value: float
    omega_uncertainty: float
    chi_squared: float
    reduced_chi_squared: float
    residuals_rmse: float
    n_points: int
    converged: bool
    flag_v_obs_lt_v_bary: bool
    method_version: str
    upsilon_disk: float
    upsilon_bulge: float
    # Arrays for plotting (not stored in DB)
    v_bary: np.ndarray = field(repr=False)
    v_model: np.ndarray = field(repr=False)
    residuals: np.ndarray = field(repr=False)

    def to_dict(self) -> dict:
        """Convert to dict suitable for database insertion (excludes arrays)."""
        return {
            "galaxy_id": self.galaxy_id,
            "omega_value": self.omega_value,
            "omega_uncertainty": self.omega_uncertainty,
            "chi_squared": self.chi_squared,
            "reduced_chi_squared": self.reduced_chi_squared,
            "residuals_rmse": self.residuals_rmse,
            "n_points": self.n_points,
            "converged": self.converged,
            "flag_v_obs_lt_v_bary": self.flag_v_obs_lt_v_bary,
            "method_version": self.method_version,
            "upsilon_disk": self.upsilon_disk,
            "upsilon_bulge": self.upsilon_bulge,
        }

    def to_summary_dataframe(self) -> "pd.DataFrame":
        """Create a single-row summary DataFrame for display and CSV export."""
        import pandas as pd

        row = {
            "galaxy": self.galaxy_id,
            "omega_km_s_kpc": self.omega_value,
            "omega_err": self.omega_uncertainty,
            "chi2_reduced": self.reduced_chi_squared,
            "rmse_km_s": self.residuals_rmse,
            "n_points": self.n_points,
            "converged": self.converged,
            "upsilon_disk": self.upsilon_disk,
            "upsilon_bulge": self.upsilon_bulge,
            "method_version": self.method_version,
            "flag_v_obs_lt_v_bary": self.flag_v_obs_lt_v_bary,
        }
        return pd.DataFrame([row])

    def to_radial_profile_dataframe(
        self,
        radius: np.ndarray,
        v_obs: np.ndarray,
        v_err: np.ndarray,
        v_gas: np.ndarray,
        v_disk: np.ndarray,
        v_bulge: np.ndarray,
    ) -> "pd.DataFrame":
        """Build an augmented radial profile DataFrame with model columns."""
        import pandas as pd

        return pd.DataFrame({
            "radius_kpc": radius,
            "v_obs": v_obs,
            "v_err": v_err,
            "v_gas": v_gas,
            "v_disk": v_disk,
            "v_bulge": v_bulge,
            "v_bary": self.v_bary,
            "v_model": self.v_model,
            "residual": self.residuals,
        })


def compute_v_bary(
    v_gas: np.ndarray,
    v_disk: np.ndarray,
    v_bulge: np.ndarray,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
) -> np.ndarray:
    """Compute total baryonic velocity using the Lelli et al. (2016) Eq. 2 convention.

    V_bary(R) = sqrt(|V_gas|*V_gas + Upsilon_d * |V_disk|*V_disk
                     + Upsilon_b * |V_bulge|*V_bulge)

    The |V|*V pattern preserves the sign: if V is negative (e.g., gas central
    depression), its squared contribution is negative, reducing V_bary.

    SPARC provides V_disk and V_bulge at Upsilon=1, so Upsilon enters as a
    direct multiplier on V^2 (equivalent to sqrt(Upsilon) scaling on V).

    Args:
        v_gas: Gas velocity component (km/s). May contain negative values.
        v_disk: Disk velocity component at Upsilon=1 (km/s).
        v_bulge: Bulge velocity component at Upsilon=1 (km/s).
        upsilon_disk: Stellar mass-to-light ratio for disk. Default 0.5.
        upsilon_bulge: Stellar mass-to-light ratio for bulge. Default 0.7.

    Returns:
        Total baryonic velocity (km/s). Always non-negative.
    """
    v_gas = np.asarray(v_gas, dtype=np.float64)
    v_disk = np.asarray(v_disk, dtype=np.float64)
    v_bulge = np.asarray(v_bulge, dtype=np.float64)

    v2_gas = np.abs(v_gas) * v_gas
    v2_disk = upsilon_disk * np.abs(v_disk) * v_disk
    v2_bulge = upsilon_bulge * np.abs(v_bulge) * v_bulge

    v2_total = v2_gas + v2_disk + v2_bulge

    # If total is negative at some radius, the baryonic model is unphysical there
    negative_mask = v2_total < 0
    if np.any(negative_mask):
        n_neg = np.sum(negative_mask)
        logger.warning(
            "Negative V^2_bary at %d point(s) — setting V_bary=0 there", n_neg
        )
        v2_total = np.where(negative_mask, 0.0, v2_total)

    return np.sqrt(v2_total)


def fit_omega(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    galaxy_id: str = "unknown",
    method_version: str = "v1_fixed_ML",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    omega_bounds: tuple = (-50.0, 50.0),
) -> OmegaFitResult:
    """Fit the omega parameter to an observed rotation curve.

    Minimizes chi^2 = sum_i [(V_obs_i - V_bary_i - omega*R_i) / V_err_i]^2
    using scipy.optimize.curve_fit with error weighting.

    Args:
        radius: Radial positions (kpc).
        v_obs: Observed velocities (km/s).
        v_err: Velocity errors (km/s). Values <= 0 are replaced.
        v_bary: Pre-computed baryonic velocity (km/s).
        galaxy_id: Galaxy identifier for labeling.
        method_version: Version string for reproducibility.
        upsilon_disk: Mass-to-light ratio used for V_bary.
        upsilon_bulge: Mass-to-light ratio used for V_bary.
        omega_bounds: (lower, upper) bounds for omega in km/s/kpc.

    Returns:
        OmegaFitResult with all fit diagnostics and arrays.
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)

    # Replace zero/negative errors with minimum nonzero error (or 1.0 km/s)
    positive_errs = v_err[v_err > 0]
    min_err = float(np.min(positive_errs)) if len(positive_errs) > 0 else 1.0
    v_err_safe = np.where(v_err > 0, v_err, min_err)
    if np.any(v_err <= 0):
        logger.warning(
            "%s: replaced %d zero/negative errors with %.2f km/s",
            galaxy_id,
            np.sum(v_err <= 0),
            min_err,
        )

    # Check for V_obs < V_bary at any radius
    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    # Define the model as a closure over v_bary
    def _model(r, omega):
        return v_bary + omega * r

    try:
        popt, pcov = curve_fit(
            _model,
            radius,
            v_obs,
            p0=[0.0],
            sigma=v_err_safe,
            absolute_sigma=True,
            bounds=([omega_bounds[0]], [omega_bounds[1]]),
        )
        omega_best = float(popt[0])
        omega_err = float(np.sqrt(np.diag(pcov))[0])
        converged = True
    except RuntimeError as e:
        logger.warning("%s: curve_fit did not converge — %s", galaxy_id, e)
        omega_best = float("nan")
        omega_err = float("nan")
        converged = False

    # Compute diagnostics
    v_model = _model(radius, omega_best)
    residuals = v_obs - v_model
    chi2 = float(np.sum((residuals / v_err_safe) ** 2))
    dof = max(n_points - 1, 1)  # 1 free parameter
    reduced_chi2 = chi2 / dof
    rmse = float(np.sqrt(np.mean(residuals**2)))

    result = OmegaFitResult(
        galaxy_id=galaxy_id,
        omega_value=omega_best,
        omega_uncertainty=omega_err,
        chi_squared=chi2,
        reduced_chi_squared=reduced_chi2,
        residuals_rmse=rmse,
        n_points=n_points,
        converged=converged,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary,
        v_model=v_model,
        residuals=residuals,
    )

    if converged:
        logger.info(
            "%s: omega=%.4f +/- %.4f km/s/kpc  chi2_r=%.2f  RMSE=%.2f km/s",
            galaxy_id,
            omega_best,
            omega_err,
            reduced_chi2,
            rmse,
        )

    return result


def fit_omega_quadrature(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    galaxy_id: str = "unknown",
    method_version: str = "v1_quadrature",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    omega_bounds: tuple = (-50.0, 50.0),
) -> OmegaFitResult:
    """Fit the omega parameter using the quadrature (force-addition) model.

    Model: V_model = sqrt(V_bary^2 + (omega * R)^2)

    This treats omega*R as a dynamical force contribution whose potential
    adds in quadrature with the baryonic potential, analogous to how
    standard dark matter halo models combine with baryonic components
    (e.g., Lelli et al. 2016).

    Physical interpretation: omega acts as a background force field
    (centrifugal-like) rather than a kinematic velocity boost.

    Args:
        radius: Radial positions (kpc).
        v_obs: Observed velocities (km/s).
        v_err: Velocity errors (km/s). Values <= 0 are replaced.
        v_bary: Pre-computed baryonic velocity (km/s).
        galaxy_id: Galaxy identifier for labeling.
        method_version: Version string for reproducibility.
        upsilon_disk: Mass-to-light ratio used for V_bary.
        upsilon_bulge: Mass-to-light ratio used for V_bary.
        omega_bounds: (lower, upper) bounds for omega in km/s/kpc.

    Returns:
        OmegaFitResult with all fit diagnostics and arrays.
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)

    # Replace zero/negative errors with minimum nonzero error (or 1.0 km/s)
    positive_errs = v_err[v_err > 0]
    min_err = float(np.min(positive_errs)) if len(positive_errs) > 0 else 1.0
    v_err_safe = np.where(v_err > 0, v_err, min_err)
    if np.any(v_err <= 0):
        logger.warning(
            "%s: replaced %d zero/negative errors with %.2f km/s",
            galaxy_id,
            np.sum(v_err <= 0),
            min_err,
        )

    # Check for V_obs < V_bary at any radius
    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    # Quadrature model: V = sqrt(V_bary^2 + (omega * R)^2)
    def _model(r, omega):
        return np.sqrt(v_bary**2 + (omega * r) ** 2)

    try:
        popt, pcov = curve_fit(
            _model,
            radius,
            v_obs,
            p0=[0.0],
            sigma=v_err_safe,
            absolute_sigma=True,
            bounds=([omega_bounds[0]], [omega_bounds[1]]),
        )
        omega_best = float(popt[0])
        omega_err = float(np.sqrt(np.diag(pcov))[0])
        converged = True
    except RuntimeError as e:
        logger.warning("%s: curve_fit did not converge — %s", galaxy_id, e)
        omega_best = float("nan")
        omega_err = float("nan")
        converged = False

    # Compute diagnostics
    v_model = _model(radius, omega_best)
    residuals = v_obs - v_model
    chi2 = float(np.sum((residuals / v_err_safe) ** 2))
    dof = max(n_points - 1, 1)  # 1 free parameter
    reduced_chi2 = chi2 / dof
    rmse = float(np.sqrt(np.mean(residuals**2)))

    result = OmegaFitResult(
        galaxy_id=galaxy_id,
        omega_value=omega_best,
        omega_uncertainty=omega_err,
        chi_squared=chi2,
        reduced_chi_squared=reduced_chi2,
        residuals_rmse=rmse,
        n_points=n_points,
        converged=converged,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary,
        v_model=v_model,
        residuals=residuals,
    )

    if converged:
        logger.info(
            "%s [quadrature]: omega=%.4f +/- %.4f km/s/kpc  chi2_r=%.2f  RMSE=%.2f km/s",
            galaxy_id,
            omega_best,
            omega_err,
            reduced_chi2,
            rmse,
        )

    return result


def fit_omega_tapered(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    galaxy_id: str = "unknown",
    method_version: str = "v1_rational_taper",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    omega_bounds: tuple = (0.0, 50.0),
    rt_bounds: tuple = (0.1, 50.0),
) -> "TaperedFitResult":
    """Fit the rational-taper (saturation) model to a rotation curve.

    Model: V_model = V_bary + omega * R / (1 + R / R_t)

    At small R the correction is linear (omega * R); at large R it
    saturates to a constant velocity (omega * R_t).

    Args:
        radius: Radial positions (kpc).
        v_obs: Observed velocities (km/s).
        v_err: Velocity errors (km/s). Values <= 0 are replaced.
        v_bary: Pre-computed baryonic velocity (km/s).
        galaxy_id: Galaxy identifier for labeling.
        method_version: Version string for reproducibility.
        upsilon_disk: Mass-to-light ratio used for V_bary.
        upsilon_bulge: Mass-to-light ratio used for V_bary.
        omega_bounds: (lower, upper) bounds for omega in km/s/kpc.
        rt_bounds: (lower, upper) bounds for R_t in kpc.

    Returns:
        TaperedFitResult with all fit diagnostics and arrays.
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)

    positive_errs = v_err[v_err > 0]
    min_err = float(np.min(positive_errs)) if len(positive_errs) > 0 else 1.0
    v_err_safe = np.where(v_err > 0, v_err, min_err)

    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    def _model(r, omega, r_t):
        return v_bary + omega * r / (1.0 + r / r_t)

    try:
        popt, pcov = curve_fit(
            _model,
            radius,
            v_obs,
            p0=[5.0, 5.0],
            sigma=v_err_safe,
            absolute_sigma=True,
            bounds=(
                [omega_bounds[0], rt_bounds[0]],
                [omega_bounds[1], rt_bounds[1]],
            ),
        )
        omega_best = float(popt[0])
        rt_best = float(popt[1])
        perr = np.sqrt(np.diag(pcov))
        omega_err = float(perr[0])
        rt_err = float(perr[1])
        converged = True
    except RuntimeError as e:
        logger.warning("%s: rational taper fit did not converge — %s", galaxy_id, e)
        omega_best = rt_best = omega_err = rt_err = float("nan")
        converged = False

    v_model = _model(radius, omega_best, rt_best)
    residuals = v_obs - v_model
    chi2 = float(np.sum((residuals / v_err_safe) ** 2))
    dof = max(n_points - 2, 1)  # 2 free parameters
    reduced_chi2 = chi2 / dof
    rmse = float(np.sqrt(np.mean(residuals**2)))

    result = TaperedFitResult(
        galaxy_id=galaxy_id,
        model_name="rational_taper",
        omega_value=omega_best,
        omega_uncertainty=omega_err,
        rt_value=rt_best,
        rt_uncertainty=rt_err,
        chi_squared=chi2,
        reduced_chi_squared=reduced_chi2,
        residuals_rmse=rmse,
        n_points=n_points,
        n_params=2,
        converged=converged,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary,
        v_model=v_model,
        residuals=residuals,
    )

    if converged:
        logger.info(
            "%s [rational taper]: omega=%.4f +/- %.4f  R_t=%.4f +/- %.4f  chi2_r=%.2f  RMSE=%.2f",
            galaxy_id, omega_best, omega_err, rt_best, rt_err, reduced_chi2, rmse,
        )

    return result


def fit_omega_tanh(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    galaxy_id: str = "unknown",
    method_version: str = "v1_tanh_taper",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    vmax_bounds: tuple = (0.0, 500.0),
    rt_bounds: tuple = (0.1, 50.0),
) -> "TaperedFitResult":
    """Fit the tanh-taper (shear) model to a rotation curve.

    Model: V_model = V_bary + V_max * tanh(R / R_t)

    At small R the correction is approximately linear (V_max * R / R_t);
    at large R it saturates to V_max.

    Args:
        radius: Radial positions (kpc).
        v_obs: Observed velocities (km/s).
        v_err: Velocity errors (km/s). Values <= 0 are replaced.
        v_bary: Pre-computed baryonic velocity (km/s).
        galaxy_id: Galaxy identifier for labeling.
        method_version: Version string for reproducibility.
        upsilon_disk: Mass-to-light ratio used for V_bary.
        upsilon_bulge: Mass-to-light ratio used for V_bary.
        vmax_bounds: (lower, upper) bounds for V_max in km/s.
        rt_bounds: (lower, upper) bounds for R_t in kpc.

    Returns:
        TaperedFitResult with all fit diagnostics and arrays.
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)

    positive_errs = v_err[v_err > 0]
    min_err = float(np.min(positive_errs)) if len(positive_errs) > 0 else 1.0
    v_err_safe = np.where(v_err > 0, v_err, min_err)

    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    def _model(r, v_max, r_t):
        return v_bary + v_max * np.tanh(r / r_t)

    try:
        popt, pcov = curve_fit(
            _model,
            radius,
            v_obs,
            p0=[80.0, 5.0],
            sigma=v_err_safe,
            absolute_sigma=True,
            bounds=(
                [vmax_bounds[0], rt_bounds[0]],
                [vmax_bounds[1], rt_bounds[1]],
            ),
        )
        vmax_best = float(popt[0])
        rt_best = float(popt[1])
        perr = np.sqrt(np.diag(pcov))
        vmax_err = float(perr[0])
        rt_err = float(perr[1])
        converged = True
    except RuntimeError as e:
        logger.warning("%s: tanh taper fit did not converge — %s", galaxy_id, e)
        vmax_best = rt_best = vmax_err = rt_err = float("nan")
        converged = False

    v_model = _model(radius, vmax_best, rt_best)
    residuals = v_obs - v_model
    chi2 = float(np.sum((residuals / v_err_safe) ** 2))
    dof = max(n_points - 2, 1)  # 2 free parameters
    reduced_chi2 = chi2 / dof
    rmse = float(np.sqrt(np.mean(residuals**2)))

    # For tanh model, the effective omega at R=0 is V_max / R_t
    omega_effective = vmax_best / rt_best if rt_best > 0 else float("nan")

    result = TaperedFitResult(
        galaxy_id=galaxy_id,
        model_name="tanh_taper",
        omega_value=omega_effective,
        omega_uncertainty=float("nan"),  # derived quantity, no simple error
        rt_value=rt_best,
        rt_uncertainty=rt_err,
        chi_squared=chi2,
        reduced_chi_squared=reduced_chi2,
        residuals_rmse=rmse,
        n_points=n_points,
        n_params=2,
        converged=converged,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary,
        v_model=v_model,
        residuals=residuals,
        vmax_value=vmax_best,
        vmax_uncertainty=vmax_err,
    )

    if converged:
        logger.info(
            "%s [tanh taper]: V_max=%.4f +/- %.4f  R_t=%.4f +/- %.4f  chi2_r=%.2f  RMSE=%.2f",
            galaxy_id, vmax_best, vmax_err, rt_best, rt_err, reduced_chi2, rmse,
        )

    return result


def fit_omega_tapered_kRd(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    r_d: float,
    galaxy_id: str = "unknown",
    method_version: str = "v2_kRd_taper",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    omega_bounds: tuple = (0.0, 200.0),
    k_bounds: tuple = (0.1, 20.0),
) -> "TaperedFitResult":
    """Fit the rational-taper model with R_t = k * R_d (universal coupling).

    Model: V_model = V_bary + omega * R / (1 + R / (k * R_d))

    Instead of fitting R_t directly, this parameterizes the transition
    radius as a multiple of the disk scale length R_d. If k is universal
    across galaxies, this implies the "dark" transition scale is set by
    the baryonic disk.

    Args:
        radius: Radial positions (kpc).
        v_obs: Observed velocities (km/s).
        v_err: Velocity errors (km/s). Values <= 0 are replaced.
        v_bary: Pre-computed baryonic velocity (km/s).
        r_d: Disk scale length (kpc), fixed from photometry.
        galaxy_id: Galaxy identifier for labeling.
        method_version: Version string for reproducibility.
        upsilon_disk: Mass-to-light ratio used for V_bary.
        upsilon_bulge: Mass-to-light ratio used for V_bary.
        omega_bounds: (lower, upper) bounds for omega in km/s/kpc.
        k_bounds: (lower, upper) bounds for k (dimensionless).

    Returns:
        TaperedFitResult with rt_value = k * R_d and k stored in metadata.
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)

    positive_errs = v_err[v_err > 0]
    min_err = float(np.min(positive_errs)) if len(positive_errs) > 0 else 1.0
    v_err_safe = np.where(v_err > 0, v_err, min_err)

    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    def _model(r, omega, k):
        r_t = k * r_d
        return v_bary + omega * r / (1.0 + r / r_t)

    try:
        popt, pcov = curve_fit(
            _model,
            radius,
            v_obs,
            p0=[5.0, 2.0],
            sigma=v_err_safe,
            absolute_sigma=True,
            bounds=(
                [omega_bounds[0], k_bounds[0]],
                [omega_bounds[1], k_bounds[1]],
            ),
        )
        omega_best = float(popt[0])
        k_best = float(popt[1])
        perr = np.sqrt(np.diag(pcov))
        omega_err = float(perr[0])
        k_err = float(perr[1])
        converged = True
    except RuntimeError as e:
        logger.warning("%s: k*R_d taper fit did not converge — %s", galaxy_id, e)
        omega_best = k_best = omega_err = k_err = float("nan")
        converged = False

    rt_best = k_best * r_d
    rt_err = k_err * r_d  # linear error propagation (R_d is fixed)

    v_model = _model(radius, omega_best, k_best)
    residuals = v_obs - v_model
    chi2 = float(np.sum((residuals / v_err_safe) ** 2))
    dof = max(n_points - 2, 1)  # 2 free parameters
    reduced_chi2 = chi2 / dof
    rmse = float(np.sqrt(np.mean(residuals**2)))

    result = TaperedFitResult(
        galaxy_id=galaxy_id,
        model_name="kRd_taper",
        omega_value=omega_best,
        omega_uncertainty=omega_err,
        rt_value=rt_best,
        rt_uncertainty=rt_err,
        chi_squared=chi2,
        reduced_chi_squared=reduced_chi2,
        residuals_rmse=rmse,
        n_points=n_points,
        n_params=2,
        converged=converged,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary,
        v_model=v_model,
        residuals=residuals,
        k_value=k_best,
        k_uncertainty=k_err,
        r_d_value=r_d,
    )

    if converged:
        logger.info(
            "%s [k*R_d taper]: omega=%.4f +/- %.4f  k=%.4f +/- %.4f  R_t=%.4f  chi2_r=%.2f  RMSE=%.2f",
            galaxy_id, omega_best, omega_err, k_best, k_err, rt_best, reduced_chi2, rmse,
        )

    return result


@dataclass
class TaperedFitResult:
    """Container for results from a tapered omega fit (2-parameter models)."""

    galaxy_id: str
    model_name: str
    omega_value: float
    omega_uncertainty: float
    rt_value: float
    rt_uncertainty: float
    chi_squared: float
    reduced_chi_squared: float
    residuals_rmse: float
    n_points: int
    n_params: int
    converged: bool
    flag_v_obs_lt_v_bary: bool
    method_version: str
    upsilon_disk: float
    upsilon_bulge: float
    v_bary: np.ndarray = field(repr=False)
    v_model: np.ndarray = field(repr=False)
    residuals: np.ndarray = field(repr=False)
    # tanh model stores V_max separately
    vmax_value: Optional[float] = None
    vmax_uncertainty: Optional[float] = None
    # k*R_d model stores coupling constant and disk scale length
    k_value: Optional[float] = None
    k_uncertainty: Optional[float] = None
    r_d_value: Optional[float] = None

    def to_dict(self) -> dict:
        """Convert to dict suitable for database insertion (excludes arrays)."""
        d = {
            "galaxy_id": self.galaxy_id,
            "model_name": self.model_name,
            "omega_value": self.omega_value,
            "omega_uncertainty": self.omega_uncertainty,
            "rt_value": self.rt_value,
            "rt_uncertainty": self.rt_uncertainty,
            "chi_squared": self.chi_squared,
            "reduced_chi_squared": self.reduced_chi_squared,
            "residuals_rmse": self.residuals_rmse,
            "n_points": self.n_points,
            "n_params": self.n_params,
            "converged": self.converged,
            "flag_v_obs_lt_v_bary": self.flag_v_obs_lt_v_bary,
            "method_version": self.method_version,
            "upsilon_disk": self.upsilon_disk,
            "upsilon_bulge": self.upsilon_bulge,
        }
        if self.vmax_value is not None:
            d["vmax_value"] = self.vmax_value
            d["vmax_uncertainty"] = self.vmax_uncertainty
        if self.k_value is not None:
            d["k_value"] = self.k_value
            d["k_uncertainty"] = self.k_uncertainty
            d["r_d_value"] = self.r_d_value
        return d


def compute_bic(n_points: int, k_params: int, chi_squared: float) -> float:
    """Compute the Bayesian Information Criterion (BIC).

    BIC = chi^2 + k * ln(n)

    where chi^2 is the total (non-reduced) chi-squared, k is the number
    of free parameters, and n is the number of data points.

    Lower BIC indicates a preferred model. The difference Delta_BIC
    between two models follows the Kass & Raftery (1995) scale:
      - |Delta_BIC| < 2:   Not worth mentioning
      - 2 < |Delta_BIC| < 6:  Positive evidence
      - 6 < |Delta_BIC| < 10: Strong evidence
      - |Delta_BIC| > 10:     Very strong evidence

    Args:
        n_points: Number of data points.
        k_params: Number of free parameters in the model.
        chi_squared: Total (non-reduced) chi-squared value.

    Returns:
        BIC value.
    """
    return chi_squared + k_params * np.log(n_points)


def compute_validation_metrics(
    radius: np.ndarray,
    v_bary: np.ndarray,
    v_bary_reference: np.ndarray,
    min_radius: float = 2.0,
    threshold_pct: float = 5.0,
) -> "pd.DataFrame":
    """Compare V_bary against a reference and flag deviations.

    Computes percentage difference at each radius above *min_radius* and
    flags points exceeding *threshold_pct*.

    Args:
        radius: Radii in kpc.
        v_bary: Our computed V_bary (km/s).
        v_bary_reference: Reference V_bary to compare against (km/s).
        min_radius: Only evaluate radii above this value (kpc).
        threshold_pct: Flag percentage differences above this value.

    Returns:
        DataFrame with columns: radius_kpc, v_bary, v_bary_ref,
        pct_diff, exceeds_threshold.
    """
    import pandas as pd

    mask = radius > min_radius
    r = radius[mask]
    vb = v_bary[mask]
    vb_ref = v_bary_reference[mask]

    pct_diff = np.where(
        vb_ref != 0,
        100.0 * np.abs(vb - vb_ref) / np.abs(vb_ref),
        0.0,
    )

    return pd.DataFrame({
        "radius_kpc": r,
        "v_bary": vb,
        "v_bary_ref": vb_ref,
        "pct_diff": pct_diff,
        "exceeds_threshold": pct_diff > threshold_pct,
    })
