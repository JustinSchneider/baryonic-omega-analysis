"""Tests for the physics module: omega model and fitting."""

import pytest
import numpy as np
from scipy.special import i0, i1, k0, k1

from src.physics import G_PC, KPC_TO_PC, circular_velocity_thin_disk, compute_v_bary, fit_omega


class TestCircularVelocityThinDisk:
    """Tests for the Casertano (1983) thin-disk velocity calculation."""

    def test_point_mass_limit(self):
        """A narrow ring of mass M should produce ~Keplerian V at R >> a."""
        # Create a narrow Gaussian-like ring centered at a=1 kpc
        a = 1.0
        width = 0.05  # kpc
        r_profile = np.linspace(a - 3 * width, a + 3 * width, 50)
        r_profile = r_profile[r_profile > 0]
        sigma_profile = 100.0 * np.exp(-0.5 * ((r_profile - a) / width) ** 2)

        # Evaluate far from the ring where V -> sqrt(G*M_ring / R)
        r_eval = np.array([10.0, 20.0, 50.0])

        v = circular_velocity_thin_disk(r_eval, r_profile, sigma_profile)

        # At large R, V should fall off as ~1/sqrt(R)
        # V(R1)/V(R2) ~ sqrt(R2/R1)
        ratio_v = v[0] / v[1]
        ratio_expected = np.sqrt(r_eval[1] / r_eval[0])
        assert ratio_v == pytest.approx(ratio_expected, rel=0.05)

    def test_exponential_disk_freeman(self):
        """Compare against Freeman (1970) analytic result for an exponential disk.

        V^2(R) = 4*pi*G*Sigma_0*R_d * y^2 * [I_0(y)*K_0(y) - I_1(y)*K_1(y)]
        where y = R / (2*R_d).
        """
        R_d_kpc = 2.0  # disk scale length
        Sigma_0 = 50.0  # central surface density M_sun/pc^2

        # Create a finely-sampled exponential profile
        r_profile = np.linspace(0.05, 12 * R_d_kpc, 500)
        sigma_profile = Sigma_0 * np.exp(-r_profile / R_d_kpc)

        # Evaluate at several radii (avoid center where resolution matters)
        r_eval = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0])

        v_numerical = circular_velocity_thin_disk(r_eval, r_profile, sigma_profile)

        # Freeman (1970) analytic formula
        R_d_pc = R_d_kpc * KPC_TO_PC
        Sigma_0_pc = Sigma_0  # already in M_sun/pc^2
        y = (r_eval * KPC_TO_PC) / (2.0 * R_d_pc)
        v2_analytic = (
            4.0 * np.pi * G_PC * Sigma_0_pc * R_d_pc
            * y ** 2
            * (i0(y) * k0(y) - i1(y) * k1(y))
        )
        v_analytic = np.sqrt(np.maximum(v2_analytic, 0.0))

        # Should agree within 5% at all evaluation radii
        for i in range(len(r_eval)):
            if v_analytic[i] > 1.0:  # skip near-zero values
                assert v_numerical[i] == pytest.approx(v_analytic[i], rel=0.05), (
                    f"Mismatch at R={r_eval[i]} kpc: "
                    f"numerical={v_numerical[i]:.2f}, analytic={v_analytic[i]:.2f}"
                )

    def test_helium_factor_scaling(self):
        """V with helium_factor should scale as sqrt(factor) relative to factor=1."""
        r_profile = np.linspace(0.5, 10.0, 30)
        sigma_profile = 5.0 * np.exp(-r_profile / 3.0)
        r_eval = np.array([2.0, 5.0, 8.0])

        v_bare = circular_velocity_thin_disk(r_eval, r_profile, sigma_profile, helium_factor=1.0)
        v_he = circular_velocity_thin_disk(r_eval, r_profile, sigma_profile, helium_factor=1.33)

        expected_ratio = np.sqrt(1.33)
        for i in range(len(r_eval)):
            if v_bare[i] > 1.0:
                actual_ratio = v_he[i] / v_bare[i]
                assert actual_ratio == pytest.approx(expected_ratio, rel=0.01)

    def test_output_units_reasonable(self):
        """For M33-like surface densities, velocities should be 10-100 km/s."""
        r_profile = np.linspace(0.5, 15.0, 50)
        # M33-like gas profile: peaks ~7 M_sun/pc^2, declines outward
        sigma_profile = 7.0 * np.exp(-r_profile / 5.0)
        r_eval = np.array([1.0, 3.0, 5.0, 10.0])

        v = circular_velocity_thin_disk(r_eval, r_profile, sigma_profile, helium_factor=1.33)

        assert np.all(v > 0), "All velocities should be positive"
        assert np.all(v < 200), "Velocities should be physical (< 200 km/s for a gas disk)"
        assert np.max(v) > 10, "Peak velocity should exceed 10 km/s for M33-like gas"

    def test_returns_nonnegative(self):
        """Output should always be non-negative."""
        r_profile = np.array([1.0, 2.0, 3.0])
        sigma_profile = np.array([10.0, 5.0, 1.0])
        r_eval = np.linspace(0.1, 5.0, 20)

        v = circular_velocity_thin_disk(r_eval, r_profile, sigma_profile)
        assert np.all(v >= 0)


class TestFitOmega:
    def test_recovers_known_omega(self, synthetic_rotation_curve):
        """Fit synthetic data and verify omega is recovered within uncertainty."""
        data = synthetic_rotation_curve
        result = fit_omega(
            data["radius"],
            data["v_obs"],
            data["v_err"],
            data["v_bary"],
            galaxy_id="SYNTHETIC",
            upsilon_disk=1.0,
            upsilon_bulge=1.0,
        )
        assert result.converged
        assert result.omega_value == pytest.approx(data["true_omega"], abs=1.0)

    def test_zero_omega_flat_curve(self):
        """If V_obs == V_bary, fitted omega should be approximately zero."""
        radius = np.linspace(1.0, 10.0, 20)
        v_bary = 100.0 * np.ones_like(radius)
        v_obs = v_bary.copy()
        v_err = np.full_like(radius, 2.0)

        result = fit_omega(radius, v_obs, v_err, v_bary, galaxy_id="FLAT")
        assert result.converged
        assert result.omega_value == pytest.approx(0.0, abs=0.5)

    def test_positive_omega_detected(self):
        """Rising curve beyond V_bary should give positive omega."""
        radius = np.linspace(1.0, 10.0, 20)
        v_bary = 80.0 * np.ones_like(radius)
        omega_true = 3.0
        v_obs = v_bary + omega_true * radius
        v_err = np.full_like(radius, 1.0)

        result = fit_omega(radius, v_obs, v_err, v_bary)
        assert result.converged
        assert result.omega_value == pytest.approx(omega_true, abs=0.1)

    def test_flag_vobs_lt_vbary(self):
        """Should flag when V_obs < V_bary at any point."""
        radius = np.array([1.0, 2.0, 3.0])
        v_bary = np.array([50.0, 60.0, 70.0])
        v_obs = np.array([45.0, 65.0, 75.0])  # First point: V_obs < V_bary
        v_err = np.array([3.0, 3.0, 3.0])

        result = fit_omega(radius, v_obs, v_err, v_bary)
        assert result.flag_v_obs_lt_v_bary is True

    def test_no_flag_when_vobs_gt_vbary(self):
        """Should not flag when V_obs >= V_bary everywhere."""
        radius = np.array([1.0, 2.0, 3.0])
        v_bary = np.array([50.0, 60.0, 70.0])
        v_obs = np.array([55.0, 65.0, 75.0])
        v_err = np.array([3.0, 3.0, 3.0])

        result = fit_omega(radius, v_obs, v_err, v_bary)
        assert result.flag_v_obs_lt_v_bary is False

    def test_zero_errors_handled(self):
        """Zero errors should be replaced, not cause division by zero."""
        radius = np.array([1.0, 2.0, 3.0])
        v_bary = np.array([50.0, 60.0, 70.0])
        v_obs = np.array([55.0, 70.0, 85.0])
        v_err = np.array([0.0, 3.0, 3.0])  # First error is zero

        result = fit_omega(radius, v_obs, v_err, v_bary)
        assert result.converged
        assert np.isfinite(result.omega_value)

    def test_to_dict_excludes_arrays(self, synthetic_rotation_curve):
        """to_dict() should not contain numpy arrays."""
        data = synthetic_rotation_curve
        result = fit_omega(
            data["radius"], data["v_obs"], data["v_err"], data["v_bary"]
        )
        d = result.to_dict()
        for key, val in d.items():
            assert not isinstance(val, np.ndarray), f"Array found in to_dict(): {key}"

    def test_rmse_is_nonnegative(self, synthetic_rotation_curve):
        data = synthetic_rotation_curve
        result = fit_omega(
            data["radius"], data["v_obs"], data["v_err"], data["v_bary"]
        )
        assert result.residuals_rmse >= 0.0

    def test_reduced_chi2_reasonable(self, synthetic_rotation_curve):
        """For well-fitting data, reduced chi2 should be near 1."""
        data = synthetic_rotation_curve
        result = fit_omega(
            data["radius"], data["v_obs"], data["v_err"], data["v_bary"]
        )
        # With noise of 2 km/s and error bars of 3 km/s, reduced chi2 < 5 is reasonable
        assert result.reduced_chi_squared < 5.0
