"""Tests for database schema, insertion, and round-trip integrity."""

import pytest
import numpy as np
import pandas as pd

from src.database import (
    Galaxy,
    OmegaFit,
    RadialProfile,
    insert_galaxy,
    insert_omega_fit,
    insert_radial_profiles,
    query_profiles_as_dataframe,
)


class TestGalaxyOperations:
    def test_insert_and_query_galaxy(self, in_memory_session):
        galaxy = insert_galaxy(
            in_memory_session, "NGC5055",
            distance_mpc=9.2, inclination=55.0, quality_flag=1,
        )
        assert galaxy.galaxy_id == "NGC5055"
        assert galaxy.distance_mpc == pytest.approx(9.2)
        assert galaxy.quality_flag == 1

    def test_upsert_updates_existing(self, in_memory_session):
        insert_galaxy(in_memory_session, "NGC5055", distance_mpc=9.0)
        insert_galaxy(in_memory_session, "NGC5055", distance_mpc=9.2)

        result = in_memory_session.get(Galaxy, "NGC5055")
        assert result.distance_mpc == pytest.approx(9.2)


class TestRadialProfiles:
    def test_insert_and_query(self, in_memory_session):
        insert_galaxy(in_memory_session, "TEST001")

        df = pd.DataFrame({
            "radius_kpc": [0.5, 1.0, 2.0],
            "v_obs": [25.0, 50.0, 80.0],
            "v_err": [3.0, 4.0, 5.0],
            "v_gas": [10.0, 20.0, 25.0],
            "v_disk": [15.0, 30.0, 40.0],
            "v_bulge": [0.0, 5.0, 8.0],
            "v_baryon_total": [12.5, 25.0, 35.0],
        })
        n = insert_radial_profiles(in_memory_session, "TEST001", df)
        assert n == 3

        result = query_profiles_as_dataframe(in_memory_session, "TEST001")
        assert len(result) == 3
        assert list(result.columns) == [
            "radius_kpc", "v_obs", "v_err", "v_gas",
            "v_disk", "v_bulge", "v_baryon_total",
        ]

    def test_round_trip_float_precision(self, in_memory_session):
        """Critical test: ingest -> SQL -> Pandas must preserve float values."""
        insert_galaxy(in_memory_session, "PRECISION_TEST")

        # Use values that stress float precision
        original = pd.DataFrame({
            "radius_kpc": [0.123456789, 1.987654321, 3.141592653],
            "v_obs": [25.123456, 50.654321, 80.111111],
            "v_err": [3.333333, 4.444444, 5.555555],
            "v_gas": [10.101010, 20.202020, 25.252525],
            "v_disk": [15.151515, 30.303030, 40.404040],
            "v_bulge": [0.0, 5.050505, 8.080808],
        })

        insert_radial_profiles(in_memory_session, "PRECISION_TEST", original)
        result = query_profiles_as_dataframe(in_memory_session, "PRECISION_TEST")

        # SQLite stores IEEE 754 doubles; should be exact to ~15 sig figs
        for col in ["radius_kpc", "v_obs", "v_err", "v_gas", "v_disk", "v_bulge"]:
            np.testing.assert_allclose(
                result[col].values,
                original[col].values,
                rtol=1e-10,
                err_msg=f"Float precision lost for column: {col}",
            )

    def test_replace_on_reingest(self, in_memory_session):
        """Second ingestion should replace, not duplicate."""
        insert_galaxy(in_memory_session, "REPLACE_TEST")

        df1 = pd.DataFrame({
            "radius_kpc": [1.0, 2.0],
            "v_obs": [50.0, 80.0],
            "v_err": [3.0, 4.0],
            "v_gas": [10.0, 20.0],
            "v_disk": [30.0, 40.0],
            "v_bulge": [0.0, 0.0],
        })
        insert_radial_profiles(in_memory_session, "REPLACE_TEST", df1)

        df2 = pd.DataFrame({
            "radius_kpc": [1.0, 2.0, 3.0],
            "v_obs": [55.0, 85.0, 100.0],
            "v_err": [3.0, 4.0, 5.0],
            "v_gas": [12.0, 22.0, 28.0],
            "v_disk": [32.0, 42.0, 48.0],
            "v_bulge": [0.0, 0.0, 0.0],
        })
        insert_radial_profiles(in_memory_session, "REPLACE_TEST", df2)

        result = query_profiles_as_dataframe(in_memory_session, "REPLACE_TEST")
        assert len(result) == 3  # Should be 3 (from df2), not 5
        assert result.iloc[0]["v_obs"] == pytest.approx(55.0)

    def test_query_empty_galaxy(self, in_memory_session):
        result = query_profiles_as_dataframe(in_memory_session, "NONEXISTENT")
        assert result.empty


class TestDataSource:
    def test_insert_with_data_source(self, in_memory_session):
        galaxy = insert_galaxy(
            in_memory_session, "M33",
            distance_mpc=0.7, data_source="Corbelli2000_VizieR",
        )
        assert galaxy.data_source == "Corbelli2000_VizieR"

    def test_data_source_nullable(self, in_memory_session):
        galaxy = insert_galaxy(in_memory_session, "NGC1234")
        assert galaxy.data_source is None

    def test_upsert_preserves_data_source(self, in_memory_session):
        insert_galaxy(in_memory_session, "M33", data_source="manual")
        insert_galaxy(in_memory_session, "M33", data_source="Corbelli2000_VizieR")
        result = in_memory_session.get(Galaxy, "M33")
        assert result.data_source == "Corbelli2000_VizieR"


class TestOmegaFits:
    def test_insert_and_query(self, in_memory_session):
        insert_galaxy(in_memory_session, "FIT_TEST")

        fit = insert_omega_fit(in_memory_session, {
            "galaxy_id": "FIT_TEST",
            "omega_value": 3.456,
            "omega_uncertainty": 0.123,
            "chi_squared": 15.0,
            "reduced_chi_squared": 1.5,
            "residuals_rmse": 4.2,
            "method_version": "v1_fixed_ML",
            "upsilon_disk": 0.5,
            "upsilon_bulge": 0.7,
            "n_points": 20,
            "converged": True,
            "flag_v_obs_lt_v_bary": False,
        })

        assert fit.fit_id is not None
        assert fit.omega_value == pytest.approx(3.456)
        assert fit.method_version == "v1_fixed_ML"
        assert fit.converged is True
