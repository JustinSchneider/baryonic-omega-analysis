"""End-to-end integration tests for the full pipeline."""

import pytest
import numpy as np

from src.database import (
    Base,
    get_session,
    init_db,
    query_profiles_as_dataframe,
)
from src.ingest import ingest_sparc_file
from src.physics import compute_v_bary, fit_omega


class TestEndToEnd:
    def test_ingest_and_fit(self, sample_sparc_file, tmp_path):
        """Full pipeline: parse file -> DB -> query -> fit -> verify."""
        db_path = str(tmp_path / "test.db")

        # Ingest
        galaxy_id = ingest_sparc_file(
            sample_sparc_file,
            upsilon_disk=0.5,
            upsilon_bulge=0.7,
            db_path=db_path,
        )
        assert galaxy_id == "TEST001"

        # Query back from DB
        engine = init_db(db_path)
        session = get_session(engine)
        df = query_profiles_as_dataframe(session, galaxy_id)
        session.close()

        assert len(df) == 5
        assert df.iloc[0]["radius_kpc"] == pytest.approx(0.5)

        # Compute V_bary and fit
        v_bary = compute_v_bary(
            df["v_gas"].values,
            df["v_disk"].values,
            df["v_bulge"].values,
            upsilon_disk=0.5,
            upsilon_bulge=0.7,
        )

        result = fit_omega(
            df["radius_kpc"].values,
            df["v_obs"].values,
            df["v_err"].values,
            v_bary,
            galaxy_id=galaxy_id,
        )

        assert result.converged
        assert np.isfinite(result.omega_value)
        assert result.n_points == 5

    def test_ingest_round_trip(self, sample_sparc_file, tmp_path):
        """Verify ingested data matches original file values."""
        from src.ingest import parse_sparc_rotmod

        db_path = str(tmp_path / "test_rt.db")

        # Parse original
        original = parse_sparc_rotmod(sample_sparc_file)

        # Ingest into DB
        galaxy_id = ingest_sparc_file(sample_sparc_file, db_path=db_path)

        # Query back
        engine = init_db(db_path)
        session = get_session(engine)
        from_db = query_profiles_as_dataframe(session, galaxy_id)
        session.close()

        # Compare key columns
        np.testing.assert_allclose(
            from_db["radius_kpc"].values, original["Rad"].values, rtol=1e-10,
        )
        np.testing.assert_allclose(
            from_db["v_obs"].values, original["Vobs"].values, rtol=1e-10,
        )
        np.testing.assert_allclose(
            from_db["v_gas"].values, original["Vgas"].values, rtol=1e-10,
        )
        np.testing.assert_allclose(
            from_db["v_disk"].values, original["Vdisk"].values, rtol=1e-10,
        )

    def test_negative_vgas_pipeline(self, negative_vgas_file, tmp_path):
        """Pipeline handles negative V_gas without errors."""
        db_path = str(tmp_path / "test_neg.db")

        galaxy_id = ingest_sparc_file(negative_vgas_file, db_path=db_path)

        engine = init_db(db_path)
        session = get_session(engine)
        df = query_profiles_as_dataframe(session, galaxy_id)
        session.close()

        # Verify negative value was preserved
        assert df.iloc[0]["v_gas"] == pytest.approx(-5.0)

        # Verify V_bary is still non-negative
        v_bary = compute_v_bary(
            df["v_gas"].values, df["v_disk"].values, df["v_bulge"].values,
        )
        assert np.all(v_bary >= 0)
