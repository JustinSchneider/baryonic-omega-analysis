"""Tests for the SPARC data parser and ingestion pipeline."""

import pytest
import numpy as np
import pandas as pd

from src.ingest import (
    extract_galaxy_name_from_filename,
    parse_sparc_rotmod,
    parse_corbelli_vizier,
    create_m33_manual_data,
    load_m33_corbelli2014_data,
    parse_things_galaxy,
    SPARC_COLUMNS,
)
from src.physics import compute_v_bary


class TestParseSparc:
    def test_returns_correct_columns(self, sample_sparc_file):
        df = parse_sparc_rotmod(sample_sparc_file)
        expected = ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
        assert list(df.columns) == expected

    def test_skips_comment_lines(self, sample_sparc_file):
        df = parse_sparc_rotmod(sample_sparc_file)
        # File has 2 comment lines + 5 data lines
        assert len(df) == 5

    def test_values_are_float(self, sample_sparc_file):
        df = parse_sparc_rotmod(sample_sparc_file)
        for col in df.columns:
            assert df[col].dtype == np.float64

    def test_first_row_values(self, sample_sparc_file):
        df = parse_sparc_rotmod(sample_sparc_file)
        assert df.iloc[0]["Rad"] == pytest.approx(0.50)
        assert df.iloc[0]["Vobs"] == pytest.approx(25.0)
        assert df.iloc[0]["Vgas"] == pytest.approx(10.0)

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_sparc_rotmod(tmp_path / "nonexistent.dat")

    def test_negative_vgas_parsed(self, negative_vgas_file):
        df = parse_sparc_rotmod(negative_vgas_file)
        assert df.iloc[0]["Vgas"] == pytest.approx(-5.0)


class TestExtractGalaxyName:
    def test_standard_name(self):
        assert extract_galaxy_name_from_filename("NGC5055_rotmod.dat") == "NGC5055"

    def test_ugc_name(self):
        assert extract_galaxy_name_from_filename("UGC02885_rotmod.dat") == "UGC02885"

    def test_path_with_directory(self, tmp_path):
        path = tmp_path / "NGC5055_rotmod.dat"
        assert extract_galaxy_name_from_filename(path) == "NGC5055"

    def test_name_without_rotmod(self):
        assert extract_galaxy_name_from_filename("M33.dat") == "M33"


class TestComputeVBaryon:
    def test_all_positive(self):
        v_gas = np.array([10.0])
        v_disk = np.array([20.0])
        v_bulge = np.array([0.0])
        result = compute_v_bary(v_gas, v_disk, v_bulge, upsilon_disk=0.5, upsilon_bulge=0.7)
        # V^2 = 10*10 + 0.5*20*20 + 0 = 100 + 200 = 300
        expected = np.sqrt(300.0)
        assert result[0] == pytest.approx(expected)

    def test_negative_vgas_reduces_vbary(self):
        v_gas_pos = np.array([10.0])
        v_gas_neg = np.array([-10.0])
        v_disk = np.array([20.0])
        v_bulge = np.array([0.0])

        v_bary_pos = compute_v_bary(v_gas_pos, v_disk, v_bulge)
        v_bary_neg = compute_v_bary(v_gas_neg, v_disk, v_bulge)

        # Negative V_gas should reduce V_bary
        assert v_bary_neg[0] < v_bary_pos[0]

    def test_zero_bulge_no_effect(self):
        v_gas = np.array([10.0])
        v_disk = np.array([20.0])
        v_bulge_zero = np.array([0.0])
        v_bulge_nonzero = np.array([5.0])

        result_zero = compute_v_bary(v_gas, v_disk, v_bulge_zero)
        result_nonzero = compute_v_bary(v_gas, v_disk, v_bulge_nonzero)

        assert result_zero[0] < result_nonzero[0]

    def test_upsilon_scaling(self):
        v_gas = np.array([0.0])
        v_disk = np.array([20.0])
        v_bulge = np.array([0.0])

        result_low = compute_v_bary(v_gas, v_disk, v_bulge, upsilon_disk=0.25)
        result_high = compute_v_bary(v_gas, v_disk, v_bulge, upsilon_disk=1.0)

        # Ratio should be sqrt(1.0/0.25) = 2.0
        assert (result_high[0] / result_low[0]) == pytest.approx(2.0)

    def test_output_always_nonnegative(self):
        # Extreme case: large negative V_gas dominating
        v_gas = np.array([-50.0])
        v_disk = np.array([5.0])
        v_bulge = np.array([0.0])
        result = compute_v_bary(v_gas, v_disk, v_bulge)
        assert result[0] >= 0.0


class TestParseCorbelli:
    def test_standard_vizier_columns(self):
        """Parse a DataFrame with typical VizieR column names."""
        df = pd.DataFrame({
            "Rmaj": [1.0, 2.0, 3.0],
            "Vrot": [50.0, 80.0, 100.0],
            "e_Vrot": [3.0, 4.0, 5.0],
            "Vgas": [10.0, 20.0, 25.0],
            "Vdisk": [30.0, 40.0, 45.0],
        })
        result = parse_corbelli_vizier(df)
        assert list(result.columns) == SPARC_COLUMNS
        assert result.iloc[0]["Rad"] == pytest.approx(1.0)
        assert result.iloc[1]["Vobs"] == pytest.approx(80.0)
        # M33 has no bulge
        assert (result["Vbul"] == 0.0).all()

    def test_alternative_column_names(self):
        """Parse with different but valid VizieR column names."""
        df = pd.DataFrame({
            "R": [1.0, 2.0],
            "Vc": [50.0, 80.0],
            "e_Vc": [3.0, 4.0],
            "VHI": [10.0, 20.0],
            "Vstar": [30.0, 40.0],
        })
        result = parse_corbelli_vizier(df)
        assert len(result) == 2
        assert result.iloc[0]["Vgas"] == pytest.approx(10.0)
        assert result.iloc[0]["Vdisk"] == pytest.approx(30.0)

    def test_missing_required_column_raises(self):
        """Raise ValueError when Rad/Vobs/errV cannot be mapped."""
        df = pd.DataFrame({
            "Vrot": [50.0, 80.0],
            "e_Vrot": [3.0, 4.0],
            # Missing radius column
        })
        with pytest.raises(ValueError, match="Rad"):
            parse_corbelli_vizier(df)

    def test_missing_vobs_raises(self):
        df = pd.DataFrame({
            "Rad": [1.0, 2.0],
            "e_Vrot": [3.0, 4.0],
            # Missing velocity column
        })
        with pytest.raises(ValueError, match="Vobs"):
            parse_corbelli_vizier(df)


class TestM33ManualData:
    def test_has_correct_columns(self):
        df = create_m33_manual_data()
        assert list(df.columns) == SPARC_COLUMNS

    def test_has_reasonable_row_count(self):
        df = create_m33_manual_data()
        assert 15 <= len(df) <= 25  # ~20 data points per Corbelli

    def test_radius_range(self):
        df = create_m33_manual_data()
        assert df["Rad"].min() <= 1.0
        assert df["Rad"].max() >= 15.0

    def test_vobs_range(self):
        """V_obs should rise to ~120-130 km/s per Corbelli."""
        df = create_m33_manual_data()
        assert df["Vobs"].max() >= 120.0
        assert df["Vobs"].max() <= 140.0

    def test_vgas_peaks_near_8kpc(self):
        """V_gas peaks ~40 km/s near R~8 kpc per Corbelli Section 4."""
        df = create_m33_manual_data()
        assert df["Vgas"].max() >= 35.0
        assert df["Vgas"].max() <= 45.0

    def test_no_bulge(self):
        df = create_m33_manual_data()
        assert (df["Vbul"] == 0.0).all()

    def test_all_values_nonnegative(self):
        df = create_m33_manual_data()
        for col in ["Rad", "Vobs", "errV", "Vgas", "Vdisk"]:
            assert (df[col] >= 0).all(), f"{col} has negative values"


class TestM33Corbelli2014Data:
    """Tests for loading Corbelli et al. (2014) Table 1 data with thin-disk conversion."""

    @pytest.fixture(autouse=True)
    def _load_data(self):
        self.df = load_m33_corbelli2014_data()

    def test_has_correct_columns(self):
        assert list(self.df.columns) == SPARC_COLUMNS

    def test_has_expected_row_count(self):
        assert 50 <= len(self.df) <= 65  # 57 rows in Table 1

    def test_radius_range(self):
        assert self.df["Rad"].min() <= 0.3
        assert self.df["Rad"].max() >= 22.0

    def test_vobs_range(self):
        """V_obs should rise to ~120-136 km/s per Corbelli Table 1."""
        assert self.df["Vobs"].max() >= 120.0
        assert self.df["Vobs"].max() <= 140.0

    def test_vgas_is_nonzero(self):
        """V_gas should be computed from surface densities, not zeros."""
        assert self.df["Vgas"].max() > 10.0

    def test_vdisk_is_nonzero(self):
        """V_disk should be computed from stellar surface density."""
        assert self.df["Vdisk"].max() > 10.0

    def test_no_bulge(self):
        assert (self.df["Vbul"] == 0.0).all()

    def test_all_values_nonnegative(self):
        for col in ["Rad", "Vobs", "errV", "Vgas", "Vdisk"]:
            assert (self.df[col] >= 0).all(), f"{col} has negative values"


class TestParseThingsGalaxy:
    def test_standard_things_format(self):
        df = pd.DataFrame({
            "Rad": [1.0, 2.0, 3.0],
            "Vobs": [50.0, 80.0, 100.0],
            "e_Vobs": [3.0, 4.0, 5.0],
            "Vgas": [10.0, 20.0, 25.0],
            "Vdisk": [30.0, 40.0, 45.0],
        })
        result = parse_things_galaxy(df, "NGC2403")
        assert list(result.columns) == SPARC_COLUMNS
        assert (result["Vbul"] == 0.0).all()

    def test_missing_required_raises(self):
        df = pd.DataFrame({
            "Rad": [1.0], "Vobs": [50.0],
            # Missing error column
        })
        with pytest.raises(ValueError, match="errV"):
            parse_things_galaxy(df, "NGC2403")

    def test_handles_alternative_gas_column(self):
        df = pd.DataFrame({
            "Rad": [1.0, 2.0],
            "Vrot": [50.0, 80.0],
            "e_Vrot": [3.0, 4.0],
            "VHI": [10.0, 20.0],
            "Vstar": [30.0, 40.0],
        })
        result = parse_things_galaxy(df, "NGC3198")
        assert result.iloc[0]["Vgas"] == pytest.approx(10.0)
        assert result.iloc[0]["Vdisk"] == pytest.approx(30.0)
