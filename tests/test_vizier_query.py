"""Tests for VizieR query utilities (mocked — no network access)."""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

from src.vizier_query import (
    search_vizier_catalog,
    list_vizier_tables,
    query_corbelli_2000,
    query_corbelli_2014,
    query_things_rotation_curves,
)


def _make_astropy_table(columns: dict, meta: dict | None = None):
    """Create a mock astropy Table-like object that supports to_pandas()."""
    mock_table = MagicMock()
    mock_table.colnames = list(columns.keys())
    mock_table.to_pandas.return_value = pd.DataFrame(columns)
    mock_table.meta = meta or {}
    return mock_table


class TestSearchVizierCatalog:
    @patch("src.vizier_query.Vizier")
    def test_returns_dataframe_on_success(self, mock_vizier_cls):
        table = _make_astropy_table({"Rad": [1.0, 2.0], "Vrot": [50.0, 80.0]})
        mock_vizier_cls.return_value.get_catalogs.return_value = [table]

        result = search_vizier_catalog("J/MNRAS/311/441")

        assert result is not None
        assert len(result) == 2
        assert "Rad" in result.columns

    @patch("src.vizier_query.Vizier")
    def test_returns_none_when_not_found(self, mock_vizier_cls):
        mock_vizier_cls.return_value.get_catalogs.return_value = []

        result = search_vizier_catalog("J/FAKE/000/000")
        assert result is None

    @patch("src.vizier_query.Vizier")
    def test_returns_none_on_exception(self, mock_vizier_cls):
        mock_vizier_cls.return_value.get_catalogs.side_effect = Exception("timeout")

        result = search_vizier_catalog("J/MNRAS/311/441")
        assert result is None

    @patch("src.vizier_query.Vizier")
    def test_respects_table_index(self, mock_vizier_cls):
        table0 = _make_astropy_table({"metadata": ["a"]})
        table1 = _make_astropy_table({"Rad": [1.0], "Vrot": [50.0]})
        mock_vizier_cls.return_value.get_catalogs.return_value = [table0, table1]

        result = search_vizier_catalog("J/MNRAS/311/441", table_index=1)
        assert result is not None
        assert "Rad" in result.columns

    @patch("src.vizier_query.Vizier")
    def test_invalid_table_index_returns_none(self, mock_vizier_cls):
        table = _make_astropy_table({"Rad": [1.0]})
        mock_vizier_cls.return_value.get_catalogs.return_value = [table]

        result = search_vizier_catalog("J/MNRAS/311/441", table_index=5)
        assert result is None


class TestListVizierTables:
    @patch("src.vizier_query.Vizier")
    def test_lists_table_descriptions(self, mock_vizier_cls):
        table = _make_astropy_table(
            {"Rad": [1.0]},
            meta={"description": "Rotation curve"},
        )
        mock_vizier_cls.return_value.get_catalogs.return_value = [table]

        result = list_vizier_tables("J/MNRAS/311/441")
        assert result is not None
        assert len(result) == 1
        assert "Rotation curve" in result[0]

    @patch("src.vizier_query.Vizier")
    def test_returns_none_when_not_found(self, mock_vizier_cls):
        mock_vizier_cls.return_value.get_catalogs.return_value = []
        assert list_vizier_tables("J/FAKE/000/000") is None


class TestQueryCorbelli:
    @patch("src.vizier_query.search_vizier_catalog")
    def test_corbelli_2000_delegates(self, mock_search):
        mock_search.return_value = pd.DataFrame({"Rad": [1.0]})
        result = query_corbelli_2000()
        mock_search.assert_called_once_with("J/MNRAS/311/441")
        assert result is not None

    @patch("src.vizier_query.search_vizier_catalog")
    def test_corbelli_2014_delegates(self, mock_search):
        mock_search.return_value = pd.DataFrame({"R": [1.0]})
        result = query_corbelli_2014()
        mock_search.assert_called_once_with("J/A+A/572/A23")
        assert result is not None


class TestQueryThings:
    @patch("src.vizier_query.list_vizier_tables")
    @patch("src.vizier_query.search_vizier_catalog")
    def test_groups_by_galaxy_name(self, mock_search, mock_list):
        mock_list.return_value = ["[0] Rotation curves"]
        mock_search.return_value = pd.DataFrame({
            "Galaxy": ["NGC2403", "NGC2403", "NGC3198", "NGC3198", "NGC3198"],
            "Rad": [1.0, 2.0, 1.0, 2.0, 3.0],
            "Vrot": [50, 80, 60, 90, 100],
        })

        result = query_things_rotation_curves()

        assert result is not None
        assert len(result) == 2
        assert "NGC2403" in result
        assert "NGC3198" in result
        assert len(result["NGC2403"]) == 2
        assert len(result["NGC3198"]) == 3

    @patch("src.vizier_query.list_vizier_tables")
    def test_returns_none_when_catalog_missing(self, mock_list):
        mock_list.return_value = None
        result = query_things_rotation_curves()
        assert result is None

    @patch("src.vizier_query.list_vizier_tables")
    @patch("src.vizier_query.search_vizier_catalog")
    def test_returns_none_when_no_galaxy_column(self, mock_search, mock_list):
        mock_list.return_value = ["[0] Some table"]
        mock_search.return_value = pd.DataFrame({
            "Rad": [1.0], "Vrot": [50.0],
        })

        result = query_things_rotation_curves()
        assert result is None
