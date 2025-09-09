import unittest
import pandas as pd
import json
import importlib.util
import sys
import types
from gget.gget_cellxgene import cellxgene

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_cellxgene.json") as json_file:
    cellxgene_dict = json.load(json_file)


def repr_dict(adata):
    """
    Function to convert the items/structure of an AnnData object to a dictionary.
    """
    d = {}
    for attr in (
        "n_obs",
        "n_vars",
        "obs",
        "var",
        "uns",
        "obsm",
        "varm",
        "layers",
        "obsp",
        "varp",
    ):
        got_attr = getattr(adata, attr)
        if isinstance(got_attr, int):
            d[attr] = got_attr
        else:
            keys = list(got_attr.keys())
            if keys:
                d[attr] = keys
    return d


cellxgene_census_spec = importlib.util.find_spec("cellxgene_census")


@unittest.skipUnless(cellxgene_census_spec, "cellxgene_census not installed")
class TestCellxgene(unittest.TestCase):
    def test_cellxgene_adata(self):
        test = "test1"
        expected_result = cellxgene_dict[test]["expected_result"]
        result_to_test = cellxgene(**cellxgene_dict[test]["args"])

        # Convert resulting AnnData object to dictionary
        result_to_test = repr_dict(result_to_test)

        self.assertEqual(result_to_test, expected_result)

    def test_cellxgene_metadata(self):
        test = "test2"
        expected_result = cellxgene_dict[test]["expected_result"]
        result_to_test = cellxgene(**cellxgene_dict[test]["args"])

        # Convert dataframe to list (and only keep first 25 results)
        result_to_test = result_to_test.values.tolist()[:25]

        self.assertListEqual(result_to_test, expected_result)


class TestCellxgeneGeneInput(unittest.TestCase):
    def setUp(self):
        class DummyCtx:
            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                pass

        def open_soma(*args, **kwargs):
            return DummyCtx()

        def get_anndata(*, census, organism, var_value_filter, obs_value_filter, column_names):
            return {"var_value_filter": var_value_filter}

        self._orig = sys.modules.get("cellxgene_census")
        sys.modules["cellxgene_census"] = types.SimpleNamespace(
            open_soma=open_soma, get_anndata=get_anndata
        )

    def tearDown(self):
        if self._orig is not None:
            sys.modules["cellxgene_census"] = self._orig
        else:
            del sys.modules["cellxgene_census"]

    def test_gene_as_string_human(self):
        result = cellxgene(gene="ACE2")
        self.assertEqual(result["var_value_filter"], "feature_name in ['ACE2']")

    def test_gene_as_string_mouse(self):
        result = cellxgene(gene="PAX7", species="mus_musculus")
        self.assertEqual(result["var_value_filter"], "feature_name in ['Pax7']")
