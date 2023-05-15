import unittest
import pandas as pd
import json
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
