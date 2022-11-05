import unittest
import pandas as pd
import json
from gget.gget_archs4 import archs4

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_archs4.json") as json_file:
    archs4_dict = json.load(json_file)


class TestArchs4(unittest.TestCase):
    def test_archs4_defaults(self):
        test = "test1"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_mouse(self):
        test = "test2"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_mouse_json(self):
        test = "test3"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_mouse_json_ensembl(self):
        test = "test4"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue(self):
        test = "test5"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue_json(self):
        test = "test6"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue_mouse(self):
        test = "test7"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue_ensembl(self):
        test = "test8"
        expected_result = archs4_dict[test]["expected_result"]
        result_to_test = archs4(**archs4_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_bad_ensembl(self):
        test = "none_test1"
        result_to_test = archs4(**archs4_dict[test]["args"])
        self.assertIsNone(result_to_test, "Invalid argument return is not None.")

    def test_archs4_bad_gene(self):
        test = "none_test2"
        result_to_test = archs4(**archs4_dict[test]["args"])
        self.assertIsNone(result_to_test, "Invalid argument return is not None.")

    def test_archs4_bad_gene_tissue(self):
        test = "none_test3"
        result_to_test = archs4(**archs4_dict[test]["args"])
        self.assertIsNone(result_to_test, "Invalid argument return is not None.")

    def test_archs4_bad_which(self):
        test = "error_test1"
        with self.assertRaises(ValueError):
            archs4(**archs4_dict[test]["args"])

    def test_archs4_bad_species(self):
        test = "error_test2"
        with self.assertRaises(ValueError):
            archs4(**archs4_dict[test]["args"])
