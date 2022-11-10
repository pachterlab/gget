import unittest
import pandas as pd
import json
from gget.gget_info import info

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_info.json") as json_file:
    info_dict = json.load(json_file)

class TestInfo(unittest.TestCase):
    maxDiff = None

    def test_info_WB_gene(self):
        test = "test1"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_WB_transcript(self):
        test = "test2"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_FB_gene(self):
        test = "test3"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene(self):
        test = "test4"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene_list_non_model(self):
        test = "test5"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_transcript(self):
        test = "test6"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_mix(self):
        test = "test7"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_exon(self):
        test = "test8"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_bad_id(self):
        test = "none_test1"
        result_to_test = info(**info_dict[test]["args"])
        self.assertIsNone(result_to_test, "Invalid argument return is not None.")
