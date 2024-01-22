import unittest
import pandas as pd
import json
from gget.gget_search import search

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_search.json") as json_file:
    search_dict = json.load(json_file)


class TestSearch(unittest.TestCase):
    def test_search_gene_one_sw(self):
        test = "test1"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])

        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_one_sw_json(self):
        test = "test2"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertEqual(result_to_test, expected_result)

    def test_search_gene_two_sw_or(self):
        test = "test3"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_one_sw_limit(self):
        test = "test4"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_transcript_one_sw(self):
        test = "test5"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_transcript_two_sw_or(self):
        test = "test6"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_transcript_two_sw_and(self):
        test = "test7"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_two_sw_limit(self):
        test = "test8"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_two_sw_and(self):
        test = "test9"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_release(self):
        test = "test10"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_db(self):
        test = "test11"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_plant(self):
        test = "test12"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_plant_db(self):
        test = "test13"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_search_octopus(self):
        test = "test14"
        expected_result = search_dict[test]["expected_result"]
        result_to_test = search(**search_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    ## Test bad input errors
    def test_search_gene_bad_species(self):
        test = "error_test1"
        with self.assertRaises(ValueError):
            search(**search_dict[test]["args"])

    def test_search_transcript_bad_species(self):
        test = "error_test2"
        with self.assertRaises(ValueError):
            search(**search_dict[test]["args"])

    def test_search_gene_bad_andor(self):
        test = "error_test3"
        with self.assertRaises(ValueError):
            search(**search_dict[test]["args"])

    def test_search_transcript_bad_andor(self):
        test = "error_test4"
        with self.assertRaises(ValueError):
            search(**search_dict[test]["args"])
