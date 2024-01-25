import unittest
import pandas as pd
import json
import time
from gget.gget_cosmic import cosmic

# Load
#  dictionary containing arguments and expected results
with open("./tests/fixtures/test_cosmic.json") as json_file:
    cosmic_dict = json.load(json_file)


class TestCosmic(unittest.TestCase):
    def test_cosmic_defaults(self):
        test = "test1"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    # Delay to avoid server overload errors
    time.sleep(60)
    def test_cosmic_limit_and_pubmet(self):
        test = "test2"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)
    def test_cosmic_json_and_genes(self):
        test = "test3"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)
    def test_cosmic_samples(self):
        test = "test4"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)
    def test_cosmic_studies(self):
        test = "test5"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)
    def test_cosmic_cancer(self):
        test = "test6"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(120)
    def test_cosmic_tumour(self):
        test = "test7"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)


# print(TestCosmic().test_cosmic_defaults())
# print(TestCosmic().test_cosmic_limit_and_pubmet())
# print(TestCosmic().test_cosmic_json_and_genes())
# print(TestCosmic().test_cosmic_samples())
# print(TestCosmic().test_cosmic_studies())
