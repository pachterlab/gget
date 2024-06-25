import unittest
import json
import hashlib
import pandas as pd
from gget.gget_opentargets import opentargets

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_opentargets.json") as json_file:
    ot_dict = json.load(json_file)


class TestOpenTargets(unittest.TestCase):

    def test_opentargets(self):
        test = "test1"
        expected_result = ot_dict[test]["expected_result"]
        result_to_test = opentargets(**ot_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertEqual(result_to_test, expected_result)

    def test_opentargets_no_limit(self):
        test = "test2"
        expected_result = ot_dict[test]["expected_result"]
        result_to_test = opentargets(**ot_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        result_to_test = json.dumps(result_to_test)
        result_to_test = hashlib.md5(result_to_test.encode()).hexdigest()

        self.assertEqual(result_to_test, expected_result)

    def test_opentargets_no_specified_resource(self):
        test = "test3"
        expected_result = ot_dict[test]["expected_result"]
        result_to_test = opentargets(**ot_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertEqual(result_to_test, expected_result)

    def test_opentargets_drugs(self):
        test = "test4"
        expected_result = ot_dict[test]["expected_result"]
        result_to_test = opentargets(**ot_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertEqual(result_to_test, expected_result)

    def test_opentargets_drugs_no_limit(self):
        test = "test5"
        expected_result = ot_dict[test]["expected_result"]
        result_to_test = opentargets(**ot_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        result_to_test = json.dumps(result_to_test)
        result_to_test = hashlib.md5(result_to_test.encode()).hexdigest()

        self.assertEqual(result_to_test, expected_result)

    ## Test bad input errors
    def test_opentargets_bad_resource(self):
        test = "error_test1"
        with self.assertRaises(ValueError):
            opentargets(**ot_dict[test]["args"])

    def test_opentargets_bad_limit(self):
        test = "error_test2"
        with self.assertRaises(RuntimeError):
            opentargets(**ot_dict[test]["args"])

    def test_opentargets_nonexistent_id(self):
        test = "error_test3"
        with self.assertRaises(ValueError):
            opentargets(**ot_dict[test]["args"])
