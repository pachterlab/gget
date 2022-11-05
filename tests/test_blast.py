import unittest
import pandas as pd
import json
from gget.gget_blast import blast

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_blast.json") as json_file:
    blast_dict = json.load(json_file)

class TestBlast(unittest.TestCase):
    def test_blast_nt(self):
        test = "test1"
        expected_result = blast_dict[test]["expected_result"]
        result_to_test = blast(**blast_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blast_nt_json(self):
        test = "test2"
        expected_result = blast_dict[test]["expected_result"]
        result_to_test = blast(**blast_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    def test_blast_bad_seq(self):
        test = "error_test1"
        with self.assertRaises(ValueError):
            blast(**blast_dict[test]["args"])

    def test_blast_bad_fasta(self):
        test = "error_test2"
        with self.assertRaises(FileNotFoundError):
            blast(**blast_dict[test]["args"])

    def test_blast_bad_program(self):
        test = "error_test3"
        with self.assertRaises(ValueError):
            blast(**blast_dict[test]["args"])

    def test_blast_db_missing(self):
        test = "error_test4"
        with self.assertRaises(ValueError):
            blast(**blast_dict[test]["args"])

    def test_blast_bad_db1(self):
        test = "error_test5"
        with self.assertRaises(ValueError):
            blast(**blast_dict[test]["args"])

    def test_blast_bad_db2(self):
        test = "error_test6"
        with self.assertRaises(ValueError):
            blast(**blast_dict[test]["args"])
