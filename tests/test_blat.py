import unittest
import pandas as pd
import json
from gget.gget_blat import blat

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_blat.json") as json_file:
    blat_dict = json.load(json_file)

class TestBlat(unittest.TestCase):
    def test_blat_nt(self):
        test = "test1"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_json(self):
        test = "test2"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_DNA(self):
        test = "test3"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_aa(self):
        test = "test4"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_aa_protein(self):
        test = "test5"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)
    
    def test_blat_nt_fasta(self):
        test = "test6"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_txt(self):
        test = "test7"
        expected_result = blat_dict[test]["expected_result"]
        result_to_test = blat(**blat_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_protein(self):
        test = "none_test1"
        result_to_test = blat(**blat_dict[test]["args"])
        self.assertIsNone(result_to_test, "DNA search in protein database is not None.")

    def test_blat_nt_RNA(self):
        test = "none_test2"
        result_to_test = blat(**blat_dict[test]["args"])
        self.assertIsNone(result_to_test, "DNA search in RNA database is not None.")

    def test_blat_nt_transDNA(self):
        test = "none_test3"
        result_to_test = blat(**blat_dict[test]["args"])
        self.assertIsNone(result_to_test, "DNA search in translated DNA database is not None.")

    def test_blat_aa_RNA(self):
        test = "none_test4"
        result_to_test = blat(**blat_dict[test]["args"])
        self.assertIsNone(result_to_test, "DNA search in RNA database is not None.")

    def test_blat_bad_assembly(self):
        test = "none_test5"
        result_to_test = blat(**blat_dict[test]["args"])
        self.assertIsNone(result_to_test, "Invalid assembly result is not None.")

    def test_blat_shortseq(self):
        test = "none_test6"
        result_to_test = blat(**blat_dict[test]["args"])
        self.assertIsNone(result_to_test, "Sequence too short result is not None.")

    def test_blat_bad_seqtype(self):
        test = "error_test1"
        with self.assertRaises(ValueError):
            blat(**blat_dict[test]["args"])

    def test_blat_bad_fileformat(self):
        test = "error_test2"
        with self.assertRaises(ValueError):
            blat(**blat_dict[test]["args"])

    def test_blat_bad_fasta(self):
        test = "error_test3"
        with self.assertRaises(FileNotFoundError):
            blat(**blat_dict[test]["args"])

    def test_blat_bad_txt(self):
        test = "error_test4"
        with self.assertRaises(FileNotFoundError):
            blat(**blat_dict[test]["args"])
