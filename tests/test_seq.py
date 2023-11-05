import unittest
import pandas as pd
import json
from gget.gget_seq import seq

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_seq.json") as json_file:
    seq_dict = json.load(json_file)


class TestSeq(unittest.TestCase):
    def test_seq_gene(self):
        test = "test1"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript_gene_WB(self):
        test = "test2"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript_transcript_WB(self):
        test = "test3"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript_gene_WB_iso(self):
        test = "test4"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript_gene(self):
        test = "test5"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript(self):
        test = "test6"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_gene_iso(self):
        test = "test7"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_gene_transcript_iso(self):
        test = "test8"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript_gene_iso(self):
        test = "test9"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_seq_transcript_transcript_iso(self):
        test = "test10"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)
    
    def test_seq_missing_uniprot_gene_name(self):
        test = "test11"
        expected_result = seq_dict[test]["expected_result"]
        result_to_test = seq(**seq_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)
