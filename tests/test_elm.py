import unittest

# Library to test functions that have calls to print()
from unittest import mock
import os
import contextlib
import filecmp
import json

from gget.gget_elm import diamond, elm

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_elm.json") as json_file:
    elm_dict = json.load(json_file)


class TestELM(unittest.TestCase):
    def test_diamond(self):
        # File with sequences to align
        fasta = "tests/fixtures/input.fa"
        # File the results will be saved in
        out = "tests/fixtures/out.tsv"

        # Run muscle (trying to use contextlib here to silence muscle return (unsuccesfully))
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            diamond(fasta, out=out)

        # Expected result
        ref_path = "tests/fixtures/out.tsv"
        self.assertTrue(
            filecmp.cmp(out, ref_path, shallow=False),
            "The reference and diamond blast are not the same.",
        )
    
    def test_elm_uniprot_id(self):
        test = "test1"
        expected_result = elm_dict[test]["expected_result"]
        result_to_test = elm(**elm_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        time.sleep(sleep_time)
        self.assertListEqual(result_to_test, expected_result)
    
    def test_elm_seq(self):
        test = "test2"
        expected_result = elm_dict[test]["expected_result"]
        result_to_test = elm(**elm_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        time.sleep(sleep_time)
        self.assertListEqual(result_to_test, expected_result)
