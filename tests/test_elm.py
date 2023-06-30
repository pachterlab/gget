import unittest
import pandas as pd
import json
import time
from gget.gget_elm import elm

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_elm.json") as json_file:
    elm_dict = json.load(json_file)

# Sleep time (in seconds) between tests to prevent surpassing the server rate limit
sleep_time = 120


class Testelm(unittest.TestCase):
    def test_elm_aa_seq(self):
        test = "test1"
        expected_result = elm_dict[test]["expected_result"]

        time.sleep(sleep_time)
        result_to_test = elm(**elm_dict[test]["args"])

        # replace \xa0 with a space.
        result_to_test.replace("\xa0", " ", regex=True, inplace=True)
        # cast all values to str add astype
        result_to_test = result_to_test.astype(str).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_elm_aa_seq_2(self):
        test = "test2"
        expected_result = elm_dict[test]["expected_result"]

        time.sleep(sleep_time * 3)
        result_to_test = elm(**elm_dict[test]["args"])

        # replace \xa0 with a space.
        result_to_test.replace("\xa0", " ", regex=True, inplace=True)
    
        # cast all values to str add astype
        result_to_test = result_to_test.astype(str).values.tolist()

        self.assertEqual(result_to_test, expected_result)

    def test_elm_uniprot_id(self):
        test = "test3"
        expected_result = set(elm_dict[test]["expected_result"])

        time.sleep(sleep_time * 3)
        result_to_test = elm(**elm_dict[test]["args"])

        # replace \xa0 with a space.
        result_to_test.replace("\xa0", " ", regex=True, inplace=True)
        # cast all values to str add astype
        result_to_test = result_to_test.astype(str).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_elm_bad_aa_seq(self):
        test = "test5"
        time.sleep(sleep_time)
        result_to_test = elm(**elm_dict[test]["args"])
        self.assertIsNone(result_to_test, "Bad AA sequence result is not None.")
