import unittest
import pandas as pd
import json
import time
from gget.gget_elm import elm
from bs4 import BeautifulSoup

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_elm.json") as json_file:
    elm_dict = json.load(json_file)

# Sleep time (in seconds) between tests to prevent surpassing the server rate limit
sleep_time = 65

class Testelm(unittest.TestCase):
    def test_elm_aa_seq(self):
        test = "test1"
        expected_result = elm_dict[test]["expected_result"]
        time.sleep(sleep_time)
        result_to_test = elm(**elm_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            #replace \xa0 with a space.
            result_to_test.replace(u'\xa0',u' ', regex=True, inplace=True)
             # cast all values to str add astype
            result_to_test = result_to_test.astype(str).values.tolist()
        
        time.sleep(sleep_time)
        self.assertListEqual(result_to_test, expected_result)

    def test_elm_aa_seq_2(self):
        test = "test2"
        expected_result = elm_dict[test]["expected_result"]
        time.sleep(sleep_time)
        result_to_test = elm(**elm_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            #replace \xa0 with a space.
            result_to_test.replace(u'\xa0',u' ', regex=True, inplace=True)
             # cast all values to str add astype
            result_to_test = result_to_test.astype(str).values.tolist()
        
        time.sleep(sleep_time)
        self.assertListEqual(result_to_test, expected_result)


    def test_elm_uniprot_id(self):
        test = "test3"
        expected_result = elm_dict[test]["expected_result"]
        result_to_test = elm(**elm_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            #replace \xa0 with a space.
            result_to_test.replace(u'\xa0',u' ', regex=True, inplace=True)
             # cast all values to str add astype
            result_to_test = result_to_test.astype(str).values.tolist()
        
        time.sleep(sleep_time * 3 + 5)
        self.assertListEqual(result_to_test, expected_result)

    def test_elm_uniprot_id_2(self):
        test = "test4"
        expected_result = elm_dict[test]["expected_result"]
        result_to_test = elm(**elm_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            #replace \xa0 with a space.
            result_to_test.replace(u'\xa0',u' ', regex=True, inplace=True)
             # cast all values to str add astype
            result_to_test = result_to_test.astype(str).values.tolist()
        
        time.sleep(sleep_time * 3 + 5)
        self.assertListEqual(result_to_test, expected_result)

    def test_elm_bad_aa_seq(self):
        time.sleep(sleep_time)
        test = "test5"
        with self.assertRaises(ValueError):
            elm(**elm_dict[test]["args"])
        time.sleep(sleep_time * 3 + 5)

