import unittest
import json

from gget.gget_diamond import diamond

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_diamond.json") as json_file:
    diamond_dict = json.load(json_file)


class TestDiamond(unittest.TestCase):
    def test_diamond_seqs_multiple(self):
        test = "test1"
        expected_result = diamond_dict[test]["expected_result"]

        result_to_test = diamond(**diamond_dict[test]["args"])
        result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_diamond_seqs_single(self):
        test = "test2"
        expected_result = diamond_dict[test]["expected_result"]

        result_to_test = diamond(**diamond_dict[test]["args"])
        result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)
    
    def test_diamond_ref_file(self):
        test = "test3"
        expected_result = diamond_dict[test]["expected_result"]

        result_to_test = diamond(**diamond_dict[test]["args"])
        result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)
    
    def test_diamond_query_file(self):
        test = "test4"
        expected_result = diamond_dict[test]["expected_result"]

        result_to_test = diamond(**diamond_dict[test]["args"])
        result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)
    
    def test_diamond_both_files(self):
        test = "test5"
        expected_result = diamond_dict[test]["expected_result"]

        result_to_test = diamond(**diamond_dict[test]["args"])
        result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)
    
    def test_diamond_JSON_out(self):
        test = "test6"
        expected_result = diamond_dict[test]["expected_result"]

        result_to_test = diamond(**diamond_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)