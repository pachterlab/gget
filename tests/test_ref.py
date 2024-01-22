import unittest
import json
from gget.gget_ref import ref

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_ref.json") as json_file:
    ref_dict = json.load(json_file)


class TestRef(unittest.TestCase):
    def test_ref(self):
        test = "test1"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)

    def test_ref_which(self):
        test = "test2"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)

    def test_ref_rel(self):
        test = "test3"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)

    def test_ref_rel_ftp(self):
        test = "test4"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_ref_ftp(self):
        test = "test5"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_ref_plant(self):
        test = "test6"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)

    def test_ref_which_plant(self):
        test = "test7"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)

    def test_ref_rel_protist(self):
        test = "test8"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertEqual(result_to_test, expected_result)

    def test_ref_rel_ftp_octopus(self):
        test = "test9"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    def test_ref_list(self):
        test = "test10"
        expected_result = ref_dict[test]["expected_result"]
        result_to_test = ref(**ref_dict[test]["args"])

        self.assertListEqual(result_to_test, expected_result)

    ## Test bad input errors
    def test_ref_bad_species(self):
        test = "error_test1"
        with self.assertRaises(RuntimeError):
            ref(**ref_dict[test]["args"])

    def test_ref_bad_which(self):
        test = "error_test2"
        with self.assertRaises(ValueError):
            ref(**ref_dict[test]["args"])

    def test_ref_bad_rel(self):
        test = "error_test3"
        with self.assertRaises(RuntimeError):
            ref(**ref_dict[test]["args"])
