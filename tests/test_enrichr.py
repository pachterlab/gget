import unittest
import pandas as pd
import json
import matplotlib
import matplotlib.pyplot as plt
import math

# Prevent matplotlib from opening windows
matplotlib.use("Agg")
from gget.gget_enrichr import enrichr

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_enrichr.json") as json_file:
    enrichr_dict = json.load(json_file)


class TestEnrichr(unittest.TestCase):
    def test_enrichr_pathway(self):
        test = "test1"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_json(self):
        test = "test2"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_none(self):
        test = "test3"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_transcription(self):
        test = "test4"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_ensembl_ids(self):
        test = "test5"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_ontology(self):
        test = "test6"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_diseases_drugs(self):
        test = "test7"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_celltypes(self):
        test = "test8"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_kinase_interactions(self):
        test = "test9"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_bad_gene(self):
        test = "test10"
        df = enrichr(**enrichr_dict[test]["args"])

        self.assertTrue(df.empty, "Invalid gene result is not empty data frame.")

    def test_enrichr_background(self):
        test = "test11"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()[:20]
            result_to_test = [
                list(map(lambda x: x if x != math.inf else "inf", i))
                for i in result_to_test
            ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_background_ensembl(self):
        test = "test12"
        expected_result = enrichr_dict[test]["expected_result"]
        result_to_test = enrichr(**enrichr_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
            result_to_test = [
                list(map(lambda x: x if x != math.inf else "inf", i))
                for i in result_to_test
            ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_plot(self):
        # Number of plots before running enrichr plot
        num_figures_before = plt.gcf().number
        enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                "ASF1A",
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="transcription",
            plot=True,
        )
        # Number of plots after running enrichr plot
        num_figures_after = plt.gcf().number

        self.assertGreater(
            num_figures_after,
            num_figures_before,
            "No matplotlib plt object was created.",
        )
