import unittest
import os
import pathlib as pl
from parameterized import parameterized
import pandas as pd
import json
import time
from gget.gget_cosmic import cosmic
from gget.utils import get_latest_cosmic

# Test download data from the latest COSMIC release
cosmic_version = get_latest_cosmic()

# List of argument combinations to test for the COSMIC database download
# Only example will work in the automatic tests since the others require COSMIC login data
arg_combinations = [
    ("cancer_example", 37, "CancerMutationCensus", cosmic_version),
    # ("cancer", 37, "CancerMutationCensus_AllData", cosmic_version),
    # ("cell_line", 38, "CellLinesProject_GenomeScreensMutant", cosmic_version),
    # ("census", 38, "Cosmic_MutantCensus", cosmic_version),
    # ("resistance", 38, "Cosmic_ResistanceMutations", cosmic_version),
    # ("screen", 38, "Cosmic_GenomeScreensMutant", cosmic_version),
]


# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_cosmic.json") as json_file:
    cosmic_dict = json.load(json_file)


class TestCosmic(unittest.TestCase):
    def test_cosmic_defaults(self):
        test = "test1"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    # Delay to avoid server overload errors
    time.sleep(60)

    def test_cosmic_limit_and_pubmet(self):
        test = "test2"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)

    def test_cosmic_json_and_genes(self):
        test = "test3"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)

    def test_cosmic_samples(self):
        test = "test4"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)

    def test_cosmic_studies(self):
        test = "test5"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(60)

    def test_cosmic_cancer(self):
        test = "test6"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)

    time.sleep(120)

    def test_cosmic_tumour(self):
        test = "test7"
        expected_result = cosmic_dict[test]["expected_result"]
        result_to_test = cosmic(**cosmic_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.values.tolist()
        self.assertListEqual(result_to_test, expected_result)


class TestCaseBase(unittest.TestCase):
    def assertIsFile(self, path):
        if not pl.Path(path).resolve().is_file():
            raise AssertionError("File does not exist: %s" % str(path))


class TestCosmicDownload(TestCaseBase):
    @parameterized.expand(arg_combinations)
    def test_cosmic_download(self, mut_class, grch_version, file_name, cosmic_version):
        cosmic(
            searchterm=None,
            download_cosmic=True,
            mutation_class=mut_class,
            cosmic_version=cosmic_version,
            grch_version=grch_version,
            gget_mutate=True,
        )
        tarred_folder = f"{file_name}_Tsv_v{cosmic_version}_GRCh{grch_version}"
        contained_file = f"{file_name}_v{cosmic_version}_GRCh{grch_version}.tsv"

        if mut_class == "cancer_example":
            path = pl.Path(contained_file)
        else:
            path = pl.Path(os.path.join(tarred_folder, contained_file))

        self.assertIsFile(path)

    def tearDown(self):
        # Get the folder name from the previous test's parameters
        test_case = self._testMethodName
        if "test_cosmic_download_" in test_case:
            params = test_case.replace("test_cosmic_download_", "").split("_")
            grch_version = params[1]
            file_name = params[2]
            cosmic_version = params[3]

            tarred_folder = f"{file_name}_Tsv_v{cosmic_version}_GRCh{grch_version}"

            if os.path.exists(tarred_folder):
                for root, dirs, files in os.walk(tarred_folder, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                    for name in dirs:
                        os.rmdir(os.path.join(root, name))
                os.rmdir(tarred_folder)

        super(TestCosmicDownload, self).tearDown()
