import unittest
import os
import pathlib as pl
from parameterized import parameterized
import pandas as pd
import json
import time
from gget.gget_cosmic import cosmic
from gget.utils import get_latest_cosmic
from .from_json import from_json

# Test download data from the latest COSMIC release
cosmic_version = get_latest_cosmic()

# List of argument combinations to test for the COSMIC database download
# Only example will work in the automatic tests since the others require COSMIC login data
arg_combinations = [
    ("cancer_example", 37, "CancerMutationCensus_AllData", cosmic_version),
    # ("cancer", 37, "CancerMutationCensus_AllData", cosmic_version),
    # ("cell_line", 38, "CellLinesProject_GenomeScreensMutant", cosmic_version),
    # ("census", 38, "Cosmic_MutantCensus", cosmic_version),
    # ("resistance", 38, "Cosmic_ResistanceMutations", cosmic_version),
    # ("genome_screen", 38, "Cosmic_GenomeScreensMutant", cosmic_version),
    # ("targeted_screen", 38, "Cosmic_CompleteTargetedScreensMutant", cosmic_version),
]


# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_cosmic.json") as json_file:
    cosmic_dict = json.load(json_file)

# Sleep time in seconds (wait [sleep_time] seconds between server requests to avoid 429 errors)
sleep_time = 10


class TestCosmic(unittest.TestCase, metaclass=from_json(cosmic_dict, cosmic, pre_test=lambda: time.sleep(sleep_time))):
    pass  # all tests are loaded from json


class TestCaseBase(unittest.TestCase):
    def assertIsFile(self, path):
        if not pl.Path(path).resolve().is_file():
            raise AssertionError("File does not exist: %s" % str(path))


class TestCosmicDownload(TestCaseBase):
    @parameterized.expand(arg_combinations)
    def test_cosmic_download(self, mut_class, grch_version, file_name, cosmic_version):
        time.sleep(sleep_time)

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
        gm_file = f"{file_name}_v{cosmic_version}_GRCh{grch_version}_mutation_workflow.csv"

        if mut_class == "cancer_example":
            tarred_folder = f"example_GRCh{grch_version}"

        path1 = pl.Path(os.path.abspath(os.path.join(tarred_folder, contained_file)))
        path2 = pl.Path(os.path.abspath(os.path.join(tarred_folder, gm_file)))

        self.assertIsFile(path1)
        self.assertIsFile(path2)

    def tearDown(self):
        # Get the folder name from the previous test's parameters
        test_case = self._testMethodName
        if "test_cosmic_download_" in test_case:
            mc = test_case.replace("test_cosmic_download_", "").split("_")[1:]
            mc = "_".join(mc)

            for arg_tuple in arg_combinations:
                if arg_tuple[0] == mc:
                    grch_version = arg_tuple[1]
                    file_name = arg_tuple[2]
                    cosmic_version = arg_tuple[3]

            if mc == "cancer_example":
                tarred_folder = f"example_GRCh{grch_version}"
            else:
                tarred_folder = f"{file_name}_Tsv_v{cosmic_version}_GRCh{grch_version}"

            # Remove tarred COSMIC database folder
            os.remove(tarred_folder + ".tar")

            # Delete untarred COSMIC database folder and all files within
            if os.path.exists(tarred_folder):
                for root, dirs, files in os.walk(tarred_folder, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                    for name in dirs:
                        os.rmdir(os.path.join(root, name))
                os.rmdir(tarred_folder)

        super(TestCosmicDownload, self).tearDown()
