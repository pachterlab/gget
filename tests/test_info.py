import unittest
import pandas as pd
import json
from gget.gget_info import info

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_info.json") as json_file:
    info_dict = json.load(json_file)


class TestInfo(unittest.TestCase):
    maxDiff = None

    def test_info_WB_transcript(self):
        test = "test2"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_FB_gene(self):
        test = "test3"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene(self):
        test = "test4"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_transcript(self):
        test = "test6"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_mix(self):
        test = "test7"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_exon(self):
        test = "test8"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    # def test_info_pdb(self):
    #     test = "test9"
    #     expected_result = info_dict[test]["expected_result"]
    #     result_to_test = info(**info_dict[test]["args"])
    #     # If result is a DataFrame, convert to list
    #     if isinstance(result_to_test, pd.DataFrame):
    #         result_to_test = result_to_test.dropna(axis=1).values.tolist()

    #     self.assertListEqual(result_to_test, expected_result)

    def test_info_ncbifalse_uniprottrue(self):
        test = "test10"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_ncbitrue_uniprotfalse(self):
        test = "test11"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_ncbifalse_uniprotfalse(self):
        test = "test12"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_ensembl_only(self):
        test = "test13"
        expected_result = info_dict[test]["expected_result"]
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        self.assertListEqual(result_to_test, expected_result)

    def test_info_bad_id(self):
        test = "none_test1"
        result_to_test = info(**info_dict[test]["args"])
        self.assertIsNone(result_to_test, "Invalid argument return is not None.")

    # Expected result not part of the unittest dictionary because of the unittest.mock.ANY entries
    def test_info_WB_gene(self):
        test = "test1"
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        expected_result = [
            [
                "WBGene00043981",
                "Q5WRS0",
                "caenorhabditis_elegans",
                "WBcel235",
                "aaim-1",
                "T14E8.4",
                [],
                "Protein aaim-1",
                "Uncharacterized protein [Source:NCBI gene;Acc:3565421]",
                "(Microbial infection) Promotes infection by microsporidian pathogens such as N.parisii in the early larval stages of development (PubMed:34994689). Involved in ensuring the proper orientation and location of the spore proteins of N.parisii during intestinal cell invasion (PubMed:34994689) Plays a role in promoting resistance to bacterial pathogens such as P.aeruginosa by inhibiting bacterial intestinal colonization",
                ["Secreted"],
                "Gene",
                "protein_coding",
                "T14E8.4.1.",
                "X",
                -1,
                6559466,
                6562428,
                ["T14E8.4.1"],
                ["protein_coding"],
                [unittest.mock.ANY],
                [-1],
                [6559466],
                [6562428],
            ]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene_list_non_model(self):
        test = "test5"
        result_to_test = info(**info_dict[test]["args"])
        # If result is a DataFrame, convert to list
        if isinstance(result_to_test, pd.DataFrame):
            result_to_test = result_to_test.dropna(axis=1).values.tolist()

        expected_result = [
            [
                "ENSMMUG00000054106.1",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "lncRNA",
                "ENSMMUT00000080640.1",
                "8",
                1,
                64990191,
                65000159,
                ["ENSMMUT00000080640.1", "ENSMMUT00000100253.1"],
                ["lncRNA", "lncRNA"],
                [unittest.mock.ANY, unittest.mock.ANY],
                [1, 1],
                [64990191, 64993715],
                [65000159, 64999625],
            ],
            [
                "ENSMMUG00000053116.1",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "protein_coding",
                "ENSMMUT00000091015.1",
                "3",
                -1,
                111461994,
                111475279,
                ["ENSMMUT00000091015.1"],
                ["protein_coding"],
                [unittest.mock.ANY],
                [-1],
                [111461994],
                [111475279],
            ],
            [
                "ENSMMUG00000021246.4",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "protein_coding",
                "ENSMMUT00000029894.4",
                "2",
                -1,
                98646979,
                98755023,
                [
                    "ENSMMUT00000029893.4",
                    "ENSMMUT00000053619.2",
                    "ENSMMUT00000104418.1",
                    "ENSMMUT00000087615.1",
                    "ENSMMUT00000103912.1",
                    "ENSMMUT00000086824.1",
                    "ENSMMUT00000029894.4",
                    "ENSMMUT00000104481.1",
                    "ENSMMUT00000090481.1",
                    "ENSMMUT00000026408.4",
                ],
                [
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                ],
                [
                    "CCDC13-201",
                    "CCDC13-202",
                    "CCDC13-203",
                    "CCDC13-204",
                    "CCDC13-205",
                    "CCDC13-206",
                    "CCDC13-207",
                    "CCDC13-208",
                    "CCDC13-209",
                    "CCDC13-210",
                ],
                [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                [
                    98646979,
                    98646982,
                    98646984,
                    98646984,
                    98646985,
                    98662368,
                    98663259,
                    98690412,
                    98737900,
                    98738870,
                ],
                [
                    98655536,
                    98656674,
                    98688669,
                    98656968,
                    98656865,
                    98713982,
                    98746548,
                    98713982,
                    98755023,
                    98754684,
                ],
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)
