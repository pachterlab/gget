import unittest
import math
from gget.gget_info import info


class TestRef(unittest.TestCase):
    def test_info_gene(self):
        df = info("ENSTGUG00000006139")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1",
                ["FUNDC1"],
                "Uncharacterized protein",
                "FUN14 domain containing 1 [Source:NCBI gene;Acc:100228946]",
                "Gene",
                "protein_coding",
                "ENSTGUT00000027003.1",
                "1",
                -1,
                107513786,
                107528106,
            ]
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_gene_expand(self):
        df = info("ENSTGUG00000006139", expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1",
                ["FUNDC1"],
                "Uncharacterized protein",
                "FUN14 domain containing 1 [Source:NCBI gene;Acc:100228946]",
                "Gene",
                "protein_coding",
                "ENSTGUT00000027003.1",
                "1",
                -1,
                107513786,
                107528106,
                ["ENSTGUT00000006367", "ENSTGUT00000027003"],
                ["protein_coding", "protein_coding"],
                ["FUNDC1-201", "FUNDC1-202"],
            ]
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_transcript(self):
        df = info("ENSTGUT00000027003.1")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1-202",
                ["FUNDC1"],
                "ENSTGUG00000006139",
                "Uncharacterized protein",
                "Transcript",
                "protein_coding",
                "1",
                -1,
                107513786,
                107526965,
            ]
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_transcript_expand(self):
        df = info("ENSTGUT00000027003.1", expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1-202",
                ["FUNDC1"],
                "ENSTGUG00000006139",
                "Uncharacterized protein",
                "Transcript",
                "protein_coding",
                "1",
                -1,
                107513786,
                107526965,
                [
                    "ENSTGUEE00000179311",
                    "ENSTGUE00000062507",
                    "ENSTGUE00000062517",
                    "ENSTGUE00000062555",
                    "ENSTGUE00000062600",
                    "ENSTGUEE00000242055",
                ],
                [107526792, 107526298, 107524396, 107518467, 107515279, 107513786],
                [107526965, 107526454, 107524471, 107518595, 107516951, 107515185],
                [
                    "ENSTGUEE00000179311",
                    "ENSTGUE00000062507",
                    "ENSTGUE00000062517",
                    "ENSTGUE00000062555",
                    "ENSTGUE00000062600",
                    "ENSTGUEE00000242055",
                ],
                [107526792, 107526298, 107524396, 107518467, 107515279, 107513786],
                [107526965, 107526454, 107524471, 107518595, 107516951, 107515185],
            ]
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_mix(self):
        df = info(["ENSTGUT00000027003.1", "ENSTGUG00000006139"])
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1-202",
                ["FUNDC1"],
                "Uncharacterized protein",
                "Transcript",
                "protein_coding",
                "1",
                -1,
                107513786,
                107526965,
            ],
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1",
                ["FUNDC1"],
                "Uncharacterized protein",
                "Gene",
                "protein_coding",
                "1",
                -1,
                107513786,
                107528106,
            ],
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_exon(self):
        df = info(["ENSTGUEE00000179311"])
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "Exon",
                "1",
                -1,
                107526792,
                107526965,
            ]
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_exon_expand(self):
        df = info(["ENSTGUEE00000179311"], expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "Exon",
                "1",
                -1,
                107526792,
                107526965,
            ]
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_info_bad_id(self):
        result = info(["banana"])
        self.assertIsNone(result, "Invalid ID output is not None.")

    def test_info_bad_id_expand(self):
        result = info(["banana"], expand=True)
        self.assertIsNone(result, "Invalid ID output is not None.")
