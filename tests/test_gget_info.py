import unittest
from gget.gget_info import info


class TestRef(unittest.TestCase):
    def test_ref(self):
        df = info("ENSTGUG00000006139")

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "A0A674GVD2",
                "FUNDC1",
                "FUNDC1",
                ["FUNDC1"],
                math.nan,
                "Uncharacterized protein",
                math.nan,
                "FUN14 domain containing 1 [Source:NCBI gene;Acc:100228946]",
                "Gene",
                "protein_coding",
                "ENSTGUT00000027003.1",
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "1",
                -1,
                107513786,
                107528106,
            ]
        ]

        self.assertEqual(result_to_test, expected_result)
