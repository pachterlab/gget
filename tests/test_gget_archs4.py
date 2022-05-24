import unittest
from gget.gget_archs4 import archs4

from .fixtures import (
    ARCHS4_FUNDC1,
    ARCHS4_FUNDC1_TISSUE,
    ARCHS4_FUNDC1_TISSUE_MOUSE,
    ARCHS4_FUNDC1_TISSUE_JSON,
)


class TestArchs4(unittest.TestCase):
    def test_archs4(self):
        df = archs4("FUNDC1")
        result_to_test = df.values.tolist()
        expected_result = ARCHS4_FUNDC1

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_mouse(self):
        df = archs4("FUNDC1", species="mouse")
        result_to_test = df.values.tolist()
        expected_result = ARCHS4_FUNDC1

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_mouse_json(self):
        result_to_test = archs4("FUNDC1", gene_count=5, species="mouse", json=True)
        expected_result = [
            {"gene_symbol": "VBP1", "pearson_correlation": 0.537270844},
            {"gene_symbol": "SNRPB2", "pearson_correlation": 0.5164885521},
            {"gene_symbol": "PCMT1", "pearson_correlation": 0.5021503568},
            {"gene_symbol": "MRPL47", "pearson_correlation": 0.5008688569},
            {"gene_symbol": "NFU1", "pearson_correlation": 0.498568058},
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue(self):
        df = archs4("fuNdC1", which="tissue")
        result_to_test = df.values.tolist()
        expected_result = ARCHS4_FUNDC1_TISSUE

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue_json(self):
        result_to_test = archs4("fuNdC1", which="tissue", json=True)
        expected_result = ARCHS4_FUNDC1_TISSUE_JSON

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_tissue_mouse(self):
        df = archs4("fuNdC1", which="tissue", species="mouse")
        result_to_test = df.values.tolist()
        expected_result = ARCHS4_FUNDC1_TISSUE_MOUSE

        self.assertListEqual(result_to_test, expected_result)

    def test_archs4_bad_gene(self):
        result = archs4("BANANA")
        self.assertIsNone(result, "Invalid gene result is not None.")

    def test_archs4_bad_gene_tissue(self):
        result = archs4("BANANA", which="tissue", species="mouse")
        self.assertIsNone(result, "Invalid gene result is not None.")

    def test_archs4_bad_which(self):
        with self.assertRaises(ValueError):
            archs4("OAS1", which="banana")

    def test_archs4_bad_species(self):
        with self.assertRaises(ValueError):
            archs4("OAS1", species="banana")
