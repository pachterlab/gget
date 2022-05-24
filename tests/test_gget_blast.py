import unittest
from gget.gget_blast import blast


class TestBlast(unittest.TestCase):
    def test_blast_nt(self):
        df = blast(
            "ATACTCAGTCACACAAGCCATAGCAGGAAACAGCGAGCTTGCAGCCTCACCGACGAGTCTCAACTAAAAGGGACTCCCGGAGCTAGGGGTGGGGACTCGGCCTCACACAGTGAGTGCCGG",
            limit=1,
        )
        result_to_test = df.values.tolist()
        expected_result = [
            [
                "PREDICTED: Homo sapiens uncharacterized LOC105373836 (LOC105373836), transcript variant X2, ncRNA",
                "Homo sapiens",
                "human",
                9606,
                222,
                222,
                "100%",
                5.999999999999999e-54,
                "100.00%",
                2137,
                "XR_923785.3",
            ]
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_blast_nt_json(self):
        result_to_test = blast(
            "CACACATCCGGTTCTTCCGGGAGCTAGGGG",
            json=True,
        )
        expected_result = {
            "0": {
                "Description": "PREDICTED: Mus caroli G protein subunit alpha i3 (Gnai3), mRNA",
                "Scientific Name": "Mus caroli",
                "Common Name": "Ryukyu mouse",
                "Taxid": 10089,
                "Max Score": 56.5,
                "Total Score": 56.5,
                "Query Cover": "100%",
                "E value": 2e-05,
                "Per. Ident": "100.00%",
                "Acc. Len": 3253,
                "Accession": "XM_021157518.1",
            },
            "1": {
                "Description": "Mus musculus guanine nucleotide binding protein (G protein), alpha inhibiting 3 (Gnai3), mRNA",
                "Scientific Name": "Mus musculus",
                "Common Name": "house mouse",
                "Taxid": 10090,
                "Max Score": 56.5,
                "Total Score": 56.5,
                "Query Cover": "100%",
                "E value": 2e-05,
                "Per. Ident": "100.00%",
                "Acc. Len": 3294,
                "Accession": "NM_010306.3",
            },
            "2": {
                "Description": "Mus musculus adult male lung cDNA, RIKEN full-length enriched library, clone:1200003G10 product:guanine nucleotide binding protein, alpha inhibiting 3, full insert sequence",
                "Scientific Name": "Mus musculus",
                "Common Name": "house mouse",
                "Taxid": 10090,
                "Max Score": 56.5,
                "Total Score": 56.5,
                "Query Cover": "100%",
                "E value": 2e-05,
                "Per. Ident": "100.00%",
                "Acc. Len": 3261,
                "Accession": "AK004566.1",
            },
            "3": {
                "Description": "Mouse DNA sequence from clone RP23-89M15 on chromosome 3, complete sequence",
                "Scientific Name": "Mus musculus",
                "Common Name": "house mouse",
                "Taxid": 10090,
                "Max Score": 56.5,
                "Total Score": 56.5,
                "Query Cover": "100%",
                "E value": 2e-05,
                "Per. Ident": "100.00%",
                "Acc. Len": 207265,
                "Accession": "AL671854.14",
            },
        }
        self.assertEqual(result_to_test, expected_result)

    ## Turning this test off to reduce the number of BLASTs to run
    # def test_blast_aa_txt(self):
    #     df = blast(
    #         "tests/fixtures/muscle_aa_test.txt",
    #         limit=1,
    #     )
    #     result_to_test = df.values.tolist()
    #     expected_result = [
    #         [
    #             "FUN14 domain-containing protein 1 isoform X2 [Taeniopygia guttata]",
    #             "Taeniopygia guttata",
    #             "zebra finch",
    #             59729,
    #             180,
    #             180,
    #             "100%",
    #             3.9999999999999997e-56,
    #             "100.00%",
    #             156,
    #             "XP_002190216.1",
    #         ]
    #     ]
    #     self.assertListEqual(result_to_test, expected_result)

    def test_blast_bad_seq(self):
        with self.assertRaises(ValueError):
            blast("BANANA123")

    def test_blast_bad_fasta(self):
        with self.assertRaises(FileNotFoundError):
            blast("banana.fa")

    def test_blast_bad_txt(self):
        with self.assertRaises(FileNotFoundError):
            blast("banana.txt")

    def test_blast_bad_fileformat(self):
        with self.assertRaises(ValueError):
            blast("tests/fixtures/muscle_nt_test.banana")

    def test_blast_bad_program(self):
        with self.assertRaises(ValueError):
            blast("tests/fixtures/muscle_aa_test.txt", limit=3, program="banana")

    def test_blast_db_missing(self):
        with self.assertRaises(ValueError):
            blast(
                "tests/fixtures/muscle_aa_test.txt",
                limit=3,
                program="blastn",
            )

    def test_blast_bad_db1(self):
        with self.assertRaises(ValueError):
            blast(
                "tests/fixtures/muscle_aa_test.txt",
                limit=3,
                program="blastp",
                database="banana",
            )

    def test_blast_bad_db2(self):
        with self.assertRaises(ValueError):
            blast("tests/fixtures/muscle_aa_test.txt", limit=3, database="banana")
