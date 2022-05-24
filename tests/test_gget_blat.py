import unittest
from gget.gget_blat import blat


class TestBlat(unittest.TestCase):
    def test_blat_nt(self):
        df = blat("ATGCTGAATTTATGCTGAATTTATGCTGAATTTATGCTGAATTT")
        result_to_test = df.values.tolist()
        expected_result = [
            ["hg38", 44, 17, 44, 28, 0, 63.64, 100.0, "chr18", "+", 1952614, 1952647]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_json(self):
        result_to_test = blat(
            "CACACATCCGGTTCTTCCGGGAGCTAGGGG", assembly="mouse", json=True
        )
        expected_result = [
            {
                "genome": "mm39",
                "query_size": 30,
                "aligned_start": 1,
                "aligned_end": 30,
                "matches": 30,
                "mismatches": 0,
                "%_aligned": 100.0,
                "%_matched": 100.0,
                "chromosome": "chr3",
                "strand": "-",
                "start": 108053433,
                "end": 108053462,
            }
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_DNA(self):
        df = blat("ATGCTGAATTTATGCTGAATTTATGCTGAATTTATGCTGAATTT", seqtype="DNA")
        result_to_test = df.values.tolist()
        expected_result = [
            ["hg38", 44, 17, 44, 28, 0, 63.64, 100.0, "chr18", "+", 1952614, 1952647]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_protein(self):
        result = blat(
            "ATGCTGAATTTATGCTGAATTTATGCTGAATTTATGCTGAATTT",
            seqtype="protein",
            assembly="zebrafinch",
        )
        self.assertIsNone(result, "DNA search in protein database is not None.")

    def test_blat_nt_RNA(self):
        result = blat(
            "ATGCTGAATTTATGCTGAATTTATGCTGAATTTATGCTGAATTT", seqtype="translated%20RNA"
        )
        self.assertIsNone(result, "DNA search in RNA database is not None.")

    def test_blat_nt_transDNA(self):
        result = blat(
            "ATGCTGAATTTATGCTGAATTTATGCTGAATTTATGCTGAATTT", seqtype="translated%20DNA"
        )
        self.assertIsNone(result, "DNA search in translated DNA database is not None.")

    def test_blat_aa(self):
        df = blat(
            "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
            assembly="mouse",
        )
        result_to_test = df.values.tolist()
        expected_result = [
            ["mm39", 64, 25, 64, 35, 5, 62.5, 87.5, "chrX", "+-", 17437571, 17437690]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_blat_aa_protein(self):
        df = blat(
            "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
            seqtype="protein",
            assembly="mouse",
        )
        result_to_test = df.values.tolist()
        expected_result = [
            ["mm39", 64, 25, 64, 35, 5, 62.5, 87.5, "chrX", "+-", 17437571, 17437690]
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_blat_aa_DNA(self):
        result = blat(
            "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
            seqtype="DNA",
        )
        self.assertIsNone(result, "DNA search in protein database is not None.")

    def test_blat_aa_RNA(self):
        result = blat(
            "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
            seqtype="translated%20RNA",
        )
        self.assertIsNone(result, "DNA search in RNA database is not None.")

    def test_blat_aa_transDNA(self):
        result = blat(
            "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
            seqtype="translated%20DNA",
        )
        self.assertIsNone(result, "DNA search in translated DNA database is not None.")

    def test_blat_bad_seqtype(self):
        with self.assertRaises(ValueError):
            blat(
                "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
                seqtype="banana",
            )

    def test_blat_bad_assembly(self):
        result = blat(
            "MLMPGPLRRALGQKFSIFPSVDHDSDDDSYEVLDLTEYARRHHWWNRLFGRNSGPVVEKYSVAT",
            seqtype="DNA",
            assembly="banana",
        )
        self.assertIsNone(result, "Invalid assembly result is not None.")

    def test_blat_shortseq(self):
        result = blat("MLMPGPLRRALGQ")
        self.assertIsNone(result, "Sequence too short result is not None.")

    def test_blat_bad_fileformat(self):
        with self.assertRaises(ValueError):
            blat("tests/fixtures/muscle_nt_test.banana", assembly="zebrafinch")

    def test_blat_nt_fasta(self):
        df = blat("tests/fixtures/muscle_nt_test.fa", assembly="zebrafinch")
        result_to_test = df.values.tolist()
        expected_result = [
            ["taeGut2", 63, 1, 63, 63, 0, 100.0, 100.0, "chr1", "+", 5648870, 5648932],
            [
                "taeGut2",
                63,
                19,
                40,
                21,
                1,
                34.92,
                95.45,
                "chr27",
                "-",
                3975325,
                3975346,
            ],
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_blat_nt_txt(self):
        df = blat("tests/fixtures/muscle_nt_test.txt", assembly="zebrafinch")
        result_to_test = df.values.tolist()
        expected_result = [
            ["taeGut2", 63, 1, 63, 63, 0, 100.0, 100.0, "chr1", "+", 5648870, 5648932],
            [
                "taeGut2",
                63,
                19,
                40,
                21,
                1,
                34.92,
                95.45,
                "chr27",
                "-",
                3975325,
                3975346,
            ],
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_blat_bad_fasta(self):
        with self.assertRaises(FileNotFoundError):
            blat("banana.fa")

    def test_blat_bad_txt(self):
        with self.assertRaises(FileNotFoundError):
            blat("banana.txt")
