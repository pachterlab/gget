import unittest
import os
import contextlib
import filecmp

from gget.gget_muscle import muscle


class TestMuscle(unittest.TestCase):
    def test_muscle_nt(self):
        # File with sequences to align
        fasta = "tests/fixtures/muscle_nt_test.fa"
        # File the results will be saved in
        out = "tests/fixtures/tmp1.afa"

        # Run muscle
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            muscle(fasta, out=out)

        # Expected result
        ref_path = "tests/fixtures/muscle_nt_test.afa"
        self.assertTrue(
            filecmp.cmp(out, ref_path, shallow=False),
            "The reference and muscle nucleotide alignment are not the same.",
        )

    def tearDown(self):
        super(TestMuscle, self).tearDown()
        # Delete temporary result file
        os.remove("tests/fixtures/tmp1.afa")


class TestMuscleSuper(unittest.TestCase):
    def test_muscle_nt_super5(self):
        # File with sequences to align
        fasta = "tests/fixtures/muscle_nt_test.fa"
        # File the results will be saved in
        out = "tests/fixtures/tmp2.afa"

        # Run muscle
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            muscle(fasta, out=out, super5=True)

        # Expected result
        ref_path = "tests/fixtures/muscle_super5_nt_test.afa"
        self.assertTrue(
            filecmp.cmp(out, ref_path, shallow=False),
            "The reference and muscle super5 nucleotide alignment are not the same.",
        )

    def tearDown(self):
        super(TestMuscleSuper, self).tearDown()
        # Delete temporary result file
        os.remove("tests/fixtures/tmp2.afa")


class TestMuscleAA(unittest.TestCase):
    def test_muscle_aa(self):
        # File with sequences to align
        fasta = "tests/fixtures/muscle_aa_test.fa"
        # File the results will be saved in
        out = "tests/fixtures/tmp3.afa"

        # Run muscle
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            muscle(fasta, out=out)

        # Expected result
        ref_path = "tests/fixtures/muscle_aa_test.afa"
        self.assertTrue(
            filecmp.cmp(out, ref_path, shallow=False),
            "The reference and muscle amino acid alignment are not the same.",
        )

    def tearDown(self):
        super(TestMuscleAA, self).tearDown()
        # Delete temporary result file
        os.remove("tests/fixtures/tmp3.afa")


class TestMuscleAASuper(unittest.TestCase):
    def test_muscle_aa_super5(self):
        # File with sequences to align
        fasta = "tests/fixtures/muscle_aa_test.fa"
        # File the results will be saved in
        out = "tests/fixtures/tmp4.afa"

        # Run muscle
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            muscle(fasta, out=out, super5=True)

        # Expected result
        ref_path = "tests/fixtures/muscle_super5_aa_test.afa"
        self.assertTrue(
            filecmp.cmp(out, ref_path, shallow=False),
            "The reference and muscle super5 amino acid alignment are not the same.",
        )

    def tearDown(self):
        super(TestMuscleAASuper, self).tearDown()
        # Delete temporary result file
        os.remove("tests/fixtures/tmp4.afa")
