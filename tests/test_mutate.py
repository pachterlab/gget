import unittest

from gget.gget_mutate import (
    create_mutant_sequence,
    substitution_mutation,
    deletion_mutation,
    insertion_mutation,
    delins_mutation,
    duplication_mutation,
    inversion_mutation,
)

# Length of flanking sequences
k = 30


# Test sequences
def alphabet_sequence():
    return "ABCDEFG"


def alphabet_sequence_long():
    return (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"
    )


def create_test_row(mutation, sequence):
    return {"mutation": mutation, "full_sequence": sequence}


class TestMutate(unittest.TestCase):
    def setUp(self):
        self.alphabet_sequence = alphabet_sequence()
        self.alphabet_sequence_long = alphabet_sequence_long()

    def test_substitution(self):
        test_row = create_test_row("c.3C>X", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, substitution_mutation, kmer_flanking_length=k
        )
        assert result == "ABXDEFG", f"Expected ABXDEFG, got {result}"

    def test_substitution_long(self):
        test_row = create_test_row("c.35I>X", self.alphabet_sequence_long)
        result = create_mutant_sequence(
            test_row, substitution_mutation, kmer_flanking_length=k
        )
        assert (
            result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHXJKLMNOPQRSTUVWXYZABCDEFGHIJKLM"
        ), f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHXJKLMNOPQRSTUVWXYZABCDEFGHIJKLM, got {result}"

    def test_multi_deletion(self):
        test_row = create_test_row("c.3_6del", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, deletion_mutation, kmer_flanking_length=k
        )
        assert result == "ABG", f"Expected ABG, got {result}"

    def test_multi_deletion_long(self):
        test_row = create_test_row("c.35_37del", self.alphabet_sequence_long)
        result = create_mutant_sequence(
            test_row, deletion_mutation, kmer_flanking_length=k
        )
        assert (
            result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHLMNOPQRSTUVWXYZABCDEFGHIJKLMNO"
        ), f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHLMNOPQRSTUVWXYZABCDEFGHIJKLMNO, got {result}"

    def test_single_deletion(self):
        test_row = create_test_row("c.3del", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, deletion_mutation, kmer_flanking_length=k
        )
        assert result == "ABDEFG", f"Expected ABDEFG, got {result}"

    def test_multi_delins(self):
        test_row = create_test_row("c.3_6delinsXYZ", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, delins_mutation, kmer_flanking_length=k
        )
        assert result == "ABXYZG", f"Expected ABXYZG, got {result}"

    def test_multi_delins_long(self):
        test_row = create_test_row("c.35_37delinsWXYZ", self.alphabet_sequence_long)
        result = create_mutant_sequence(
            test_row, delins_mutation, kmer_flanking_length=k
        )
        assert (
            result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHWXYZLMNOPQRSTUVWXYZABCDEFGHIJKLMNO"
        ), f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHWXYZLMNOPQRSTUVWXYZABCDEFGHIJKLMNO, got {result}"

    def test_single_delins(self):
        test_row = create_test_row("c.3delinsXYZ", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, delins_mutation, kmer_flanking_length=k
        )
        assert result == "ABXYZDEFG", f"Expected ABXYZDEFG, got {result}"

    def test_ins(self):
        test_row = create_test_row("c.3_4insXYZ", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, insertion_mutation, kmer_flanking_length=k
        )
        assert result == "ABCXYZDEFG", f"Expected ABCXYZDEFG, got {result}"

    def test_ins_long(self):
        test_row = create_test_row("c.35_36insPPP", self.alphabet_sequence_long)
        result = create_mutant_sequence(
            test_row, insertion_mutation, kmer_flanking_length=k
        )
        assert (
            result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHIPPPJKLMNOPQRSTUVWXYZABCDEFGHIJKLM"
        ), f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHIPPPJKLMNOPQRSTUVWXYZABCDEFGHIJKLM, got {result}"

    def test_multi_dup(self):
        test_row = create_test_row("c.3_5dup", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, duplication_mutation, kmer_flanking_length=k
        )
        assert result == "ABCDECDEFG", f"Expected ABCDECDEFG, got {result}"

    def test_multi_dup_long(self):
        test_row = create_test_row("c.39_42dup", self.alphabet_sequence_long)
        result = create_mutant_sequence(
            test_row, duplication_mutation, kmer_flanking_length=k
        )
        assert (
            result
            == "IJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST"
        ), f"Expected IJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST, got {result}"

    def test_single_dup(self):
        test_row = create_test_row("c.5dup", self.alphabet_sequence)
        result = create_mutant_sequence(
            test_row, duplication_mutation, kmer_flanking_length=k
        )
        assert result == "ABCDEEFG", f"Expected ABCDEEFG, got {result}"

    # def test_inv(self):
    #     test_row = create_test_row("c.3_4inv", self.alphabet_sequence)
    #     result = create_mutant_sequence(
    #         test_row, inversion_mutation, kmer_flanking_length=k
    #     )
    #     assert result == "ABDCEFG", f"Expected ABDCEFG, got {result}"

    # def test_inv_long(self):
    #     test_row = create_test_row("c.41_42inv", self.alphabet_sequence_long)
    #     result = create_mutant_sequence(
    #         test_row, inversion_mutation, kmer_flanking_length=k
    #     )
    #     assert (
    #         result == "KLMNOPQRSTUVWXYZABCDEFGHIJKLMNPOQRSTUVWXYZABCDEFGHIJKLMNOPQRST"
    #     ), f"Expected KLMNOPQRSTUVWXYZABCDEFGHIJKLMNPOQRSTUVWXYZABCDEFGHIJKLMNOPQRST, got {result}"
