import pytest
from src.cosmic_cancer_mutation_fasta_creation import create_mutant_sequence, substitution_mutation, deletion_mutation, insertion_mutation, delins_mutation, duplication_mutation, inversion_mutation, unknown_mutation

@pytest.fixture
def alphabet_sequence():
    return 'ABCDEFG'

@pytest.fixture
def alphabet_sequence_long():
    return 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ'

def create_test_row(mutation_cds, sequence):
    return {
        'Mutation CDS': mutation_cds,
        'full_sequence': sequence
    }

def test_substitution(alphabet_sequence):
    test_row = create_test_row('c.3C>X', alphabet_sequence)
    result = create_mutant_sequence(test_row, substitution_mutation)
    assert result == "ABXDEFG", f"Expected ABXDEFG, got {result}"

def test_substitution_long(alphabet_sequence_long):
    test_row = create_test_row('c.35I>X', alphabet_sequence_long)
    result = create_mutant_sequence(test_row, substitution_mutation)
    assert result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHXJKLMNOPQRSTUVWXYZABCDEFGHIJKLM", f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHXJKLMNOPQRSTUVWXYZABCDEFGHIJKLM, got {result}"

def test_multi_deletion(alphabet_sequence):
    test_row = create_test_row('c.3_6del', alphabet_sequence)
    result = create_mutant_sequence(test_row, deletion_mutation)
    assert result == "ABG", f"Expected ABG, got {result}"

def test_multi_deletion_long(alphabet_sequence_long):
    test_row = create_test_row('c.35_37del', alphabet_sequence_long)
    result = create_mutant_sequence(test_row, deletion_mutation)
    assert result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHLMNOPQRSTUVWXYZABCDEFGHIJKLMNO", f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHLMNOPQRSTUVWXYZABCDEFGHIJKLMNO, got {result}"

def test_single_deletion(alphabet_sequence):
    test_row = create_test_row('c.3del', alphabet_sequence)
    result = create_mutant_sequence(test_row, deletion_mutation)
    assert result == "ABDEFG", f"Expected ABDEFG, got {result}"

def test_multi_delins(alphabet_sequence):
    test_row = create_test_row('c.3_6delinsXYZ', alphabet_sequence)
    result = create_mutant_sequence(test_row, delins_mutation)
    assert result == "ABXYZG", f"Expected ABXYZG, got {result}"

def test_multi_delins_long(alphabet_sequence_long):
    test_row = create_test_row('c.35_37delinsWXYZ', alphabet_sequence_long)
    result = create_mutant_sequence(test_row, delins_mutation)
    assert result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHWXYZLMNOPQRSTUVWXYZABCDEFGHIJKLMNO", f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHWXYZLMNOPQRSTUVWXYZABCDEFGHIJKLMNO, got {result}"

def test_single_delins(alphabet_sequence):
    test_row = create_test_row('c.3delinsXYZ', alphabet_sequence)
    result = create_mutant_sequence(test_row, delins_mutation)
    assert result == "ABXYZDEFG", f"Expected ABXYZDEFG, got {result}"

def test_ins(alphabet_sequence):
    test_row = create_test_row('c.3_4insXYZ', alphabet_sequence)
    result = create_mutant_sequence(test_row, insertion_mutation)
    assert result == "ABCXYZDEFG", f"Expected ABCXYZDEFG, got {result}"

def test_ins_long(alphabet_sequence_long):
    test_row = create_test_row('c.35_36insPPP', alphabet_sequence_long)
    result = create_mutant_sequence(test_row, insertion_mutation)
    assert result == "EFGHIJKLMNOPQRSTUVWXYZABCDEFGHIPPPJKLMNOPQRSTUVWXYZABCDEFGHIJKLM", f"Expected EFGHIJKLMNOPQRSTUVWXYZABCDEFGHIPPPJKLMNOPQRSTUVWXYZABCDEFGHIJKLM, got {result}"

def test_multi_dup(alphabet_sequence):
    test_row = create_test_row('c.3_5dup', alphabet_sequence)
    result = create_mutant_sequence(test_row, duplication_mutation)
    assert result == "ABCDECDEFG", f"Expected ABCDECDEFG, got {result}"

def test_multi_dup_long(alphabet_sequence_long):
    test_row = create_test_row('c.39_42dup', alphabet_sequence_long)
    result = create_mutant_sequence(test_row, duplication_mutation)
    assert result == "IJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST", f"Expected IJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST, got {result}"

def test_single_dup(alphabet_sequence):
    test_row = create_test_row('c.5dup', alphabet_sequence)
    result = create_mutant_sequence(test_row, duplication_mutation)
    assert result == "ABCDEEFG", f"Expected ABCDEEFG, got {result}"

def test_inv(alphabet_sequence):
    test_row = create_test_row('c.3_4inv', alphabet_sequence)
    result = create_mutant_sequence(test_row, inversion_mutation)
    assert result == "ABDCEFG", f"Expected ABDCEFG, got {result}"

def test_inv_long(alphabet_sequence_long):
    test_row = create_test_row('c.41_42inv', alphabet_sequence_long)
    result = create_mutant_sequence(test_row, inversion_mutation)
    assert result == "KLMNOPQRSTUVWXYZABCDEFGHIJKLMNPOQRSTUVWXYZABCDEFGHIJKLMNOPQRST", f"Expected KLMNOPQRSTUVWXYZABCDEFGHIJKLMNPOQRSTUVWXYZABCDEFGHIJKLMNOPQRST, got {result}"







