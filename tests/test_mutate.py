import json

import pytest
import unittest
import gget
import pandas as pd
import os
import tempfile
from typing import Callable
from .from_json import from_json, do_call

LONG_SEQUENCE = 'CCCCGCCCCACCCCGCCCCTCCCCGCCCCACCCCGCCCCTCCCCGCCCCACCCCGCCCCTCCCCGCCCCACCCCG'


@pytest.fixture
def long_sequence():
    return LONG_SEQUENCE


@pytest.fixture
def create_temp_files(long_sequence):
    # Create a temporary CSV file
    temp_csv_file = tempfile.NamedTemporaryFile(delete=False, suffix='.csv')
    
    # Data to write to CSV
    mutation_list = ["c.35G>A", "c.65G>A", "c.35del", "c.4_5insT"]
    protein_mutation_list = ['A12T', 'A22T', 'A12del', 'A4_5insT']
    mut_ID_list = ['GENE1_MUT1A_MUT1B', 'GENE1_MUT2A_MUT2B', 'GENE2_MUT1A_MUT1B', 'GENE3_MUT1A_MUT1B']
    seq_ID_list = ['ENST1', 'ENST2', 'ENST3', 'ENST4']
    
    data = {
        'mutation': mutation_list,
        'mutation_aa': protein_mutation_list,
        'mut_ID': mut_ID_list,
        'seq_ID': seq_ID_list
    }

    df = pd.DataFrame(data)
    df.to_csv(temp_csv_file.name, index=False)

    # Create a temporary FASTA file
    sequence_list = [long_sequence for _ in range(len(mutation_list))]
    temp_fasta_file = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta')
    
    with open(temp_fasta_file.name, 'w') as fasta_file:
        for seq_id, sequence in zip(seq_ID_list, sequence_list):
            fasta_file.write(f">{seq_id}\n")
            fasta_file.write(f"{sequence}\n")
    
    yield temp_csv_file.name, temp_fasta_file.name
    
    # Cleanup
    os.remove(temp_csv_file.name)
    os.remove(temp_fasta_file.name)


def assert_global_variables_zero(number_intronic_position_mutations = 0, number_posttranslational_region_mutations = 0, number_uncertain_mutations = 0, number_ambiguous_position_mutations = 0, number_index_errors = 0):
    assert gget.gget_mutate.intronic_mutations == number_intronic_position_mutations
    assert gget.gget_mutate.posttranslational_region_mutations == number_posttranslational_region_mutations
    assert gget.gget_mutate.uncertain_mutations == number_uncertain_mutations
    assert gget.gget_mutate.ambiguous_position_mutations == number_ambiguous_position_mutations
    assert gget.gget_mutate.mut_idx_outside_seq == number_index_errors


def _recursive_replace(v: str|dict|list, old: str, new: str, exact: bool = False) -> str|dict|list:
    if isinstance(v, str):
        if exact:
            if v == old:
                return new
            else:
                return v
        else:
            return v.replace(old, new)
    elif isinstance(v, dict):
        return {_recursive_replace(k, old, new, exact=exact): _recursive_replace(v, old, new, exact=exact) for k, v in v.items()}
    elif isinstance(v, list):
        return [_recursive_replace(v, old, new, exact=exact) for v in v]
    else:
        return v


def _assert_mutate(name: str, td: dict[str, dict[str, ...]], func: Callable) -> Callable:
    ls = LONG_SEQUENCE

    def assert_mutate(self: unittest.TestCase):
        test = name
        expected_result = td[test].get("expected_result", None)
        global_variables = td[test].get("global_variables", {})

        args = td[test]["args"]

        args = _recursive_replace(args, "<long_sequence>", ls, exact=True)

        result = do_call(func, args)

        if expected_result:
            self.assertEqual(result[0], expected_result)

        assert_global_variables_zero(**global_variables)

    return assert_mutate


with open("./tests/fixtures/test_mutate.json") as json_file:
    mutate_dict = json.load(json_file)


class TestMutate(unittest.TestCase, metaclass=from_json(mutate_dict, gget.mutate, {"assert_mutate": _assert_mutate})):
    pass  # all tests are loaded from json


def test_single_substitution(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.35G>A"
    )

    assert result[0] == "GCCCCACCCCGCCCCTCCCCGCCCCACCCCACCCCTCCCCGCCCCACCCCGCCCCTCCCCG"

    assert_global_variables_zero()

def test_single_substitution_near_right_end(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.65G>A"
    )

    assert result[0] == "GCCCCTCCCCGCCCCACCCCGCCCCTCCCCACCCCACCCCG"

    assert_global_variables_zero()


def test_single_substitution_near_left_end(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.5G>A"
    )

    assert result[0] == "CCCCACCCCACCCCGCCCCTCCCCGCCCCACCCCG"

    assert_global_variables_zero()


def test_single_deletion(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.35del"  # del the G
    )

    assert result[0] == "GCCCCACCCCGCCCCTCCCCGCCCCACCCCCCCCTCCCCGCCCCACCCCGCCCCTCCCCG"

    assert_global_variables_zero()


def test_multi_deletion(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.35_40del"
    )

    assert result[0] == "GCCCCACCCCGCCCCTCCCCGCCCCACCCCCCCCGCCCCACCCCGCCCCTCCCCGCCCCA"

    assert_global_variables_zero()

def test_single_deletion_with_right_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.31del"
    )

    assert result[0] == "CGCCCCACCCCGCCCCTCCCCGCCCCACCCGCCCCTCCCCGCCCCACCCCGCCCCTC"

    assert_global_variables_zero()

def test_single_deletion_with_left_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.34del"
    )

    assert result[0] == "CGCCCCACCCCGCCCCTCCCCGCCCCACCCGCCCCTCCCCGCCCCACCCCGCCCCTC"

    assert_global_variables_zero()

def test_multi_deletion_with_right_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.31_32del"
    )

    assert result[0] == "CCGCCCCACCCCGCCCCTCCCCGCCCCACCGCCCCTCCCCGCCCCACCCCGCCCCTCC"

    assert_global_variables_zero()

def test_single_insertion(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.4_5insT"
    )

    assert result[0] == "CCCCTGCCCCACCCCGCCCCTCCCCGCCCCACCCC"

    assert_global_variables_zero()


def test_multi_insertion(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.65_66insTTTTT"
    )

    assert result[0] == "CCCCTCCCCGCCCCACCCCGCCCCTCCCCGTTTTTCCCCACCCCG"

    assert_global_variables_zero()


def test_multi_insertion_with_left_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.20_21insCCAAA"
    )

    assert result[0] == "CCCCGCCCCACCCCGCCCCTCCAAACCCCGCCCCACCCCGCCCCTCCCCGCCCCA"

    assert_global_variables_zero()


def test_single_delins(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.38delinsAAA"
    )

    assert result[0] == "CCACCCCGCCCCTCCCCGCCCCACCCCGCCAAACTCCCCGCCCCACCCCGCCCCTCCCCGCCC"

    assert_global_variables_zero()


def test_multi_delins(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.38_40delinsAAA"
    )

    assert result[0] == "CCACCCCGCCCCTCCCCGCCCCACCCCGCCAAACCCCGCCCCACCCCGCCCCTCCCCGCCCCA"

    assert_global_variables_zero()


def test_multi_delins_with_psuedo_left_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.36_37delinsAG"
    )

    assert result[0] == "CCCCACCCCGCCCCTCCCCGCCCCACCCCGAGCCTCCCCGCCCCACCCCGCCCCTCCCCGCC"

    assert_global_variables_zero()

def test_multi_delins_with_true_left_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.36_37delinsAC"
    )

    assert result[0] == "CCCCACCCCGCCCCTCCCCGCCCCACCCCGACCCTCCCCGCCCCACCCCGCCCCTCCCCGC"

    assert_global_variables_zero()


def test_multi_delins_with_true_right_repeats(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.36_37delinsCA"
    )

    assert result[0] == "CCCACCCCGCCCCTCCCCGCCCCACCCCGCACCTCCCCGCCCCACCCCGCCCCTCCCCGCC"

    assert_global_variables_zero()

def test_single_dup(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.35dup"
    )

    assert result[0] == "CCCCACCCCGCCCCTCCCCGCCCCACCCCGGCCCCTCCCCGCCCCACCCCGCCCCTCCCC"

    assert_global_variables_zero()

def test_multi_dup(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.35_37dup"
    )

    assert result[0] == "CCACCCCGCCCCTCCCCGCCCCACCCCGCCGCCCCTCCCCGCCCCACCCCGCCCCTCC"

    assert_global_variables_zero()

def test_inversion_with_overlaps(long_sequence):
    result = gget.mutate(
        sequences=long_sequence,
        mutations="c.35_38inv"
    )

    assert result[0] == "CCCCACCCCGCCCCTCCCCGCCCCACCCCGGGCCTCCCCGCCCCACCCCGCCCCTCCCCGCC"

    assert_global_variables_zero()




def test_list_of_mutations(long_sequence):
    mutation_list = ["c.35G>A", "c.65G>A", "c.35del", "c.4_5insT"]
    sequence_list = [long_sequence for _ in range(len(mutation_list))]
    
    result = gget.mutate(
        sequences=sequence_list,
        mutations=mutation_list
    )

    assert result == ["GCCCCACCCCGCCCCTCCCCGCCCCACCCCACCCCTCCCCGCCCCACCCCGCCCCTCCCCG", "GCCCCTCCCCGCCCCACCCCGCCCCTCCCCACCCCACCCCG", "GCCCCACCCCGCCCCTCCCCGCCCCACCCCCCCCTCCCCGCCCCACCCCGCCCCTCCCCG", "CCCCTGCCCCACCCCGCCCCTCCCCGCCCCACCCC"]

    assert_global_variables_zero()


# fixme special
def test_csv_of_mutations(create_temp_files):
    mutation_temp_csv_file, sequence_temp_fasta_path = create_temp_files

    result = gget.mutate(
        sequences=sequence_temp_fasta_path,
        mutations=mutation_temp_csv_file
    )

    assert result == ["GCCCCACCCCGCCCCTCCCCGCCCCACCCCACCCCTCCCCGCCCCACCCCGCCCCTCCCCG", "GCCCCTCCCCGCCCCACCCCGCCCCTCCCCACCCCACCCCG", "GCCCCACCCCGCCCCTCCCCGCCCCACCCCCCCCTCCCCGCCCCACCCCGCCCCTCCCCG", "CCCCTGCCCCACCCCGCCCCTCCCCGCCCCACCCC"]

    assert_global_variables_zero()



def test_intron_mutation_plus(long_sequence):
    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.20+3T>A")
    
    assert_global_variables_zero(number_intronic_position_mutations=1)

def test_intron_mutation_minus(long_sequence):
    

    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.20-3T>A")

    assert_global_variables_zero(number_intronic_position_mutations=1)


def test_posttranslational_mutation(long_sequence):
    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.20*5T>A")

    assert_global_variables_zero(number_posttranslational_region_mutations=1)


def test_uncertain_mutation(long_sequence):
    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.?")

    assert_global_variables_zero(number_uncertain_mutations=1)


def test_ambiguous_mutation(long_sequence):
    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.(20_28)del")

    assert_global_variables_zero(number_ambiguous_position_mutations=1)


def test_index_error(long_sequence):
    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.99999999C>A")

    assert_global_variables_zero(number_index_errors=1)


# fixme special
def test_mismatch_error(long_sequence):
    gget.gget_mutate.mutate(
        sequences=long_sequence,
        mutations="c.2G>A")
    
    assert gget.gget_mutate.cosmic_incorrect_wt_base == 1

    assert_global_variables_zero()
