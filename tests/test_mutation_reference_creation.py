import pytest
import os
import src.cosmic_cancer_mutation_fasta_creation
import sys
from pdb import set_trace as st

# # uncomment these if I move global variables back into the main function in src.cosmic_cancer_mutation_fasta_creation
# intronic_mutations = 0
# posttranslational_region_mutations = 0
# unknown_mutations = 0
# uncertain_mutations = 0
# ambiguous_position_mutations = 0
# cosmic_incorrect_wt_base = 0

# List of argument combinations to test
argument_combinations = [
    ('cancer', '37'),
    ('cell_line', '38'),
    ('census', '38'),
    ('resistance', '38'),
    ('screen', '38'),
]

@pytest.mark.parametrize("mutation_class_name, grch_number", argument_combinations)
def test_mutation_reference_creation_main(mutation_class_name, grch_number, monkeypatch):
    # Set the TEST_MODE environment variable
    monkeypatch.setenv('TEST_MODE', '1')
    
    # Mock command-line arguments
    monkeypatch.setattr(sys, 'argv', ['my_script.py', '-m', mutation_class_name, '-g', grch_number])
    
    src.cosmic_cancer_mutation_fasta_creation.main()

    # access one of the global variables with referencing eg src.cosmic_cancer_mutation_fasta_creation.intronic_mutations

    # st()
    
# run with pytest tests/test_mutation_reference_creation.py -v -s