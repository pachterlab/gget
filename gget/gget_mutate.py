import pandas as pd
import re
from tqdm import tqdm

tqdm.pandas()

from .utils import read_fasta, set_up_logger

logger = set_up_logger()

# Define global variables to count occurences of weird mutations
intronic_mutations = 0
posttranslational_region_mutations = 0
unknown_mutations = 0
uncertain_mutations = 0
ambiguous_position_mutations = 0
cosmic_incorrect_wt_base = 0
mut_idx_outside_seq = 0

mutation_pattern = r"c\.([0-9_\-\+\*]+)([a-zA-Z>]+)"  # more complex: r'c\.([0-9_\-\+\*\(\)\?]+)([a-zA-Z>\(\)0-9]+)'


def extract_mutation_type(mutation):
    match = re.match(mutation_pattern, mutation)

    if match:
        letters = match.group(2)  # The letter sequence
        mutation_type = re.findall(r"[a-z>-]+", letters)
        if mutation_type[0] == ">":  # substitution - eg c.1252C>T
            return "substitution"

        elif mutation_type[0] == "del":  # deletion - eg c.2126_2128del OR c.1627del
            return "deletion"

        elif (
            mutation_type[0] == "delins"
        ):  # insertion-deletion - eg c.2239_2253delinsAAT OR c.646delins
            return "delins"

        elif (
            mutation_type[0] == "ins"
        ):  # insertion - eg c.2239_2240insAAT   (interval is always 1)
            return "insertion"

        elif mutation_type[0] == "dup":  # duplication - eg c.2239_2253dup OR c.646dup
            return "duplication"

        elif (
            mutation_type[0] == "inv"
        ):  # inversion - eg c.2239_2253inv    (always has an interval)
            return "inversion"

        else:
            logger.debug(
                f"mutation {mutation} matches re but is not a known mutation type"
            )
            # raise_pytest_error()
            return "unknown"
    else:
        logger.debug(f"mutation {mutation} does not match re")
        # raise_pytest_error()
        return "unknown"


def create_mutant_sequence(
    row,
    mutation_function,
    kmer_flanking_length=30,
    mut_column="mutation",
    seq_id_column="seq_ID",
):
    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations

    if "?" in row[mut_column]:
        logger.debug(
            f"Uncertain mutation found in mutation {row[mut_column]} - mutation is ignored"
        )
        uncertain_mutations += 1
        return ""

    if "(" in row[mut_column]:
        logger.debug(
            f"Ambiguous mutational position found in mutation {row[mut_column]} - mutation is ignored"
        )
        ambiguous_position_mutations += 1
        return ""

    match = re.match(mutation_pattern, row[mut_column])
    if match:
        nucleotide_position = match.group(1)  # The number sequence
        letters = match.group(2)  # The letter sequence

        try:
            if "-" in nucleotide_position or "+" in nucleotide_position:
                logger.debug(
                    f"Intronic nucleotide position found in mutation {row[mut_column]} - {nucleotide_position}{letters} - mutation is ignored"
                )
                intronic_mutations += 1
                return ""

            elif "*" in nucleotide_position:
                logger.debug(
                    f"Posttranslational region nucleotide position found in mutation {row[mut_column]} - {nucleotide_position}{letters} - mutation is ignored"
                )
                posttranslational_region_mutations += 1
                return ""

            else:
                if "_" in nucleotide_position:
                    starting_nucleotide_position_index_0 = (
                        int(nucleotide_position.split("_")[0]) - 1
                    )
                    ending_nucleotide_position_index_0 = (
                        int(nucleotide_position.split("_")[1]) - 1
                    )
                else:
                    starting_nucleotide_position_index_0 = int(nucleotide_position) - 1
                    ending_nucleotide_position_index_0 = (
                        starting_nucleotide_position_index_0
                    )
        except:
            logger.debug(f"Error with mutation {row[mut_column]}]")
            unknown_mutations += 1
            # raise_pytest_error()
            return ""

        mutant_sequence, adjusted_end_position = mutation_function(
            row,
            letters,
            starting_nucleotide_position_index_0,
            ending_nucleotide_position_index_0,
            mut_column,
            seq_id_column,
        )

        if not mutant_sequence:
            return ""

        kmer_start = max(0, starting_nucleotide_position_index_0 - kmer_flanking_length)
        kmer_end = min(
            len(mutant_sequence), adjusted_end_position + kmer_flanking_length + 1
        )
        return str(mutant_sequence[kmer_start:kmer_end])
    else:
        logger.debug(f"Error with mutation {row[mut_column]}")
        unknown_mutations += 1
        # raise_pytest_error()
        return ""


def substitution_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    # assert letters[0] == row['full_sequence'][starting_nucleotide_position_index_0], f"Transcript has {row['full_sequence'][starting_nucleotide_position_index_0]} at position {starting_nucleotide_position_index_0} but mutation is {letters[-1]} at position {starting_nucleotide_position_index_0} in {row[mut_column]}"
    global cosmic_incorrect_wt_base
    global mut_idx_outside_seq
    try:
        if letters[0] != row["full_sequence"][starting_nucleotide_position_index_0]:
            logger.debug(
                f"Sequence {row[seq_id_column]} has nucleotide '{row['full_sequence'][starting_nucleotide_position_index_0]}' at position {starting_nucleotide_position_index_0}, but mutation {row[mut_column]} expected '{letters[0]}'."
            )
            cosmic_incorrect_wt_base += 1
        mutant_sequence = (
            row["full_sequence"][:starting_nucleotide_position_index_0]
            + letters[-1]
            + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
        )
        adjusted_end_position = ending_nucleotide_position_index_0

        return mutant_sequence, adjusted_end_position

    # Mutation index is outside sequence length
    except IndexError:
        mut_idx_outside_seq += 1
        return "", ""


def deletion_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    global mut_idx_outside_seq
    try:
        # Accessing the sequence to trigger IndexError if out of bounds
        _ = row["full_sequence"][starting_nucleotide_position_index_0]

        mutant_sequence = (
            row["full_sequence"][:starting_nucleotide_position_index_0]
            + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
        )
        adjusted_end_position = starting_nucleotide_position_index_0 - 1
        return mutant_sequence, adjusted_end_position

    # Mutation index is outside sequence length
    except IndexError:
        mut_idx_outside_seq += 1
        return "", ""


def delins_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    global mut_idx_outside_seq
    try:
        # Accessing the sequence to trigger IndexError if out of bounds
        _ = row["full_sequence"][starting_nucleotide_position_index_0]

        # yields 60+(length of insertion)-mer
        insertion_string = "".join(re.findall(r"[A-Z]+", letters))
        mutant_sequence = (
            row["full_sequence"][:starting_nucleotide_position_index_0]
            + insertion_string
            + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
        )
        adjusted_end_position = (
            starting_nucleotide_position_index_0 + len(insertion_string) - 1
        )
        return mutant_sequence, adjusted_end_position

    # Mutation index is outside sequence length
    except IndexError:
        mut_idx_outside_seq += 1
        return "", ""


def insertion_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    global mut_idx_outside_seq
    try:
        # Accessing the sequence to trigger IndexError if out of bounds
        _ = row["full_sequence"][starting_nucleotide_position_index_0]

        # yields 61+(length of insertion)-mer  - k before, length of insertion, k-1 after (k before and k-1 after is a little unconventional, but if I want to make it k-1 before then I have to return an adjusted start as well as an adjusted end)
        insertion_string = "".join(re.findall(r"[A-Z]+", letters))
        mutant_sequence = (
            row["full_sequence"][: starting_nucleotide_position_index_0 + 1]
            + insertion_string
            + row["full_sequence"][starting_nucleotide_position_index_0 + 1 :]
        )
        adjusted_end_position = (
            ending_nucleotide_position_index_0 + len(insertion_string) - 1
        )
        return mutant_sequence, adjusted_end_position

    # Mutation index is outside sequence length
    except IndexError:
        mut_idx_outside_seq += 1
        return "", ""


def duplication_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    global mut_idx_outside_seq
    try:
        # Accessing the sequence to trigger IndexError if out of bounds
        _ = row["full_sequence"][starting_nucleotide_position_index_0]
        insertion_string = row["full_sequence"][
            starting_nucleotide_position_index_0 : ending_nucleotide_position_index_0
            + 1
        ]
        mutant_sequence = (
            row["full_sequence"][:starting_nucleotide_position_index_0]
            + insertion_string * 2
            + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
        )
        adjusted_end_position = ending_nucleotide_position_index_0 + len(
            insertion_string
        )
        return mutant_sequence, adjusted_end_position

    # Mutation index is outside sequence length
    except IndexError:
        mut_idx_outside_seq += 1
        return "", ""


def inversion_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    global mut_idx_outside_seq
    try:
        # Accessing the sequence to trigger IndexError if out of bounds
        _ = row["full_sequence"][starting_nucleotide_position_index_0]
        insertion_string = row["full_sequence"][
            starting_nucleotide_position_index_0 : ending_nucleotide_position_index_0
            + 1
        ]

        # Reverse
        reverse_insertion_string = insertion_string[::-1]

        # Get complement
        complement = {
            "A": "T",
            "T": "A",
            "U": "A",
            "C": "G",
            "G": "C",
            "N": "N",
            "a": "t",
            "t": "a",
            "u": "a",
            "c": "g",
            "g": "c",
            "n": "n",
            ".": ".",  # annotation for gaps
            "-": "-",  # annotation for gaps
        }
        mutated_string = "".join(
            complement.get(nucleotide, "N") for nucleotide in reverse_insertion_string
        )

        # Create final sequence
        mutant_sequence = (
            row["full_sequence"][:starting_nucleotide_position_index_0]
            + mutated_string
            + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
        )
        adjusted_end_position = ending_nucleotide_position_index_0
        return mutant_sequence, adjusted_end_position

    # Mutation index is outside sequence length
    except IndexError:
        mut_idx_outside_seq += 1
        return "", ""


def unknown_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
    mut_column,
    seq_id_column,
):
    return "", ""


def mutate(
    sequences,
    mutations,
    k=30,
    mut_column="mutation",
    mut_id_column="mut_ID",
    seq_id_column="seq_ID",
    out=None,
    verbose=True,
):
    """
    Takes in nucleotide sequences and mutations (in standard mutation annotation - see below)
    and returns mutated versions of the input sequences according to the provided mutations.

    Args:
    - sequences     (str) Path to the fasta file containing the sequences to be mutated, e.g., 'seqs.fa'.
                    Sequence identifiers following the '>' character must correspond to the identifiers
                    in the seq_ID column of 'mutations'.
                    NOTE: Only string until first space or dot will be used as sequence identifier
                    - Version numbers of Ensembl IDs will be ignored.

                    Example:
                    >seq1 (or ENSG00000106443)
                    ACTGCGATAGACT
                    >seq2
                    AGATCGCTAG

                    Alternatively: Input sequence(s) as a string or list, e.g. 'AGCTAGCT' or ['ACTGCTAGCT', 'AGCTAGCT'].
    - mutations     Path to csv or tsv file (str) (e.g., 'mutations.csv') or data frame (DataFrame object)
                    containing information about the mutations in the following format:

                    | mutation         | mut_ID | seq_ID |
                    | c.2C>T           | mut1   | seq1   | -> Apply mutation 1 to sequence 1
                    | c.9_13inv        | mut2   | seq2   | -> Apply mutation 2 to sequence 2
                    | c.9_13inv        | mut2   | seq3   | -> Apply mutation 2 to sequence 3
                    | c.9_13delinsAAT  | mut3   | seq3   | -> Apply mutation 3 to sequence 3
                    | ...              | ...    | ...    |

                    'mutation' = Column containing the mutations to be performed written in standard mutation annotation (see below)
                    'mut_ID' = Column containing an identifier for each mutation
                    'seq_ID' = Column containing the identifiers of the sequences to be mutated (must correspond to the string following
                    the > character in the 'sequences' fasta file; do NOT include spaces or dots)

                    Alternatively: Input mutation(s) as a string or list, e.g., 'c.2C>T' or ['c.2C>T', 'c.1A>C'].
                    If a list is provided, the number of mutations must equal the number of input sequences.
    - k             (int) Length of sequences flanking the mutation. Default: 30.
                    If k > total length of the sequence, the entire sequence will be kept.
    - mut_column    (str) Name of the column containing the mutations to be performed in 'mutations'. Default: 'mutation'.
    - mut_id_column (str) Name of the column containing the IDs of each mutation in 'mutations'. Default: 'mut_ID'.
    - seq_id_column (str) Name of the column containing the IDs of the sequences to be mutated in 'mutations'. Default: 'seq_ID'.
    - out           (str) Path to output fasta file containing the mutated sequences, e.g., 'path/to/output_fasta.fa'.
                    Default: None -> returns a list of the mutated sequences to standard out.
                    The identifiers (following the '>') of the mutated sequences in the output fasta will be '>[seq_ID]_[mut_ID]'.
    - verbose       (True/False) whether to print progress information. Default: True

    For more information on the standard mutation annotation, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1867422/.

    Saves mutated sequences in fasta format (or returns a list containing the mutated sequences if out=None).
    """

    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations, cosmic_incorrect_wt_base, mut_idx_outside_seq

    # Load input sequences and their identifiers from fasta file
    if "." in sequences:
        titles, seqs = read_fasta(sequences)

    # Handle input sequences passed as a list
    elif isinstance(sequences, list):
        titles = [f"seq{i+1}" for i in range(len(sequences))]
        seqs = sequences

    # Handle a single sequence passed as a string
    elif isinstance(sequences, str):
        titles = ["seq1"]
        seqs = [sequences]

    else:
        raise ValueError(
            """
            Format of the input to the 'sequences' argument not recognized. 
            'sequences' must be one of the following:
            - Path to the fasta file containing the sequences to be mutated (e.g. 'seqs.fa')
            - A list of sequences to be mutated (e.g. ['ACTGCTAGCT', 'AGCTAGCT'])
            - A single sequence to be mutated passed as a string (e.g. 'AGCTAGCT')
            """
        )

    # Read in 'mutations' if passed as filepath to comma-separated csv
    if isinstance(mutations, str) and ".csv" in mutations:
        mutations = pd.read_csv(mutations)

    elif isinstance(mutations, str) and ".tsv" in mutations:
        mutations = pd.read_csv(mutations, sep="\t")

    # Handle mutations passed as a list
    elif isinstance(mutations, list):
        if len(mutations) > 1:
            if len(mutations) != len(seqs):
                raise ValueError(
                    "If a list is passed, the number of mutations must equal the number of input sequences."
                )

            temp = pd.DataFrame()
            temp["mutation"] = mutations
            temp["mut_ID"] = [f"mut{i+1}" for i in range(len(mutations))]
            temp["seq_ID"] = [f"seq{i+1}" for i in range(len(mutations))]
            mutations = temp
        else:
            temp = pd.DataFrame()
            temp["mutation"] = [mutations[0]] * len(seqs)
            temp["mut_ID"] = [f"mut{i+1}" for i in range(len(seqs))]
            temp["seq_ID"] = [f"seq{i+1}" for i in range(len(seqs))]
            mutations = temp

    # Handle single mutation passed as a string
    elif isinstance(mutations, str):
        # This will work for one mutation for one sequence as well as one mutation for multiple sequences
        temp = pd.DataFrame()
        temp["mutation"] = [mutations] * len(seqs)
        temp["mut_ID"] = [f"mut{i+1}" for i in range(len(seqs))]
        temp["seq_ID"] = [f"seq{i+1}" for i in range(len(seqs))]
        mutations = temp

    elif isinstance(mutations, pd.DataFrame):
        pass

    else:
        raise ValueError(
            """
            Format of the input to the 'mutations' argument not recognized. 
            'mutations' must be one of the following:
            - Path to comma-separated csv file (e.g. 'mutations.csv')
            - A pandas DataFrame object
            - A single mutation to be applied to all input sequences (e.g. 'c.2C>T')
            - A list of mutations (the number of mutations must equal the number of input sequences) (e.g. ['c.2C>T', 'c.1A>C'])
            """
        )

    # Set of possible nucleotides (- and . are gap annotations)
    nucleotides = set("ATGCUNatgcun.-")

    seq_dict = {}
    non_nuc_seqs = 0
    for title, seq in zip(titles, seqs):
        # Check that sequences are nucleotide sequences
        if not set(seq) <= nucleotides:
            non_nuc_seqs += 1

        # Keep text following the > until the first space/dot as the sequence identifier
        # Dots are removed so Ensembl version numbers are removed
        seq_dict[title.split(" ")[0].split(".")[0]] = seq

    if non_nuc_seqs > 0:
        logger.warning(
            f"""
            Non-nucleotide characters detected in {non_nuc_seqs} input sequences. gget mutate is currently only optimized for mutating nucleotide sequences.
            Specifically inversion mutations might not be performed correctly. 
            """
        )

    # Get all mutation types
    if verbose:
        tqdm.pandas(desc="Extracting mutation types")
        mutations["mutation_type"] = mutations[mut_column].progress_apply(
            extract_mutation_type
        )
    else:
        mutations["mutation_type"] = mutations[mut_column].apply(extract_mutation_type)

    # Link sequences to their mutations using the sequence identifiers
    mutations["full_sequence"] = mutations[seq_id_column].map(seq_dict)

    # Handle sequences that were not found based on their sequence IDs
    seqs_not_found = mutations[mutations["full_sequence"].isnull()]
    if 0 < len(seqs_not_found) < 20:
        logger.warning(
            f"""
            The sequences with the following {len(seqs_not_found)} sequence ID(s) were not found: {", ".join(seqs_not_found["seq_ID"].values)}  
            These sequences and their corresponding mutations will not be included in the output.  
            Ensure that the sequence IDs correspond to the string following the > character in the 'sequences' fasta file (do NOT include spaces or dots).
            """
        )
    elif len(seqs_not_found) > 0:
        logger.warning(
            f"""
            The sequences corresponding to {len(seqs_not_found)} sequence IDs were not found.  
            These sequences and their corresponding mutations will not be included in the output.  
            Ensure that the sequence IDs correspond to the string following the > character in the 'sequences' fasta file (do NOT include spaces or dots).
            """
        )

    # Drop inputs for sequences that were not found
    mutations = mutations.dropna()
    if len(mutations) < 1:
        raise ValueError(
            """
            None of the input sequences match the sequence IDs provided in 'mutations'. 
            Ensure that the sequence IDs correspond to the string following the > character in the 'sequences' fasta file (do NOT include spaces or dots).
            """
        )

    # Split data frame by mutation type
    mutation_types = [
        "substitution",
        "deletion",
        "delins",
        "insertion",
        "duplication",
        "inversion",
        "unknown",
    ]
    mutation_dict = {mutation_type: [] for mutation_type in mutation_types}
    for mutation_type in mutation_types:
        df_mutation_type = mutations[mutations["mutation_type"] == mutation_type].copy()
        df_mutation_type["mutant_sequence_kmer"] = ""
        # Define header for mutated sequences in output fasta
        df_mutation_type["header"] = (
            ">"
            + df_mutation_type[seq_id_column]
            + "_"
            + df_mutation_type[mut_id_column]
        )
        mutation_dict[mutation_type] = df_mutation_type

    # Create mutated sequences
    if verbose:
        logger.info("Mutating sequences...")
    if not mutation_dict["substitution"].empty:
        if verbose:
            tqdm.pandas(desc="Performing substitutions")
            mutation_dict["substitution"]["mutant_sequence_kmer"] = mutation_dict[
                "substitution"
            ].progress_apply(
                create_mutant_sequence,
                args=(substitution_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["substitution"]["mutant_sequence_kmer"] = mutation_dict[
                "substitution"
            ].apply(
                create_mutant_sequence,
                args=(substitution_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
    if not mutation_dict["deletion"].empty:
        if verbose:
            tqdm.pandas(desc="Performing deletions")
            mutation_dict["deletion"]["mutant_sequence_kmer"] = mutation_dict[
                "deletion"
            ].progress_apply(
                create_mutant_sequence,
                args=(deletion_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["deletion"]["mutant_sequence_kmer"] = mutation_dict[
                "deletion"
            ].apply(
                create_mutant_sequence,
                args=(deletion_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
    if not mutation_dict["delins"].empty:
        if verbose:
            tqdm.pandas(desc="Performing delins")
            mutation_dict["delins"]["mutant_sequence_kmer"] = mutation_dict[
                "delins"
            ].progress_apply(
                create_mutant_sequence,
                args=(delins_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["delins"]["mutant_sequence_kmer"] = mutation_dict[
                "delins"
            ].apply(
                create_mutant_sequence,
                args=(delins_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
    if not mutation_dict["insertion"].empty:
        if verbose:
            tqdm.pandas(desc="Performing insertions")
            mutation_dict["insertion"]["mutant_sequence_kmer"] = mutation_dict[
                "insertion"
            ].progress_apply(
                create_mutant_sequence,
                args=(insertion_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["insertion"]["mutant_sequence_kmer"] = mutation_dict[
                "insertion"
            ].apply(
                create_mutant_sequence,
                args=(insertion_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
    if not mutation_dict["duplication"].empty:
        if verbose:
            tqdm.pandas(desc="Performing duplications")
            mutation_dict["duplication"]["mutant_sequence_kmer"] = mutation_dict[
                "duplication"
            ].progress_apply(
                create_mutant_sequence,
                args=(duplication_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["duplication"]["mutant_sequence_kmer"] = mutation_dict[
                "duplication"
            ].apply(
                create_mutant_sequence,
                args=(duplication_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
    if not mutation_dict["inversion"].empty:
        if verbose:
            tqdm.pandas(desc="Performing inversions")
            mutation_dict["inversion"]["mutant_sequence_kmer"] = mutation_dict[
                "inversion"
            ].progress_apply(
                create_mutant_sequence,
                args=(inversion_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["inversion"]["mutant_sequence_kmer"] = mutation_dict[
                "inversion"
            ].apply(
                create_mutant_sequence,
                args=(inversion_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
    if not mutation_dict["unknown"].empty:
        if verbose:
            tqdm.pandas(desc="Unknown mutations")
            mutation_dict["unknown"]["mutant_sequence_kmer"] = mutation_dict[
                "unknown"
            ].progress_apply(
                create_mutant_sequence,
                args=(unknown_mutation, k, mut_column, seq_id_column),
                axis=1,
            )
        else:
            mutation_dict["unknown"]["mutant_sequence_kmer"] = mutation_dict[
                "unknown"
            ].apply(
                create_mutant_sequence,
                args=(unknown_mutation, k, mut_column, seq_id_column),
                axis=1,
            )

    # Report status of mutations back to user
    total_mutations = mutations.shape[0]
    good_mutations = (
        total_mutations
        - intronic_mutations
        - posttranslational_region_mutations
        - unknown_mutations
        - uncertain_mutations
        - ambiguous_position_mutations
        - cosmic_incorrect_wt_base
        - mut_idx_outside_seq
    )

    if good_mutations != total_mutations:
        logger.warning(
            f"""
            {good_mutations} mutations correctly recorded ({good_mutations/total_mutations*100:.2f}%)
            {intronic_mutations} intronic mutations found ({intronic_mutations/total_mutations*100:.2f}%)
            {posttranslational_region_mutations} posttranslational region mutations found ({posttranslational_region_mutations/total_mutations*100:.2f}%)
            {unknown_mutations} unknown mutations found ({unknown_mutations/total_mutations*100:.2f}%)
            {uncertain_mutations} mutations with uncertain mutation found ({uncertain_mutations/total_mutations*100:.2f}%)
            {ambiguous_position_mutations} mutations with ambiguous position found ({ambiguous_position_mutations/total_mutations*100:.2f}%)
            {cosmic_incorrect_wt_base} mutations with incorrect wildtype base found ({cosmic_incorrect_wt_base/total_mutations*100:.2f}%)
            {mut_idx_outside_seq} mutations with indeces outside of the sequence length found ({mut_idx_outside_seq/total_mutations*100:.2f}%)
            """
        )
    else:
        if verbose:
            logger.info(
                f"""
                {good_mutations} mutations correctly recorded ({good_mutations/total_mutations*100:.2f}%)
                {intronic_mutations} intronic mutations found ({intronic_mutations/total_mutations*100:.2f}%)
                {posttranslational_region_mutations} posttranslational region mutations found ({posttranslational_region_mutations/total_mutations*100:.2f}%)
                {unknown_mutations} unknown mutations found ({unknown_mutations/total_mutations*100:.2f}%)
                {uncertain_mutations} mutations with uncertain mutation found ({uncertain_mutations/total_mutations*100:.2f}%)
                {ambiguous_position_mutations} mutations with ambiguous position found ({ambiguous_position_mutations/total_mutations*100:.2f}%)
                {cosmic_incorrect_wt_base} mutations with incorrect wildtype base found ({cosmic_incorrect_wt_base/total_mutations*100:.2f}%)
                {mut_idx_outside_seq} mutations with indeces outside of the sequence length found ({mut_idx_outside_seq/total_mutations*100:.2f}%)
                """
            )

    if out:
        # Save mutated sequences in new fasta file
        with open(out, "w") as fasta_file:
            for mutation in mutation_dict:
                if not mutation_dict[mutation].empty:
                    df = mutation_dict[mutation]
                    df["fasta_format"] = (
                        df["header"] + "\n" + df["mutant_sequence_kmer"] + "\n"
                    )

                    # Do not include an empty string as a mutated sequence
                    # (these are introduced when unknown mutations are encountered)
                    df = df[df["mutant_sequence_kmer"] != ""]

                    fasta_file.write("".join(df["fasta_format"].values))

        # with open(out, "r") as fasta_file:
        #     lines = [line for line in fasta_file if line.strip()]

        # with open(out, "w") as fasta_file:
        #     fasta_file.writelines(lines)

        if verbose:
            logger.info(f"FASTA file containing mutated sequences created at {out}.")

    # When out=None, return list of mutated seqs
    else:
        all_mut_seqs = []
        for mutation in mutation_dict:
            df = mutation_dict[mutation]
            all_mut_seqs.extend(df["mutant_sequence_kmer"].values)

        # Remove empty strings from final list of mutated sequences
        # (these are introduced when unknown mutations are encountered)
        while "" in all_mut_seqs:
            all_mut_seqs.remove("")

        if len(all_mut_seqs) > 0:
            return all_mut_seqs
