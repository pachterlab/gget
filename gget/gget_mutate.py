import os
import pandas as pd
from Bio import SeqIO
import re
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

# Define global variables to count occurences of weird mutations
intronic_mutations = 0
posttranslational_region_mutations = 0
unknown_mutations = 0
uncertain_mutations = 0
ambiguous_position_mutations = 0
cosmic_incorrect_wt_base = 0

mutation_pattern = r"c\.([0-9_\-\+\*]+)([a-zA-Z>]+)"  # more complex: r'c\.([0-9_\-\+\*\(\)\?]+)([a-zA-Z>\(\)0-9]+)'


def raise_pytest_error():
    if os.getenv("TEST_MODE"):
        raise ValueError()


def find_sequence(fasta_path, accession_number):
    for record in SeqIO.parse(fasta_path, "fasta"):
        if accession_number in record.description:
            return record.seq
    return "No sequence found"  # Return None if no match is found


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
            logging.debug(
                f"mutation {mutation} matches re but is not a known mutation type"
            )
            # raise_pytest_error()
            return "unknown"
    else:
        logging.debug(f"mutation {mutation} does not match re")
        # raise_pytest_error()
        return "unknown"


def load_sequences(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        key = record.description.split()[1]
        sequences[key] = str(record.seq)
    return sequences


def create_mutant_sequence(row, mutation_function, kmer_flanking_length):
    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations

    if "?" in row["Mutation CDS"]:
        logging.debug(
            f"Uncertain mutation found in {row['mut_ID']} - mutation is ignored"
        )
        uncertain_mutations += 1
        return ""

    if "(" in row["Mutation CDS"]:
        logging.debug(
            f"Ambiguous mutational position found in {row['mut_ID']} - mutation is ignored"
        )
        ambiguous_position_mutations += 1
        return ""

    match = re.match(mutation_pattern, row["Mutation CDS"])
    if match:
        nucleotide_position = match.group(1)  # The number sequence
        letters = match.group(2)  # The letter sequence

        try:
            if "-" in nucleotide_position or "+" in nucleotide_position:
                logging.debug(
                    f"Intronic nucleotide position found in {row['mut_ID']} - {nucleotide_position}{letters} - mutation is ignored"
                )
                intronic_mutations += 1
                return ""

            elif "*" in nucleotide_position:
                logging.debug(
                    f"Posttranslational region nucleotide position found in {row['mut_ID']} - {nucleotide_position}{letters} - mutation is ignored"
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
            logging.debug(f"Error with {row['mut_ID']} - row['Mutation CDS']")
            unknown_mutations += 1
            raise_pytest_error()
            return ""

        mutant_sequence, adjusted_end_position = mutation_function(
            row,
            letters,
            starting_nucleotide_position_index_0,
            ending_nucleotide_position_index_0,
        )

        if not mutant_sequence:
            return ""

        kmer_start = max(0, starting_nucleotide_position_index_0 - kmer_flanking_length)
        kmer_end = min(
            len(mutant_sequence), adjusted_end_position + kmer_flanking_length + 1
        )
        return str(mutant_sequence[kmer_start:kmer_end])
    else:
        logging.debug(f"Error with {row['mut_ID']} - {row['Mutation CDS']}")
        unknown_mutations += 1
        raise_pytest_error()
        return ""


def substitution_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):  # yields 61-mer
    # assert letters[0] == row['full_sequence'][starting_nucleotide_position_index_0], f"Transcript has {row['full_sequence'][starting_nucleotide_position_index_0]} at position {starting_nucleotide_position_index_0} but mutation is {letters[-1]} at position {starting_nucleotide_position_index_0} in {row['mut_ID']}"
    global cosmic_incorrect_wt_base
    if letters[0] != row["full_sequence"][starting_nucleotide_position_index_0]:
        logging.debug(
            f"Transcript has {row['full_sequence'][starting_nucleotide_position_index_0]} at position {starting_nucleotide_position_index_0} but mutation is {letters[-1]} at position {starting_nucleotide_position_index_0} in {row['mut_ID']}"
        )
        cosmic_incorrect_wt_base += 1
    mutant_sequence = (
        row["full_sequence"][:starting_nucleotide_position_index_0]
        + letters[-1]
        + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
    )
    adjusted_end_position = ending_nucleotide_position_index_0
    return mutant_sequence, adjusted_end_position


def deletion_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):  # yields 60-mer
    mutant_sequence = (
        row["full_sequence"][:starting_nucleotide_position_index_0]
        + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
    )
    adjusted_end_position = starting_nucleotide_position_index_0 - 1
    return mutant_sequence, adjusted_end_position


def delins_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):  # yields 60+(length of insertion)-mer
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


def insertion_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):  # yields 61+(length of insertion)-mer  - k before, length of insertion, k-1 after (k before and k-1 after is a little unconventional, but if I want to make it k-1 before then I have to return an adjusted start as well as an adjusted end)
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


def duplication_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):
    insertion_string = row["full_sequence"][
        starting_nucleotide_position_index_0 : ending_nucleotide_position_index_0 + 1
    ]
    mutant_sequence = (
        row["full_sequence"][:starting_nucleotide_position_index_0]
        + insertion_string * 2
        + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
    )
    adjusted_end_position = ending_nucleotide_position_index_0 + len(insertion_string)
    return mutant_sequence, adjusted_end_position


def inversion_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):
    insertion_string = row["full_sequence"][
        starting_nucleotide_position_index_0 : ending_nucleotide_position_index_0 + 1
    ]
    mutant_sequence = (
        row["full_sequence"][:starting_nucleotide_position_index_0]
        + insertion_string[::-1]
        + row["full_sequence"][ending_nucleotide_position_index_0 + 1 :]
    )
    adjusted_end_position = ending_nucleotide_position_index_0
    return mutant_sequence, adjusted_end_position


def unknown_mutation(
    row,
    letters,
    starting_nucleotide_position_index_0,
    ending_nucleotide_position_index_0,
):
    return "", ""


def mutate(
    mutation_df,
    input_fasta,
    k=31,
    mut_column="mutation",
    mut_id_column="mut_ID",
    seq_id_column="seq_ID",
    output_fasta="gget_mutate_out.fa",
):
    """
    Takes in a list of nucleotide mutations and sequences and returns mutated versions of the input sequences
    according to the provided mutations.

    Args:
    - mutation_df   Path to csv file (str) or data frame (DataFrame object) in the following format:

                    | mutation             | mut_ID | seq_ID |
                    | c.1252C>T            | mut1   | seq1   |
                    | c.2239_2253inv       | mut2   | seq2   |
                    | c.2239_2253inv       | mut2   | seq4   |
                    | c.2239_2253delinsAAT | mut3   | seq4   |
                    | ...                  | ...    | ...    |

                    'mutation' = Column containing the mutations to be performed written in standard mutation annotation (see below)
                    'mut_ID' = Column containing the IDs of each mutation
                    'seq_ID' = Column containing the IDs of the sequences to be mutated (must correspond to the string following the > character in the input_fasta;
                    do not include spaces)

    - input_fasta   (str) Path to fasta file containing the sequences to be mutated.
    - k             (int) Length of sequences flanking the mutated base/amino acid. Default: 31.
                    If k > total length of the sequence, the entire sequence will be kept.
    - mut_column    (str) Name of the column containing the mutations to be performed in mutation_df. Default: 'mutation'.
    - mut_id_column (str) Name of the column containing the IDs of each mutation in mutation_df. Default: 'mut_ID'.
    - seq_id_column (str) Name of the column containing the IDs of the sequences to be mutated in mutation_df. Default: 'seq_ID'.
    - output_fasta  (str) Path to output fasta file containing the mutated sequences. Default: 'gget_mutate_out.fa'.

    For more on the standard mutation annotation, also see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1867422/.

    Saves mutated sequences in a fasta file as specified in the output_fasta argument.
    """

    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations, cosmic_incorrect_wt_base

    intronic_mutations = 0
    posttranslational_region_mutations = 0
    unknown_mutations = 0
    uncertain_mutations = 0
    ambiguous_position_mutations = 0
    cosmic_incorrect_wt_base = 0

    # Read in mutation_df if passed as filepath instead of df
    if isinstance(mutation_df, str):
        mutation_df = pd.read_csv(mutation_df)

    logging.info("Extracting mutation types...")
    # Get all mutation types
    mutation_df["mutation_type"] = mutation_df[mut_column].progress_apply(
        extract_mutation_type
    )

    # Load input sequences and link sequences to their mutations using the sequence identifier
    sequences = load_sequences(input_fasta)
    mutation_df["full_sequence"] = mutation_df[seq_id_column].map(sequences)

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
        df_mutation_type = mutation_df[
            mutation_df["mutation_type"] == mutation_type
        ].copy()
        df_mutation_type["mutant_sequence_kmer"] = ""
        df_mutation_type["header"] = (
            ">"
            + df_mutation_type[seq_id_column]
            + "_"
            + df_mutation_type[mut_id_column]
        )
        mutation_dict[mutation_type] = df_mutation_type

    # Create mutated sequences
    logging.info("Mutating sequences...")
    if not mutation_dict["substitution"].empty:
        mutation_dict["substitution"]["mutant_sequence_kmer"] = mutation_dict[
            "substitution"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                substitution_mutation,
                k,
            ),
            axis=1,
        )
    if not mutation_dict["deletion"].empty:
        mutation_dict["deletion"]["mutant_sequence_kmer"] = mutation_dict[
            "deletion"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                deletion_mutation,
                k,
            ),
            axis=1,
        )
    if not mutation_dict["delins"].empty:
        mutation_dict["delins"]["mutant_sequence_kmer"] = mutation_dict[
            "delins"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                delins_mutation,
                k,
            ),
            axis=1,
        )
    if not mutation_dict["insertion"].empty:
        mutation_dict["insertion"]["mutant_sequence_kmer"] = mutation_dict[
            "insertion"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                insertion_mutation,
                k,
            ),
            axis=1,
        )
    if not mutation_dict["duplication"].empty:
        mutation_dict["duplication"]["mutant_sequence_kmer"] = mutation_dict[
            "duplication"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                duplication_mutation,
                k,
            ),
            axis=1,
        )
    if not mutation_dict["inversion"].empty:
        mutation_dict["inversion"]["mutant_sequence_kmer"] = mutation_dict[
            "inversion"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                inversion_mutation,
                k,
            ),
            axis=1,
        )
    if not mutation_dict["unknown"].empty:
        mutation_dict["unknown"]["mutant_sequence_kmer"] = mutation_dict[
            "unknown"
        ].progress_apply(
            create_mutant_sequence,
            args=(
                unknown_mutation,
                k,
            ),
            axis=1,
        )

    # Report status of mutations back to user
    total_mutations = df.shape[0]
    good_mutations = (
        total_mutations
        - intronic_mutations
        - posttranslational_region_mutations
        - unknown_mutations
        - uncertain_mutations
        - ambiguous_position_mutations
        - cosmic_incorrect_wt_base
    )

    logging.warning(
        f"""
    {good_mutations} mutations correctly recorded ({good_mutations/total_mutations*100:.2f}%)
    {intronic_mutations} intronic mutations found ({intronic_mutations/total_mutations*100:.2f}%)
    {posttranslational_region_mutations} posttranslational region mutations found ({posttranslational_region_mutations/total_mutations*100:.2f}%)
    {unknown_mutations} unknown mutations found ({unknown_mutations/total_mutations*100:.2f}%)
    {uncertain_mutations} mutations with uncertain mutation found ({uncertain_mutations/total_mutations*100:.2f}%)
    {ambiguous_position_mutations} mutations with ambiguous position found ({ambiguous_position_mutations/total_mutations*100:.2f}%)
    {cosmic_incorrect_wt_base} mutations with incorrect wildtype base found ({cosmic_incorrect_wt_base/total_mutations*100:.2f}%)
    """
    )

    # Save mutated sequences in new fasta file
    logging.info("Creating FASTA file containing mutated sequences...")
    with open(output_fasta, "w") as fasta_file:
        for mutation in mutation_dict:
            if not mutation_dict[mutation].empty:
                df = mutation_dict[mutation]
                df["fasta_format"] = (
                    df["header"] + "\n" + df["mutant_sequence_kmer"] + "\n"
                )
                # df = df.dropna(subset=['GENOMIC_MUTATION_ID', 'fasta_format'])
                fasta_file.write("\n".join(df["fasta_format"]))

    with open(output_fasta, "r") as fasta_file:
        lines = [line for line in fasta_file if line.strip()]

    with open(output_fasta, "w") as fasta_file:
        fasta_file.writelines(lines)

    logging.info(f"FASTA file containing mutated sequences created at {output_fasta}.")
