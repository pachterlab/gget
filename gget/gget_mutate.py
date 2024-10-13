import pandas as pd
import re
from tqdm import tqdm
import numpy as np
import os
from typing import Union, List, Optional

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

mutation_pattern = r"(?:c|g)\.([0-9_\-\+\*]+)([a-zA-Z>]+)"  # more complex: r'c\.([0-9_\-\+\*\(\)\?]+)([a-zA-Z>\(\)0-9]+)'

def reverse_complement(seq):
    if pd.isna(seq):  # Check if the sequence is NaN
        return np.nan
    complement = str.maketrans("ATCGNatcgn.*", "TAGCNtagcn.*")
    return seq.translate(complement)[::-1]

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


codon_to_amino_acid = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def convert_chromosome_value_to_int_when_possible(val):
    try:
        # Try to convert the value to a float, then to an int, and finally to a string
        return str(int(float(val)))
    except ValueError:
        # If conversion fails, keep the value as it is
        return val


def merge_gtf_transcript_locations_into_cosmic_csv(mutations, gtf_path, gtf_transcript_id_column):
    gtf_df = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, names=[
    'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    if 'strand' in mutations.columns:
        mutations.rename(columns={'strand': 'strand_original'}, inplace=True)

    gtf_df = gtf_df[gtf_df['feature'] == 'transcript']

    gtf_df['transcript_id'] = gtf_df['attribute'].str.extract('transcript_id "([^"]+)"')

    assert len(gtf_df['transcript_id']) == len(set(gtf_df['transcript_id'])), "Duplicate transcript_id values found!"

    # Filter out rows where transcript_id is NaN
    gtf_df = gtf_df.dropna(subset=['transcript_id'])

    gtf_df = gtf_df[['transcript_id', 'start', 'end', 'strand']].rename(
        columns={'transcript_id': gtf_transcript_id_column, 'start': 'start_transcript_position', 'end': 'end_transcript_position'})
    
    merged_df = pd.merge(mutations, gtf_df, on=gtf_transcript_id_column, how='left')

    # Fill NaN values
    merged_df['start_transcript_position'] = merged_df['start_transcript_position'].fillna(0)
    merged_df['end_transcript_position'] = merged_df['end_transcript_position'].fillna(9999999)
    merged_df['strand'] = merged_df['strand'].fillna('.')

    return merged_df

def get_sequence_length(seq_id, seq_dict):
    return len(seq_dict.get(seq_id, ""))

def get_nucleotide_at_position(seq_id, pos, seq_dict):
    full_seq = seq_dict.get(seq_id, "")
    if pos < len(full_seq):
        return full_seq[pos]
    return None

def translate_sequence(sequence, start, end):
    amino_acid_sequence = ""
    for i in range(start, end, 3):
        codon = sequence[i : i + 3].upper()
        amino_acid = codon_to_amino_acid.get(
            codon, "X"
        )  # Use 'X' for unknown or incomplete codons
        amino_acid_sequence += amino_acid

    return amino_acid_sequence


# def remove_all_but_first_gt(line):
#     return line[:1] + line[1:].replace(">", "")

def remove_gt_after_semicolon(line):
    parts = line.split(';')
    # Remove '>' from the beginning of each part except the first part
    parts = [parts[0]] + [part.lstrip('>') for part in parts[1:]]
    return ';'.join(parts)


def wt_fragment_and_mutant_fragment_share_kmer(mutated_fragment: str, wildtype_fragment: str, k: int) -> bool:
    if len(mutated_fragment) <= k:
        if mutated_fragment in wildtype_fragment:
            return True
        else:
            return False
    else:
        for mutant_position in range(len(mutated_fragment) - k):
            mutant_kmer = mutated_fragment[mutant_position:mutant_position + k]
            if mutant_kmer in wildtype_fragment:
                # wt_position = wildtype_fragment.find(mutant_kmer)
                return True
        return False


def add_mutation_type(mutations, mut_column):
    mutations["mutation_type_id"] = mutations[mut_column].str.extract(mutation_pattern)[
        1
    ]

    # Define conditions and choices for the mutation types
    conditions = [
        mutations["mutation_type_id"].str.contains(">", na=False),
        mutations["mutation_type_id"].str.contains("delins", na=False),
        mutations["mutation_type_id"].str.contains("del", na=False)
        & ~mutations["mutation_type_id"].str.contains("delins", na=False),
        mutations["mutation_type_id"].str.contains("ins", na=False)
        & ~mutations["mutation_type_id"].str.contains("delins", na=False),
        mutations["mutation_type_id"].str.contains("dup", na=False),
        mutations["mutation_type_id"].str.contains("inv", na=False),
    ]

    choices = [
        "substitution",
        "delins",
        "deletion",
        "insertion",
        "duplication",
        "inversion",
    ]

    # Assign the mutation types
    mutations["mutation_type"] = np.select(conditions, choices, default="unknown")

    # Drop the temporary mutation_type_id column
    mutations.drop(columns=["mutation_type_id"], inplace=True)

    return mutations


def extract_sequence(row, seq_dict, seq_id_column = "seq_ID"):
    if pd.isna(row["start_mutation_position"]) or pd.isna(row["end_mutation_position"]):
        return None    
    seq = seq_dict[row[seq_id_column]][
        int(row["start_mutation_position"]) : int(row["end_mutation_position"]) + 1
    ]
    return seq


def common_prefix_length(s1, s2):
    min_len = min(len(s1), len(s2))
    for i in range(min_len):
        if s1[i] != s2[i]:
            return i
    return min_len


# Function to find the length of the common suffix with the prefix
def common_suffix_length(s1, s2):
    min_len = min(len(s1), len(s2))
    for i in range(min_len):
        if s1[-(i + 1)] != s2[-(i + 1)]:
            return i
    return min_len


def count_repeat_right_flank(mut_nucleotides, right_flank_region):
    total_overlap_len = 0
    while right_flank_region.startswith(mut_nucleotides):
        total_overlap_len += len(mut_nucleotides)
        right_flank_region = right_flank_region[len(mut_nucleotides) :]
    total_overlap_len += common_prefix_length(mut_nucleotides, right_flank_region)
    return total_overlap_len


def count_repeat_left_flank(mut_nucleotides, left_flank_region):
    total_overlap_len = 0
    while left_flank_region.endswith(mut_nucleotides):
        total_overlap_len += len(mut_nucleotides)
        left_flank_region = left_flank_region[: -len(mut_nucleotides)]
    total_overlap_len += common_suffix_length(mut_nucleotides, left_flank_region)
    return total_overlap_len


def beginning_mut_nucleotides_with_right_flank(mut_nucleotides, right_flank_region):
    if mut_nucleotides == right_flank_region[: len(mut_nucleotides)]:
        return count_repeat_right_flank(mut_nucleotides, right_flank_region)
    else:
        return common_prefix_length(mut_nucleotides, right_flank_region)


# Comparing end of mut_nucleotides to the end of left_flank_region
def end_mut_nucleotides_with_left_flank(mut_nucleotides, left_flank_region):
    if mut_nucleotides == left_flank_region[-len(mut_nucleotides) :]:
        return count_repeat_left_flank(mut_nucleotides, left_flank_region)
    else:
        return common_suffix_length(mut_nucleotides, left_flank_region)


def calculate_beginning_mutation_overlap_with_right_flank(row):
    if row["mutation_type"] == "deletion":
        sequence_to_check = row["wt_nucleotides_ensembl"]
    else:
        sequence_to_check = row["mut_nucleotides"]

    if row["mutation_type"] == "delins" or row["mutation_type"] == "inversion":
        original_sequence = row["wt_nucleotides_ensembl"] + row["right_flank_region"]
    else:
        original_sequence = row["right_flank_region"]

    return beginning_mut_nucleotides_with_right_flank(
        sequence_to_check, original_sequence
    )


def calculate_end_mutation_overlap_with_left_flank(row):
    if row["mutation_type"] == "deletion":
        sequence_to_check = row["wt_nucleotides_ensembl"]
    else:
        sequence_to_check = row["mut_nucleotides"]

    if row["mutation_type"] == "delins" or row["mutation_type"] == "inversion":
        original_sequence = row["left_flank_region"] + row["wt_nucleotides_ensembl"]
    else:
        original_sequence = row["left_flank_region"]

    return end_mut_nucleotides_with_left_flank(sequence_to_check, original_sequence)


def mutate(
    sequences: Union[str, List[str]],
    mutations: Union[str, List[str]],
    mut_column: str = "mutation",
    seq_id_column: str = "seq_ID",
    mut_id_column: Optional[str] = None,
    k: Optional[int] = None,
    out: Optional[str] = None,
    verbose: bool = True,
    **kwargs,
):
    """
    Takes in nucleotide sequences and mutations (in standard mutation annotation - see below)
    and returns mutated versions of the input sequences according to the provided mutations.

    Reuiqred input argument:
    - sequences     (str) Path to the fasta file containing the sequences to be mutated, e.g., 'seqs.fa'.
                    Sequence identifiers following the '>' character must correspond to the identifiers
                    in the seq_ID column of 'mutations'.

                    Example:
                    >seq1 (or ENSG00000106443)
                    ACTGCGATAGACT
                    >seq2
                    AGATCGCTAG

                    Alternatively: Input sequence(s) as a string or list, e.g. 'AGCTAGCT' or ['ACTGCTAGCT', 'AGCTAGCT'].

                    NOTE: Only the letters until the first space or dot will be used as sequence identifiers
                    - Version numbers of Ensembl IDs will be ignored.
                    NOTE: When 'sequences' input is a genome, also see 'gtf' argument below.
                    
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
                    
                    For more information on the standard mutation annotation, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1867422/.
   
    Additional input arguments:
    - mut_column                   (str) Name of the column containing the mutations to be performed in 'mutations'. Default: 'mutation'.
    - seq_id_column                (str) Name of the column containing the IDs of the sequences to be mutated in 'mutations'. Default: 'seq_ID'.
    - mut_id_column                (str) Name of the column containing the IDs of each mutation in 'mutations'. Default: Will use mut_column.
    - gtf                          (str) Path to .gtf file. When providing a genome fasta file as input for 'sequences', you can provide a .gtf file here 
                                   and the input sequences will be defined according to the transcript boundaries. Default: None
    - gtf_transcript_id_column     (str) Column name in the input 'mutations' file containing the transcript ID. In this case, column seq_id_column should contain the chromosome number.
                                   Required when 'gtf' is provided. Default: None

    Mutant sequence generation/filtering options:
    - k                            (int) Length of sequences flanking the mutation. Default: None (take entire sequence).
                                   If k > total length of the sequence, the entire sequence will be kept.
    - min_seq_len                  (int) Minimum length of the mutant output sequence. Mutant sequences smaller than this will be dropped. 
                                   Default: None
    - optimize_flanking_regions    (True/False) Whether to remove nucleotides from either end of the mutant sequence to ensure (when possible) 
                                   that the mutant sequence does not contain any k-mers also found in the wildtype/input sequence. Default: False
    - remove_seqs_with_wt_kmers    (True/False) Removes output sequences where at least one (k+1)-mer is also present in the wildtype/input sequence in the same region. 
                                   If optimize_flanking_regions=True, only sequences for which a wildtpye kmer is still present after optimization will be removed.
                                   Default: False
    - max_ambiguous                (int) Maximum number of 'N' characters allowed in the output sequence. Default: None (no 'N' filter will be applied)
    - merge_identical              (True/False) Whether to merge identical mutant sequences in the output (identical sequences will be merged by concatenating the sequence 
                                   headers for all identical sequences). Default: True
    - merge_identical_rc           (True/False) Whether to merge identical sequences and their reverse complements in the output. Only effective when merge_identical is also True. Default: True

    # Optional arguments to generate additional output stored in a copy of the 'mutations' DataFrame
    - update_df                    (True/False) Whether to update the input 'mutations' DataFrame to include additional columns with the mutation type, 
                                   wildtype nucleotide sequence, and mutant nucleotide sequence (only valid if 'mutations' is a csv or tsv file). Default: False
    - update_df_out                (str) Path to output csv file containing the updated DataFrame. Only valid if update_df=True.
                                   Default: None -> the new DataFrame will be saved in the same directory as the 'mutations' DataFrame with appendix '_updated'
    - store_full_sequences         (True/False) Whether to also include the complete wildtype and mutant sequences in the updated 'mutations' DataFrame (not just the sub-sequence with 
                                   k-length flanks). Only valid if update_df=True. Default: False
    - translate                    (True/False) Add additional columns to the 'mutations' DataFrame containing the wildtype and mutant amino acid sequences. 
                                   Only valid if store_full_sequences=True. Default: False
    - translate_start              (int | str | None) The position in the input nucleotide sequence to start translating. If a string is provided, it should correspond 
                                   to a column name in 'mutations' containing the open reading frame start positions for each sequence/mutation.
                                   Only valid if translate=True. Default: None (translate from the beginning of the sequence)
    - translate_end                (int | str | None) The position in the input nucleotide sequence to end translating. If a string is provided, it should correspond 
                                   to a column name in 'mutations' containing the open reading frame end positions for each sequence/mutation.
                                   Only valid if translate=True. Default: None (translate from to the end of the sequence)
    
    # General arguments:
    - out                          (str) Path to output fasta file containing the mutated sequences, e.g., 'path/to/output_fasta.fa'.
                                   Default: None -> returns a list of the mutated sequences to standard out.
                                   The identifiers (following the '>') of the mutated sequences in the output fasta will be '>[seq_ID]_[mut_ID]'.
    - verbose                      (True/False) whether to print progress information. Default: True

    Saves mutated sequences in fasta format (or returns a list containing the mutated sequences if out=None).
    """

    if kwargs.get("gtf") or kwargs.get("gtf_transcript_id_column") or kwargs.get("optimize_flanking_regions") or kwargs.get("remove_seqs_with_wt_kmers") or kwargs.get("min_seq_len") or kwargs.get("max_ambiguous") or kwargs.get("merge_identical") or kwargs.get("merge_identical_rc") or kwargs.get("update_df") or kwargs.get("update_df_out") or kwargs.get("store_full_sequences") or kwargs.get("translate") or kwargs.get("translate_start") or kwargs.get("translate_end"):
        # print a log message and raise exception
        logger.critical(
            """
            It appears that you are passing in arguments that are not supported anymore in gget mutate. For use of these arguments, please check out https://github.com/pachterlab/kvar.
            """
        )
        raise NotImplementedError
 
    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations, cosmic_incorrect_wt_base, mut_idx_outside_seq

    columns_to_keep = ["header", seq_id_column, mut_column, "mutation_type", "wt_sequence", "mutant_sequence"]

    if k is None:
        k = 999999999  # take entire sequence by default

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

    mutations_path = None

    # Read in 'mutations' if passed as filepath to comma-separated csv
    if isinstance(mutations, str) and mutations.endswith(".csv"):
        mutations_path = mutations
        mutations = pd.read_csv(mutations)
        for col in mutations.columns:
            if col not in columns_to_keep:
                columns_to_keep.append(col)  # append "mutation_aa", "gene_name", "mutation_id"

    elif isinstance(mutations, str) and mutations.endswith(".tsv"):
        mutations_path = mutations
        mutations = pd.read_csv(mutations, sep="\t")
        for col in mutations.columns:
            if col not in columns_to_keep:
                columns_to_keep.append(col)  # append "mutation_aa", "gene_name", "mutation_id"

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

    number_of_missing_seq_ids = mutations[seq_id_column].isna().sum()

    if number_of_missing_seq_ids > 0:
        logger.warning(
            f"""
            {number_of_missing_seq_ids} rows in 'mutations' are missing sequence IDs. These rows will be dropped from the analysis.
            """
        )

        # Drop rows with missing sequence IDs
        mutations = mutations.dropna(subset=[seq_id_column])

    # ensure seq_ID column is string type, and chromosome numbers don't have decimals
    mutations[seq_id_column] = mutations[seq_id_column].apply(convert_chromosome_value_to_int_when_possible)

    mutations = add_mutation_type(mutations, mut_column)

    # Handle sequences that were not found based on their sequence IDs
    seqs_not_found = mutations[~mutations[seq_id_column].isin(seq_dict.keys())]
    if 0 < len(seqs_not_found) < 20:
        logger.warning(
            f"""
            The sequences with the following {len(seqs_not_found)} sequence ID(s) were not found: {", ".join(seqs_not_found[seq_id_column].values)}  
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
    mutations = mutations.dropna(subset=[seq_id_column, mut_column])
    if len(mutations) < 1:
        raise ValueError(
            """
            None of the input sequences match the sequence IDs provided in 'mutations'. 
            Ensure that the sequence IDs correspond to the string following the > character in the 'sequences' fasta file (do NOT include spaces or dots).
            """
        )

    total_mutations = mutations.shape[0]
    mutations["mutant_sequence"] = ""

    if mut_id_column is not None:
        mutations["header"] = (
            ">" + mutations[mut_id_column]
        )
    else:
        mutations["header"] = (
            ">" + mutations[seq_id_column] + ":" + mutations[mut_column]
        )

    # Calculate number of bad mutations
    uncertain_mutations = mutations[mut_column].str.contains(r"\?").sum()

    ambiguous_position_mutations = mutations[mut_column].str.contains(r"\(|\)").sum()

    intronic_mutations = mutations[mut_column].str.contains(r"\+|\-").sum()

    posttranslational_region_mutations = mutations[mut_column].str.contains(r"\*").sum()

    # Filter out bad mutations
    combined_pattern = re.compile(r"(\?|\(|\)|\+|\-|\*)")
    mask = mutations[mut_column].str.contains(combined_pattern)
    mutations = mutations[~mask]

    # Extract nucleotide positions and mutation info from Mutation CDS
    mutations[["nucleotide_positions", "actual_mutation"]] = mutations[
        mut_column
    ].str.extract(mutation_pattern)

    # Filter out mutations that did not match the re
    unknown_mutations = mutations["nucleotide_positions"].isna().sum()
    mutations = mutations.dropna(subset=["nucleotide_positions", "actual_mutation"])

    if mutations.empty:
        logger.warning("No valid mutations found in the input.")
        return []

    # Split nucleotide positions into start and end positions
    split_positions = mutations["nucleotide_positions"].str.split("_", expand=True)

    mutations["start_mutation_position"] = split_positions[0]
    if split_positions.shape[1] > 1:
        mutations["end_mutation_position"] = split_positions[1].fillna(
            split_positions[0]
        )
    else:
        mutations["end_mutation_position"] = mutations["start_mutation_position"]

    mutations.loc[
        mutations["end_mutation_position"].isna(), "end_mutation_position"
    ] = mutations["start_mutation_position"]

    mutations[["start_mutation_position", "end_mutation_position"]] = mutations[
        ["start_mutation_position", "end_mutation_position"]
    ].astype(int)

    # Adjust positions to 0-based indexing
    mutations["start_mutation_position"] -= 1
    mutations["end_mutation_position"] -= 1  # don't forget to increment by 1 later

    # Calculate sequence length
    mutations["sequence_length"] = mutations[seq_id_column].apply(lambda x: get_sequence_length(x, seq_dict))

    # Filter out mutations with positions outside the sequence
    index_error_mask = (
        mutations["start_mutation_position"] > mutations["sequence_length"]
    ) | (mutations["end_mutation_position"] > mutations["sequence_length"])

    mut_idx_outside_seq = index_error_mask.sum()

    mutations = mutations[~index_error_mask]

    if mutations.empty:
        logger.warning("No valid mutations found in the input.")
        return []

    # Create masks for each type of mutation
    mutations["wt_nucleotides_ensembl"] = None
    substitution_mask = mutations["mutation_type"] == "substitution"
    deletion_mask = mutations["mutation_type"] == "deletion"
    delins_mask = mutations["mutation_type"] == "delins"
    insertion_mask = mutations["mutation_type"] == "insertion"
    duplication_mask = mutations["mutation_type"] == "duplication"
    inversion_mask = mutations["mutation_type"] == "inversion"

    # Create a mask for all non-substitution mutations
    non_substitution_mask = (
        deletion_mask | delins_mask | insertion_mask | duplication_mask | inversion_mask
    )

    # Extract the WT nucleotides for the substitution rows from reference fasta (i.e., Ensembl)
    start_positions = mutations.loc[substitution_mask, "start_mutation_position"].values
    
    # Get the nucleotides at the start positions
    wt_nucleotides_substitution = np.array([
        get_nucleotide_at_position(seq_id, pos, seq_dict)
        for seq_id, pos in zip(mutations.loc[substitution_mask, seq_id_column], start_positions)
    ])

    mutations.loc[substitution_mask, "wt_nucleotides_ensembl"] = (
        wt_nucleotides_substitution
    )

    # Extract the WT nucleotides for the substitution rows from the Mutation CDS (i.e., COSMIC)
    mutations["wt_nucleotides_cosmic"] = None
    mutations.loc[substitution_mask, "wt_nucleotides_cosmic"] = mutations[
        "actual_mutation"
    ].str[0]

    congruent_wt_bases_mask = (
        (mutations["wt_nucleotides_cosmic"] == mutations["wt_nucleotides_ensembl"]) |
        mutations[["wt_nucleotides_cosmic", "wt_nucleotides_ensembl"]].isna().any(axis=1)
    )

    cosmic_incorrect_wt_base = (~congruent_wt_bases_mask).sum()

    mutations = mutations[congruent_wt_bases_mask]

    if mutations.empty:
        logger.warning("No valid mutations found in the input.")
        return []

    # Adjust the start and end positions for insertions
    mutations.loc[
        insertion_mask, "start_mutation_position"
    ] += 1  # in other cases, we want left flank to exclude the start of mutation site; but with insertion, the start of mutation site as it is denoted still belongs in the flank region
    mutations.loc[
        insertion_mask, "end_mutation_position"
    ] -= 1  # in this notation, the end position is one before the start position

    # Extract the WT nucleotides for the non-substitution rows from the Mutation CDS (i.e., COSMIC)
    mutations.loc[non_substitution_mask, "wt_nucleotides_ensembl"] = mutations.loc[
        non_substitution_mask
    ].apply(lambda row: extract_sequence(row, seq_dict, seq_id_column), axis=1)

    # Apply mutations to the sequences
    mutations["mut_nucleotides"] = None
    mutations.loc[substitution_mask, "mut_nucleotides"] = mutations.loc[
        substitution_mask, "actual_mutation"
    ].str[-1]
    mutations.loc[deletion_mask, "mut_nucleotides"] = ""
    mutations.loc[delins_mask, "mut_nucleotides"] = mutations.loc[
        delins_mask, "actual_mutation"
    ].str.extract(r"delins([A-Z]+)")[0]
    mutations.loc[insertion_mask, "mut_nucleotides"] = mutations.loc[
        insertion_mask, "actual_mutation"
    ].str.extract(r"ins([A-Z]+)")[0]
    mutations.loc[duplication_mask, "mut_nucleotides"] = mutations.loc[
        duplication_mask
    ].apply(lambda row: row["wt_nucleotides_ensembl"], axis=1)
    mutations.loc[inversion_mask, "mut_nucleotides"] = mutations.loc[
        inversion_mask
    ].apply(
        lambda row: "".join(
            complement.get(nucleotide, "N")
            for nucleotide in row["wt_nucleotides_ensembl"][::-1]
        ),
        axis=1,
    )

    # Adjust the nucleotide positions of duplication mutations to mimic that of insertions (since duplications are essentially just insertions)
    mutations.loc[duplication_mask, "start_mutation_position"] = (
        mutations.loc[duplication_mask, "end_mutation_position"] + 1
    )  # in the case of duplication, the "mutant" site is still in the left flank as well

    mutations.loc[duplication_mask, "wt_nucleotides_ensembl"] = ""

    # Calculate the kmer bounds
    mutations["start_kmer_position_min"] = mutations["start_mutation_position"] - k
    mutations["start_kmer_position"] = mutations["start_kmer_position_min"].combine(
        0, max
    )

    mutations["end_kmer_position_max"] = mutations["end_mutation_position"] + k
    mutations["end_kmer_position"] = mutations[
        ["end_kmer_position_max", "sequence_length"]
    ].min(
        axis=1
    )  # don't forget to increment by 1 later on

    mut_apply = (lambda *args, **kwargs: mutations.progress_apply(*args, **kwargs)) if verbose else mutations.apply

    if verbose:
        tqdm.pandas(desc="Extracting k-mer left flank sequences")

    mutations["left_flank_region"] = mut_apply(
        lambda row: seq_dict[row[seq_id_column]][
            row["start_kmer_position"] : row["start_mutation_position"]
        ],
        axis=1,
    )  # ? vectorize

    if verbose:
        tqdm.pandas(desc="Extracting k-mer right flank sequences")

    mutations["right_flank_region"] = mut_apply(
        lambda row: seq_dict[row[seq_id_column]][
            row["end_mutation_position"] + 1 : row["end_kmer_position"] + 1
        ],
        axis=1,
    )  # ? vectorize

    mutations["beginning_mutation_overlap_with_right_flank"] = 0
    mutations["end_mutation_overlap_with_left_flank"] = 0

    # Rules for shaving off kmer ends - r1 = left flank, r2 = right flank, d = deleted portion, i = inserted portion
    # Substitution: N/A
    # Deletion:
    # To what extend the beginning of d overlaps with the beginning of r2 --> shave up to that many nucleotides off the beginning of r1 until k - len(r1) ≥ extent of overlap
    # To what extend the end of d overlaps with the beginning of r1 --> shave up to that many nucleotides off the end of r2 until k - len(r2) ≥ extent of overlap
    # Insertion, Duplication:
    # To what extend the beginning of i overlaps with the beginning of r2 --> shave up to that many nucleotides off the beginning of r1 until k - len(r1) ≥ extent of overlap
    # To what extend the end of i overlaps with the beginning of r1 --> shave up to that many nucleotides off the end of r2 until k - len(r2) ≥ extent of overlap
    # Delins, inversion:
    # To what extend the beginning of i overlaps with the beginning of d --> shave up to that many nucleotides off the beginning of r1 until k - len(r1) ≥ extent of overlap
    # To what extend the end of i overlaps with the beginning of d --> shave up to that many nucleotides off the end of r2 until k - len(r2) ≥ extent of overlap

    mutations["updated_left_flank_start"] = 0
    mutations["updated_right_flank_end"] = 0

    # Create WT substitution k-mer sequences
    mutations.loc[substitution_mask, "wt_sequence"] = (
        mutations.loc[substitution_mask, "left_flank_region"]
        + mutations.loc[substitution_mask, "wt_nucleotides_ensembl"]
        + mutations.loc[substitution_mask, "right_flank_region"]
    )

    # Create WT non-substitution k-mer sequences
    mutations.loc[non_substitution_mask, "wt_sequence"] = mutations.loc[
        non_substitution_mask
    ].apply(
        lambda row: row["left_flank_region"][row["updated_left_flank_start"] :]
        + row["wt_nucleotides_ensembl"]
        + row["right_flank_region"][: len(row['right_flank_region']) - row["updated_right_flank_end"]],
        axis=1,
    )

    # Create mutant substitution k-mer sequences
    mutations.loc[substitution_mask, "mutant_sequence"] = (
        mutations.loc[substitution_mask, "left_flank_region"]
        + mutations.loc[substitution_mask, "mut_nucleotides"]
        + mutations.loc[substitution_mask, "right_flank_region"]
    )

    # Create mutant non-substitution k-mer sequences
    mutations.loc[non_substitution_mask, "mutant_sequence"] = mutations.loc[
        non_substitution_mask
    ].apply(
        lambda row: row["left_flank_region"][row["updated_left_flank_start"] :]
        + row["mut_nucleotides"]
        + row["right_flank_region"][: len(row['right_flank_region']) - row["updated_right_flank_end"]],
        axis=1,
    )

    # Calculate k-mer lengths and report the distribution
    mutations["mutant_sequence_kmer_length"] = mutations["mutant_sequence"].apply(
        lambda x: len(x) if pd.notna(x) else 0
    )

    max_length = mutations["mutant_sequence_kmer_length"].max()

    # split_cols = mutations[mut_id_column].str.split("_", n=1, expand=True)

    # if split_cols.shape[1] == 1:
    #     split_cols[1] = None

    # # Extract gene name and mutation ID from mut_ID column (based on formatting of gget cosmic)
    # mutations["gene_name"] = split_cols[0]
    # mutations["mutation_id"] = split_cols[1].fillna(mutations[mut_id_column])

    # Report status of mutations back to user
    good_mutations = mutations.shape[0]

    # good_mutations = total_mutations - intronic_mutations - posttranslational_region_mutations - unknown_mutations - uncertain_mutations - ambiguous_position_mutations - cosmic_incorrect_wt_base - mut_idx_outside_seq

    # if remove_seqs_with_wt_kmers:
    #     good_mutations = good_mutations - long_duplications - mutations_overlapping_with_wt
        
    # if min_seq_len:
    #     good_mutations = good_mutations - rows_less_than_minimum

    # if max_ambigious:
    #     good_mutations = good_mutations - num_rows_with_N

    report = f"""
        {good_mutations} mutations correctly recorded ({good_mutations/total_mutations*100:.2f}%)
        {intronic_mutations} intronic mutations found ({intronic_mutations/total_mutations*100:.2f}%)
        {posttranslational_region_mutations} posttranslational region mutations found ({posttranslational_region_mutations/total_mutations*100:.2f}%)
        {unknown_mutations} unknown mutations found ({unknown_mutations/total_mutations*100:.2f}%)
        {uncertain_mutations} mutations with uncertain mutation found ({uncertain_mutations/total_mutations*100:.2f}%)
        {ambiguous_position_mutations} mutations with ambiguous position found ({ambiguous_position_mutations/total_mutations*100:.2f}%)
        {cosmic_incorrect_wt_base} mutations with incorrect wildtype base found ({cosmic_incorrect_wt_base/total_mutations*100:.2f}%)
        {mut_idx_outside_seq} mutations with indices outside of the sequence length found ({mut_idx_outside_seq/total_mutations*100:.2f}%)
        """

    if good_mutations != total_mutations:
        logger.warning(report)
    else:
        logger.info("All mutations correctly recorded")

    mutations = mutations[columns_to_keep]

    empty_kmer_count = (mutations["mutant_sequence"] == "").sum()
    
    if empty_kmer_count > 0 and verbose:
        logger.warning(
            f"{empty_kmer_count} mutated sequences were empty and were not included in the output."
        )

    mutations = mutations[mutations["mutant_sequence"] != ""]

    mutations['header'] = mutations['header'].str[1:]  # remove the > character

    mutations["fasta_format"] = (
        ">" + mutations["header"] + "\n" + mutations["mutant_sequence"] + "\n"
    )

    if out:
        # Save mutated sequences in new fasta file
        with open(out, "w") as fasta_file:
            fasta_file.write("".join(mutations["fasta_format"].values))

        if verbose:
            logger.info(f"FASTA file containing mutated sequences created at {out}.")

    # When out=None, return list of mutated seqs
    else:
        all_mut_seqs = []
        all_mut_seqs.extend(mutations["mutant_sequence"].values)

        # Remove empty strings from final list of mutated sequences
        # (these are introduced when unknown mutations are encountered)
        while "" in all_mut_seqs:
            all_mut_seqs.remove("")

        if len(all_mut_seqs) > 0:
            return all_mut_seqs
