import os
# import gget
import pandas as pd
from Bio import SeqIO
import re
import base64
import subprocess
import json
import tarfile
import gzip
import shutil
import argparse
from pdb import set_trace as st
import logging
from datetime import datetime
from tqdm import tqdm
tqdm.pandas()

intronic_mutations = 0
posttranslational_region_mutations = 0
unknown_mutations = 0
uncertain_mutations = 0
ambiguous_position_mutations = 0
cosmic_incorrect_wt_base = 0


logging_level = logging.INFO

# logging.basicConfig(level=logging_level, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')   # make sure to switch between logging. and logger.
logger = logging.getLogger(__name__)
logger.setLevel(logging_level)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%H:%M:%S")

console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

log_dir = 'logs'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

log_file = os.path.join(log_dir, f"log_file_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")

file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

mutation_pattern = r'c\.([0-9_\-\+\*]+)([a-zA-Z>]+)'  # more complex: r'c\.([0-9_\-\+\*\(\)\?]+)([a-zA-Z>\(\)0-9]+)'

def raise_pytest_error():
    if os.getenv('TEST_MODE'):
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
        mutation_type = re.findall(r'[a-z>-]+', letters)
        if mutation_type[0] == ">":  # substitution - eg c.1252C>T
            return "substitution"
        
        elif mutation_type[0] == "del":  # deletion - eg c.2126_2128del OR c.1627del
            return "deletion"
        
        elif mutation_type[0] == "delins":  # insertion-deletion - eg c.2239_2253delinsAAT OR c.646delins
            return "delins"
        
        elif mutation_type[0] == "ins":  # insertion - eg c.2239_2240insAAT   (interval is always 1)
            return "insertion"
        
        elif mutation_type[0] == "dup":  # duplication - eg c.2239_2253dup OR c.646dup
            return "duplication"
        
        elif mutation_type[0] == "inv":  # inversion - eg c.2239_2253inv    (always has an interval)
            return "inversion"
        
        else:
            logger.debug(f"mutation {mutation} matches re but is not a known mutation type")
            # raise_pytest_error()
            return "unknown"
    else:
        logger.debug(f"mutation {mutation} does not match re")
        # raise_pytest_error()
        return "unknown"

def load_sequences(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        key = record.description.split()[1]
        sequences[key] = str(record.seq)
    return sequences

def create_mutant_sequence(row, mutation_function):
    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations

    if '?' in row['Mutation CDS']:
        logger.debug(f"Uncertain mutation found in {row['GENOMIC_MUTATION_ID']} - mutation is ignored")
        uncertain_mutations += 1
        return ""
    
    if '(' in row['Mutation CDS']:
        logger.debug(f"Ambiguous mutational position found in {row['GENOMIC_MUTATION_ID']} - mutation is ignored")
        ambiguous_position_mutations += 1
        return ""
    
    match = re.match(mutation_pattern, row['Mutation CDS'])
    if match:
        nucleotide_position = match.group(1)  # The number sequence
        letters = match.group(2)  # The letter sequence

        try:
            if '-' in nucleotide_position or '+' in nucleotide_position:
                logger.debug(f"Intronic nucleotide position found in {row['GENOMIC_MUTATION_ID']} - {nucleotide_position}{letters} - mutation is ignored")
                intronic_mutations += 1
                return ""
            
            elif '*' in nucleotide_position:
                logger.debug(f"Posttranslational region nucleotide position found in {row['GENOMIC_MUTATION_ID']} - {nucleotide_position}{letters} - mutation is ignored")
                posttranslational_region_mutations += 1
                return ""
            
            else:
                if "_" in nucleotide_position:
                    starting_nucleotide_position_index_0 = int(nucleotide_position.split('_')[0]) - 1
                    ending_nucleotide_position_index_0 = int(nucleotide_position.split('_')[1]) - 1
                else:
                    starting_nucleotide_position_index_0 = int(nucleotide_position) - 1
                    ending_nucleotide_position_index_0 = starting_nucleotide_position_index_0
        except:
            logger.debug(f"Error with {row['GENOMIC_MUTATION_ID']} - row['Mutation CDS']")
            unknown_mutations += 1
            raise_pytest_error()
            return ""

        mutant_sequence, adjusted_end_position = mutation_function(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0)

        if not mutant_sequence:
            return ""

        kmer_flanking_length = 30
        kmer_start = max(0, starting_nucleotide_position_index_0 - kmer_flanking_length)
        kmer_end = min(len(mutant_sequence), adjusted_end_position + kmer_flanking_length + 1)
        return str(mutant_sequence[kmer_start:kmer_end])
    else:
        logger.debug(f"Error with {row['GENOMIC_MUTATION_ID']} - {row['Mutation CDS']}")
        unknown_mutations += 1
        raise_pytest_error()
        return ""
    
def substitution_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):  # yields 61-mer
    # assert letters[0] == row['full_sequence'][starting_nucleotide_position_index_0], f"Transcript has {row['full_sequence'][starting_nucleotide_position_index_0]} at position {starting_nucleotide_position_index_0} but mutation is {letters[-1]} at position {starting_nucleotide_position_index_0} in {row['GENOMIC_MUTATION_ID']}"
    global cosmic_incorrect_wt_base
    if letters[0] != row['full_sequence'][starting_nucleotide_position_index_0]:
        logger.debug(f"Transcript has {row['full_sequence'][starting_nucleotide_position_index_0]} at position {starting_nucleotide_position_index_0} but mutation is {letters[-1]} at position {starting_nucleotide_position_index_0} in {row['GENOMIC_MUTATION_ID']}")
        cosmic_incorrect_wt_base += 1    
    mutant_sequence = row['full_sequence'][:starting_nucleotide_position_index_0] + letters[-1] + row['full_sequence'][ending_nucleotide_position_index_0+1:]
    adjusted_end_position = ending_nucleotide_position_index_0
    return mutant_sequence, adjusted_end_position

def deletion_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):  # yields 60-mer
    mutant_sequence = row['full_sequence'][:starting_nucleotide_position_index_0] + row['full_sequence'][ending_nucleotide_position_index_0+1:]
    adjusted_end_position = starting_nucleotide_position_index_0 - 1
    return mutant_sequence, adjusted_end_position

def delins_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):  # yields 60+(length of insertion)-mer
    insertion_string = ''.join(re.findall(r'[A-Z]+', letters))
    mutant_sequence = row['full_sequence'][:starting_nucleotide_position_index_0] + insertion_string + row['full_sequence'][ending_nucleotide_position_index_0+1:]
    adjusted_end_position = starting_nucleotide_position_index_0 + len(insertion_string) - 1
    return mutant_sequence, adjusted_end_position

def insertion_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):  # yields 61+(length of insertion)-mer  - 31 before, length of insertion, 30 after (31 before and 30 after is a little unconventional, but if I want to make it 30 before then I have to return an adjusted start as well as an adjusted end)
    insertion_string = ''.join(re.findall(r'[A-Z]+', letters))
    mutant_sequence = row['full_sequence'][:starting_nucleotide_position_index_0+1] + insertion_string + row['full_sequence'][starting_nucleotide_position_index_0+1:]
    adjusted_end_position = ending_nucleotide_position_index_0 + len(insertion_string) - 1
    return mutant_sequence, adjusted_end_position

def duplication_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):
    insertion_string = row['full_sequence'][starting_nucleotide_position_index_0:ending_nucleotide_position_index_0+1]
    mutant_sequence = row['full_sequence'][:starting_nucleotide_position_index_0] + insertion_string*2 + row['full_sequence'][ending_nucleotide_position_index_0+1:]
    adjusted_end_position = ending_nucleotide_position_index_0 + len(insertion_string)
    return mutant_sequence, adjusted_end_position

def inversion_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):
    insertion_string = row['full_sequence'][starting_nucleotide_position_index_0:ending_nucleotide_position_index_0+1]
    mutant_sequence = row['full_sequence'][:starting_nucleotide_position_index_0] + insertion_string[::-1] + row['full_sequence'][ending_nucleotide_position_index_0+1:]
    adjusted_end_position = ending_nucleotide_position_index_0
    return mutant_sequence, adjusted_end_position

def unknown_mutation(row, letters, starting_nucleotide_position_index_0, ending_nucleotide_position_index_0):
    return "", ""

def download_reference(download_link, tar_folder_path, file_path):
    email = input("Please enter your COSMIC email: ")
    password = input("Please enter your COSMIC password: ")

    # Concatenate the email and password with a colon
    input_string = f"{email}:{password}\n"

    encoded_bytes = base64.b64encode(input_string.encode('utf-8'))
    encoded_string = encoded_bytes.decode('utf-8')
    curl_command = ["curl", "-H", f"Authorization: Basic {encoded_string}", download_link]
    result = subprocess.run(curl_command, capture_output=True, text=True)
    
    response_data = json.loads(result.stdout)
    try:
        true_download_url = response_data.get("url")
    except AttributeError:
        raise AttributeError("Invalid username or password.")

    curl_command2 = ["curl", true_download_url, "--output", f"{tar_folder_path}.tar"]
    result2 = subprocess.run(curl_command2, capture_output=True, text=True)
    
    if result2.returncode != 0:
        raise RuntimeError(f"Failed to download file. Return code: {result.returncode}\n{result.stderr}")
    
    with tarfile.open(f"{tar_folder_path}.tar", 'r') as tar:
        tar.extractall(path=tar_folder_path)
        print(f"Extracted tar file to {tar_folder_path}")

    with gzip.open(f"{file_path}.gz", 'rb') as f_in:
        with open(file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped file to {file_path}")

def select_reference(mutation_class_name, reference_dir, grch_number='37'):
    if mutation_class_name == "transcriptome":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_number}/cosmic/v99/Cosmic_Genes_Fasta_v99_GRCh{grch_number}.tar&bucket=downloads"
        tarred_folder = f"Cosmic_Genes_Fasta_v99_GRCh{grch_number}"
        contained_file = f"Cosmic_Genes_v99_GRCh{grch_number}.fasta"
    
    elif mutation_class_name == "cancer":
        assert grch_number == "37", "CancerMutationCensus data is only available for GRCh37"
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=GRCh{grch_number}/cmc/v99/CancerMutationCensus_AllData_Tsv_v99_GRCh{grch_number}.tar&bucket=downloads"
        tarred_folder = f"CancerMutationCensus_AllData_Tsv_v99_GRCh{grch_number}"
        contained_file = f"CancerMutationCensus_AllData_v99_GRCh{grch_number}.tsv"
        
    elif mutation_class_name == "cell_line":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_number}/cell_lines/v99/CellLinesProject_GenomeScreensMutant_Tsv_v99_GRCh{grch_number}.tar&bucket=downloads"
        tarred_folder = f"CellLinesProject_GenomeScreensMutant_Tsv_v99_GRCh{grch_number}"
        contained_file = f"CellLinesProject_GenomeScreensMutant_v99_GRCh{grch_number}.tsv"

    elif mutation_class_name == "census":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_number}/cosmic/v99/Cosmic_MutantCensus_Tsv_v99_GRCh{grch_number}.tar&bucket=downloads"
        tarred_folder = f"Cosmic_MutantCensus_Tsv_v99_GRCh{grch_number}"
        contained_file = f"Cosmic_MutantCensus_v99_GRCh{grch_number}.tsv"

    elif mutation_class_name == "resistance":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_number}/cosmic/v99/Cosmic_ResistanceMutations_Tsv_v99_GRCh{grch_number}.tar&bucket=downloads"
        tarred_folder = f"Cosmic_ResistanceMutations_Tsv_v99_GRCh{grch_number}"
        contained_file = f"Cosmic_ResistanceMutations_v99_GRCh{grch_number}.tsv"

    elif mutation_class_name == "screen":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_number}/cosmic/v99/Cosmic_GenomeScreensMutant_Tsv_v99_GRCh{grch_number}.tar&bucket=downloads"
        tarred_folder = f"Cosmic_GenomeScreensMutant_Tsv_v99_GRCh{grch_number}"
        contained_file = f"Cosmic_GenomeScreensMutant_v99_GRCh{grch_number}.tsv"

    elif mutation_class_name == "cancer_example":
        download_link = f"https://cog.sanger.ac.uk/cosmic-downloads-production/taster/example_grch{grch_number}.tar"
        tarred_folder = f"example_GRCh{grch_number}"
        contained_file = f"CancerMutationCensus_AllData_v99_GRCh{grch_number}.tsv"
    else:
        raise ValueError(f"Invalid mutation class name. Current options are: transcriptome, cancer, cell_line, census, resistance, screen, cancer_example")
    
    tar_folder_path = os.path.join(reference_dir, tarred_folder)
    file_path = os.path.join(tar_folder_path, contained_file)
    
    if not os.path.exists(file_path):
        if mutation_class_name == "cancer_example":
            curl_command = ["curl", "-L", "--output", f"{tar_folder_path}.tar", download_link]
            result = subprocess.run(curl_command, capture_output=True, text=True)

            with tarfile.open(f"{tar_folder_path}.tar", 'r') as tar:
                tar.extractall(path=tar_folder_path)
                print(f"Extracted tar file to {tar_folder_path}")
        else:
            proceed = input("Reference file not found. Would you like to download it from COSMIC? (this requires an account; free for academic, license for commercial): ").strip().lower()
            if proceed in ['yes', 'y']:
                download_reference(download_link, tar_folder_path, file_path)
            else:
                raise FileNotFoundError(f"Reference file not found at {file_path}. Please download the data or check your file path. Learn more about COSMIC at https://cancer.sanger.ac.uk/cosmic/download/cosmic.")
            
    return(file_path)


def parse_args():
    parser = argparse.ArgumentParser(description="Download and process COSMIC data.")
    parser.add_argument('-m', '--mutation_class_name', type=str, default="cancer", help='Name of the mutation class (default: "cancer")')
    parser.add_argument('-g', '--grch_number', type=str, default="38", help='GRCh number (default: "38")')
    return parser.parse_args()

def main():
    args = parse_args()

    grch_number = str(args.grch_number)  # 38
    mutation_class_name = args.mutation_class_name  # "cancer"  ### Change as needed

    logger.info(f"ARGS: {mutation_class_name}, {grch_number}")
    
    base_dir = "/home/jrich/Desktop/CART_prostate_sc" # os.path.dirname(os.getcwd())   # "/home/jrich/Desktop/CART_prostate_sc"
    reference_dir = os.path.join(base_dir, "reference", f"GRCh{grch_number}")

    if not os.path.exists(reference_dir):
        os.makedirs(reference_dir)

    normal_fasta_file = select_reference("transcriptome", reference_dir, grch_number)  # 5225814 rows
    # normal_fasta_file = os.path.join(reference_dir, f"example_GRCh{grch_number}/fasta/Cosmic_Genes_v99_GRCh{grch_number}.fasta")

    mutation_tsv_file = select_reference(mutation_class_name, reference_dir, grch_number)



    output_fasta_dir = os.path.join(base_dir, 'mutant_reference', f'GRCh{grch_number}')
    if os.getenv('TEST_MODE'):
        output_fasta_dir = os.path.join(output_fasta_dir, 'test')
    output_fasta_file_path = os.path.join(output_fasta_dir, f'{mutation_class_name}_mutant_reference.fasta')

    if not os.path.exists(output_fasta_dir):
        os.makedirs(output_fasta_dir)

    global intronic_mutations, posttranslational_region_mutations, unknown_mutations, uncertain_mutations, ambiguous_position_mutations, cosmic_incorrect_wt_base 

    intronic_mutations = 0
    posttranslational_region_mutations = 0
    unknown_mutations = 0
    uncertain_mutations = 0
    ambiguous_position_mutations = 0
    cosmic_incorrect_wt_base = 0

    if mutation_class_name == "cancer" or mutation_class_name == "cancer_example":
        relevant_cols = ["GENE_NAME", "ACCESSION_NUMBER", "GENOMIC_MUTATION_ID", "MUTATION_URL", "Mutation CDS"]  # "Mutation Description CDS"
    else:
        relevant_cols = ["GENE_SYMBOL", "TRANSCRIPT_ACCESSION", "GENOMIC_MUTATION_ID", "MUTATION_ID", "MUTATION_CDS"]
    logger.info("Reading in mutation dataframe")
    df = pd.read_csv(mutation_tsv_file, usecols=relevant_cols, sep='\t')

    if mutation_class_name != "cancer" and mutation_class_name != "cancer_example":
        df = df.rename(columns={'GENE_SYMBOL': 'GENE_NAME', 'TRANSCRIPT_ACCESSION': 'ACCESSION_NUMBER', 'MUTATION_CDS': 'Mutation CDS'})

    logger.info("Extracting mutation types")
    df['mutation_type'] = df['Mutation CDS'].progress_apply(extract_mutation_type)

    sequences = load_sequences(normal_fasta_file)

    df['full_sequence'] = df['ACCESSION_NUMBER'].map(sequences)

    df['GENOMIC_MUTATION_ID'] = df['GENOMIC_MUTATION_ID'].fillna("NA")

    mutation_types = ["substitution", "deletion", "delins", "insertion", "duplication", "inversion", "unknown"]

    mutation_dict = {mutation_type: [] for mutation_type in mutation_types}

    for mutation_type in mutation_types:
        df_mutation_type = df[df['mutation_type'] == mutation_type].copy()
        df_mutation_type['mutant_sequence_kmer'] = ""
        if mutation_class_name == "cancer" or mutation_class_name == "cancer_example":
            df_mutation_type['MUTATION_ID'] = df_mutation_type['MUTATION_URL'].str.extract(r'id=(\d+)')
        df_mutation_type['header'] = '>' + df_mutation_type['GENE_NAME'] + '_' + df_mutation_type['ACCESSION_NUMBER'] + '_' + df_mutation_type['GENOMIC_MUTATION_ID'] + str(df_mutation_type['MUTATION_ID'])
        mutation_dict[mutation_type] = df_mutation_type

    logger.info("Creating mutation sequences (starting with substitutions)")
    if not mutation_dict['substitution'].empty:
        mutation_dict['substitution']['mutant_sequence_kmer'] = mutation_dict['substitution'].progress_apply(create_mutant_sequence, args=(substitution_mutation,), axis=1)
    if not mutation_dict['deletion'].empty:
        mutation_dict['deletion']['mutant_sequence_kmer'] = mutation_dict['deletion'].progress_apply(create_mutant_sequence, args=(deletion_mutation,), axis=1)
    if not mutation_dict['delins'].empty:
        mutation_dict['delins']['mutant_sequence_kmer'] = mutation_dict['delins'].progress_apply(create_mutant_sequence, args=(delins_mutation,), axis=1)
    if not mutation_dict['insertion'].empty:
        mutation_dict['insertion']['mutant_sequence_kmer'] = mutation_dict['insertion'].progress_apply(create_mutant_sequence, args=(insertion_mutation,), axis=1)
    if not mutation_dict['duplication'].empty:
        mutation_dict['duplication']['mutant_sequence_kmer'] = mutation_dict['duplication'].progress_apply(create_mutant_sequence, args=(duplication_mutation,), axis=1)
    if not mutation_dict['inversion'].empty:
        mutation_dict['inversion']['mutant_sequence_kmer'] = mutation_dict['inversion'].progress_apply(create_mutant_sequence, args=(inversion_mutation,), axis=1)
    if not mutation_dict['unknown'].empty:
        mutation_dict['unknown']['mutant_sequence_kmer'] = mutation_dict['unknown'].progress_apply(create_mutant_sequence, args=(unknown_mutation,), axis=1)

    total_mutations = df.shape[0]
    good_mutations = total_mutations - intronic_mutations - posttranslational_region_mutations - unknown_mutations - uncertain_mutations - ambiguous_position_mutations - cosmic_incorrect_wt_base

    logger.warning(f"""
    {good_mutations} mutations correctly recorded ({good_mutations/total_mutations*100:.2f}%)
    {intronic_mutations} intronic mutations found ({intronic_mutations/total_mutations*100:.2f}%)
    {posttranslational_region_mutations} posttranslational region mutations found ({posttranslational_region_mutations/total_mutations*100:.2f}%)
    {unknown_mutations} unknown mutations found ({unknown_mutations/total_mutations*100:.2f}%)
    {uncertain_mutations} mutations with uncertain mutation found ({uncertain_mutations/total_mutations*100:.2f}%)
    {ambiguous_position_mutations} mutations with ambiguous position found ({ambiguous_position_mutations/total_mutations*100:.2f}%)
    {cosmic_incorrect_wt_base} mutations with incorrect wildtype base found ({cosmic_incorrect_wt_base/total_mutations*100:.2f}%)
    """)

    logger.info("Creating FASTA file")

    with open(output_fasta_file_path, 'w') as fasta_file:
        for mutation in mutation_dict:
            if not mutation_dict[mutation].empty:
                df = mutation_dict[mutation]
                df['fasta_format'] = df['header'] + '\n' + df['mutant_sequence_kmer'] + '\n'  # Ensure there's a newline after each entry
                df = df.dropna(subset=['GENOMIC_MUTATION_ID', 'fasta_format'])
                fasta_file.write("\n".join(df['fasta_format']))

    with open(output_fasta_file_path, 'r') as fasta_file:
        lines = [line for line in fasta_file if line.strip()]

    with open(output_fasta_file_path, 'w') as fasta_file:
        fasta_file.writelines(lines)

    print(f"FASTA file created at {output_fasta_file_path}")

if __name__ == "__main__":
    main()
