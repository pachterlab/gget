import time
import zipfile
import os
import sys
import shutil
import platform
import subprocess
import json
import csv
import pandas as pd
from datetime import datetime
from dateutil import parser

# !!! REMOVE BIOPYTHON DEPENDENCY
from Bio import SeqIO

from .utils import set_up_logger

logger = set_up_logger()

from .compile import PACKAGE_PATH
from .gget_setup import UUID

# Path to precompiled datasets binary
if platform.system() == "Windows":
    PRECOMPILED_DATASETS_PATH = os.path.join(
        PACKAGE_PATH, f"bins/{platform.system()}/datasets.win64.exe"
    )
else:
    PRECOMPILED_DATASETS_PATH = os.path.join(
        PACKAGE_PATH, f"bins/{platform.system()}/datasets"
    )


def run_datasets(
    virus,
    host,
    filename,
    geographic_location,
    annotated,
    complete,
    min_release_date,
    accession,
):
    args_dict = {
        "host": host,
        "filename": filename,
        "geo-location": geographic_location,
        "annotated": annotated,
        "complete-only": complete,
        "released-after": min_release_date,
    }

    # Make datasets binary executable
    if platform.system() != "Windows":
        command = f"chmod +x {PRECOMPILED_DATASETS_PATH}"
        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process_2:
            stderr_2 = process_2.stderr.read().decode("utf-8")
            # Log the standard error if it is not empty
            if stderr_2:
                sys.stderr.write(stderr_2)

        # Return None if the subprocess returned with an error
        if process_2.wait() != 0:
            raise RuntimeError(
                "Making the NCBI 'datasets' binary executable has failed."
            )

    # Initialize the base datasets command
    if accession:
        command = f"{PRECOMPILED_DATASETS_PATH} download virus genome accession {virus} --no-progressbar"
    else:
        command = f"{PRECOMPILED_DATASETS_PATH} download virus genome taxon {virus} --no-progressbar"

    # Loop through the dictionary and construct the command
    for key, value in args_dict.items():
        if value:
            command += f" --{key} '{value}'"

    # Run datasets command and write command output
    start_time = time.time()
    with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process_2:
        stderr_2 = process_2.stderr.read().decode("utf-8")
        # Log the standard error if it is not empty
        if stderr_2:
            sys.stderr.write(stderr_2)

    # Return None if the subprocess returned with an error
    if process_2.wait() != 0:
        raise RuntimeError("NCBI dataset download failed.")
    else:
        logger.debug(
            f"NCBI dataset download complete. Download time: {round(time.time() - start_time, 2)} seconds"
        )


def unzip_file(zip_file_path, extract_to_path):
    """
    Unzips a ZIP file to a specified directory.

    zip_file_path: Path to the ZIP file to unzip.
    extract_to_path: Directory where the contents should be extracted.
    """
    os.makedirs(extract_to_path, exist_ok=True)
    with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
        zip_ref.extractall(extract_to_path)


def load_metadata(jsonl_file):
    """Load metadata from a .jsonl file into a dictionary."""
    metadata_dict = {}
    with open(jsonl_file, "r") as file:
        for line in file:
            metadata = json.loads(line.strip())
            accession = metadata["accession"]
            metadata_dict[accession] = metadata
    return metadata_dict


def parse_date(date_str):
    """Parse various date formats into a standardized datetime object."""
    try:
        date = parser.parse(date_str, default=datetime(1000, 1, 1))
        return date

    except (ValueError, TypeError):
        return None


def filter_sequences(
    fna_file,
    metadata_dict,
    min_seq_length=None,
    max_seq_length=None,
    min_gene_count=None,
    max_gene_count=None,
    nuc_completeness=None,
    host=None,
    host_taxid=None,
    lab_passaged=None,
    geographic_region=None,
    # geographic_location=None,
    submitter_country=None,
    min_collection_date=None,
    max_collection_date=None,
    # annotated=None,
    source_database=None,
    # min_release_date=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
    max_ambiguous_chars=None,
):
    """Filter sequences based on various metadata criteria."""

    # Convert date filters to datetime objects
    min_collection_date = (
        parse_date(min_collection_date) if min_collection_date else None
    )
    max_collection_date = (
        parse_date(max_collection_date) if max_collection_date else None
    )
    # min_release_date = parse_date(min_release_date) if min_release_date else None
    max_release_date = parse_date(max_release_date) if max_release_date else None

    # Read sequences from the .fna file
    filtered_sequences = []
    filtered_metadata = []
    for record in SeqIO.parse(fna_file, "fasta"):
        accession = record.id.split(" ")[0]

        # Check if metadata exists for this accession number
        metadata = metadata_dict.get(accession)
        if metadata is None:
            logger.warning(
                f"No metadata found for sequence {accession}. Sequence will be dropped."
            )
            continue

        # Apply filters
        if min_seq_length is not None or max_seq_length is not None:
            sequence_length = metadata.get("length")
            if sequence_length is None:
                continue  # Skip if metadata is missing
            if min_seq_length is not None and sequence_length < min_seq_length:
                continue
            if max_seq_length is not None and sequence_length > max_seq_length:
                continue

        if min_gene_count is not None or max_gene_count is not None:
            gene_count = metadata.get("geneCount")
            if gene_count is None:
                continue
            if min_gene_count is not None and gene_count < min_gene_count:
                continue
            if max_gene_count is not None and gene_count > max_gene_count:
                continue

        if nuc_completeness is not None:
            completeness_status = metadata.get("completeness")
            if completeness_status is None:
                continue
            if completeness_status.lower() != nuc_completeness.lower():
                continue

        if host is not None:
            host_organism = "_".join(
                metadata.get("host", {}).get("organismName", "").split(" ")
            ).lower()
            if not host_organism:
                continue
            if host_organism != host.lower():
                continue

        if host_taxid is not None:
            host_lineage = metadata.get("host", {}).get("lineage", [])
            if not host_lineage:
                continue
            host_lineage_taxids = {lineage["taxId"] for lineage in host_lineage}
            if host_taxid not in host_lineage_taxids:
                continue

        if lab_passaged is not None:
            from_lab = metadata.get("isLabHost")
            if not from_lab:
                continue

        # if geographic_location is not None:
        #     location = "_".join(
        #         metadata.get("location", {}).get("geographicLocation", "").split(" ")
        #     ).lower()
        #     if not location:
        #         continue
        #     if location != geographic_location.lower():
        #         continue

        if geographic_region is not None:
            location = "_".join(
                metadata.get("location", {}).get("geographicRegion", "").split(" ")
            ).lower()
            if not location:
                continue
            if location != geographic_region.lower():
                continue

        if submitter_country is not None:
            submitter_country_value = "_".join(
                metadata.get("submitter", {}).get("country", "").split(" ")
            ).lower()
            if not submitter_country_value:
                continue
            if submitter_country_value != submitter_country.lower():
                continue

        if min_collection_date is not None or max_collection_date is not None:
            date_str = metadata.get("isolate", {}).get("collectionDate", "")
            date = parse_date(date_str)
            if date_str is None or date is None:
                continue
            if min_collection_date and date < min_collection_date:
                continue
            if max_collection_date and date > max_collection_date:
                continue

        # if annotated is not None:
        #     annotated_value = metadata.get("isAnnotated")
        #     if not annotated_value:
        #         continue

        # if virus_taxid is not None:
        #     virus_lineage = metadata.get("virus", {}).get("lineage", [])
        #     if not virus_lineage:
        #         continue
        #     virus_lineage_taxids = {lineage["taxId"] for lineage in virus_lineage}
        #     if virus_taxid not in virus_lineage_taxids:
        #         continue

        if source_database is not None:
            source_db = metadata.get("sourceDatabase", "").lower()
            if not source_db:
                continue
            if source_db != source_database.lower():
                continue

        # if min_release_date is not None or max_release_date is not None:
        #     release_date_str = metadata.get("releaseDate")
        #     release_date_value = parse_date(release_date_str.split("T")[0])
        #     if release_date_str is None or release_date_value is None:
        #         continue
        #     if min_release_date and release_date_value < min_release_date:
        #         continue
        #     if max_release_date and release_date_value > max_release_date:
        #         continue

        if min_mature_peptide_count is not None or max_mature_peptide_count is not None:
            mature_peptide_count = metadata.get("maturePeptideCount")
            if mature_peptide_count is None:
                continue
            if (
                min_mature_peptide_count is not None
                and mature_peptide_count < min_mature_peptide_count
            ):
                continue
            if (
                max_mature_peptide_count is not None
                and mature_peptide_count > max_mature_peptide_count
            ):
                continue

        if min_protein_count is not None or max_protein_count is not None:
            protein_count = metadata.get("proteinCount")
            if protein_count is None:
                continue
            if min_protein_count is not None and protein_count < min_protein_count:
                continue
            if max_protein_count is not None and protein_count > max_protein_count:
                continue

        # Filter out sequences containing 'N' or 'n' characters
        if max_ambiguous_chars is not None:
            sequence_str = str(record.seq)
            n_count = sequence_str.upper().count("N")

            if n_count > max_ambiguous_chars:
                continue

        filtered_sequences.append(record)
        filtered_metadata.append(metadata)

    num_seqs = len(filtered_sequences)
    if num_seqs > 0:
        logger.info(f"{num_seqs} sequences passed the provided filters.")
        if num_seqs != len(filtered_metadata):
            logger.warning(
                f"Number of sequences ({num_seqs}) and number of metadata entries ({len(filtered_metadata)}) do not match."
            )
        return filtered_sequences, filtered_metadata
    else:
        logger.warning("No sequences passed the provided filters.")
        return None, None


def save_metadata_to_csv(filtered_metadata, output_metadata_file):
    """Save filtered metadata to a CSV file with a specific column order."""

    # Define the column order
    columns = [
        "Accession",
        "Organism Name",
        "GenBank/RefSeq",
        "Submitters",
        "Organization",
        "Submitter Country",
        "Release date",
        "Isolate",
        "Virus Lineage",
        "Length",
        "Nuc Completeness",
        "Geo Location",
        "Host",
        "Host Lineage",
        "Lab Host",
        "Tissue/Specimen/Source",
        "Collection Date",
        "Annotated",
        "SRA Accessions",
        "Bioprojects",
        "Biosample",
    ]

    # Prepare data for DataFrame
    data_for_df = []
    for metadata in filtered_metadata:
        # Grab info on geographic location
        location_info = metadata.get("location", {})
        location_values = [v for v in location_info.values() if v and v != ""]
        location_values.reverse()
        geo_info = ":".join(location_values) if location_values else pd.NA

        # Build table
        row = {
            "Accession": metadata.get("accession", pd.NA),
            "Organism Name": metadata.get("virus", {}).get("organismName", pd.NA),
            "GenBank/RefSeq": metadata.get("sourceDatabase", pd.NA),
            "Submitters": ", ".join(metadata.get("submitter", {}).get("names", [])),
            "Organization": metadata.get("submitter", {}).get("affiliation", pd.NA),
            "Submitter Country": metadata.get("submitter", {}).get("country", ""),
            "Release date": metadata.get("releaseDate", "").split("T")[0],
            "Isolate": metadata.get("isolate", {}).get("name", pd.NA),
            "Virus Lineage": metadata.get("virus", {}).get("lineage", []),
            "Length": metadata.get("length", pd.NA),
            "Nuc Completeness": metadata.get("completeness", pd.NA),
            "Geo Location": geo_info,
            "Host": metadata.get("host", {}).get("organismName", pd.NA),
            "Host Lineage": metadata.get("host", {}).get("lineage", []),
            "Lab Host": metadata.get("labHost", pd.NA),
            "Tissue/Specimen/Source": metadata.get("isolate", {}).get("source", pd.NA),
            "Collection Date": metadata.get("isolate", {}).get("collectionDate", pd.NA),
            "Annotated": metadata.get("isAnnotated", pd.NA),
            "SRA Accessions": metadata.get("sraAccessions", []),
            "Bioprojects": metadata.get("bioprojects", []),
            "Biosample": metadata.get("biosample", pd.NA),
        }
        data_for_df.append(row)

    # Create DataFrame with specified column order
    df = pd.DataFrame(data_for_df, columns=columns)

    # Write DataFrame to CSV
    df.to_csv(output_metadata_file, index=False)


def ncbi_virus(
    virus,
    accession=False,
    outfolder=None,
    host=None,
    min_seq_length=None,
    max_seq_length=None,
    min_gene_count=None,
    max_gene_count=None,
    nuc_completeness=None,
    host_taxid=None,
    lab_passaged=None,
    geographic_region=None,
    geographic_location=None,
    submitter_country=None,
    min_collection_date=None,
    max_collection_date=None,
    annotated=None,
    source_database=None,
    min_release_date=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
    max_ambiguous_chars=None,
    complete=None
):
    """
    Download a virus genome dataset from the NCBI Virus database (https://www.ncbi.nlm.nih.gov/labs/virus/).

    Args:
    - virus                Virus taxon or accession, e.g. 'Norovirus' or 'coronaviridae' or
                           '11320' (taxid of Influenza A virus) or 'NC_045512.2'
                           If this input is a virus NCBI accession (e.g. 'NC_045512.2'), set accession = True.
    - accession            True/False whether 'virus' is an accession. Default: False
    - outfolder            Path to folder to save the requested data in, e.g. 'path/to/norovirus_folder'.
                           Default: None (saves output into current working directory)

    Filters:
    - host                 Host organism, e.g. 'homo sapiens'. Default: None
    - min_seq_length       Min length of the returned sequences, e.g. 6252. Default: None
    - max_seq_length       Max length of the returned sequences, e.g. 7500. Default: None
    - min_gene_count       Min number of genes present in the virus genome, e.g. 1. Default: None
    - max_gene_count       Max number of genes present in the virus genome, e.g. 40. Default: None
    - nuc_completeness     Completeness status of the nucleotide sequence. Should be 'partial' or 'complete'. Default: None
    - host_taxid           NCBI Taxonomy ID of the host organism. Filters the results to only include viruses
                           associated with hosts that fall within the specified TaxID. Default: None
    - lab_passaged         True/False Indicates whether the virus sequence has been passaged in a laboratory setting.
                           Default: None
    - geographic_region    The geographic region where the virus was identified or isolated, e.g. 'Africa' or 'Europe'.
                           Default: None
    - geographic_location  The specific geographic location (e.g., country or area) where the virus was identified or
                           isolated, e.g. 'South_Africa'. Default: None
    - submitter_country    The country of the submitter of the virus sequence data. Filters results by the specified country,
                           e.g. 'South_Africa' or 'Germany'. Default: None
    - min_collection_date  The earliest collection date (in 'YYYY-MM-DD' format). Samples collected before this date will be
                           excluded. Example: '2000-01-01'.  Default: None
    - max_collection_date  The latest collection date (in 'YYYY-MM-DD' format). Samples collected after this date will be
                           excluded. Example: '2014-12-04'. Default: None
    - annotated            True/False Indicates whether the virus genome sequence is annotated. Default: None
    - source_database      The source database from which the virus sequences originate, e.g. 'GenBank'.
                           Default: None
    - min_release_date     The earliest release date (in 'YYYY-MM-DD' format) of the virus sequence data. Sequences released
                           before this date will be excluded. Default: None
    - max_release_date     The latest release date (in 'YYYY-MM-DD' format) of the virus sequence data. Sequences released
                           after this date will be excluded. Default: None
    - min_mature_peptide_count    Min number of mature peptides present in the virus genome. Default: None
    - max_mature_peptide_count    Max number of mature peptides present in the virus genome. Default: None
    - min_protein_count    Min number of proteins present in the virus genome. Default: None
    - max_protein_count    Max number of proteins present in the virus genome. Default: None
    - max_ambiguous_chars  Max number of ambiguous nucleotide characters ('n' or 'N') allowed in the sequences, e.g. 50.
                           Default: None

    Returns a fasta file containing the requested sequences and .csv and .jsonl files containing the associated metadata.
    """

    # Create out and tmp folders
    # current_date = datetime.now().strftime("%Y-%m-%d")
    if outfolder is None:
        outfolder = os.getcwd()
    temp_dir = os.path.join(outfolder, f"tmp_{UUID}")
    os.makedirs(temp_dir, exist_ok=True)

    # Download all sequences matching virus and host from NCBI
    zipped_file = os.path.join(temp_dir, "ncbi_dataset.zip")
    run_datasets(
        virus,
        host,
        zipped_file,
        geographic_location,
        annotated,
        complete,
        min_release_date,
        accession,
    )

    # Unzip NCBI database
    temp_ncbi_folder = os.path.join(temp_dir, "ncbi_dataset")
    unzip_file(zipped_file, temp_ncbi_folder)

    # Define file paths
    fna_file = os.path.join(temp_ncbi_folder, "ncbi_dataset/data/genomic.fna")
    jsonl_file = os.path.join(temp_ncbi_folder, "ncbi_dataset/data/data_report.jsonl")
    output_fasta_file = os.path.join(outfolder, f"{virus}_sequences.fasta")
    output_metadata_csv = os.path.join(outfolder, f"{virus}_metadata.csv")
    output_metadata_jsonl = os.path.join(outfolder, f"{virus}_metadata.jsonl")

    # Load metadata
    metadata_dict = load_metadata(jsonl_file)

    # Define filter criteria (customize these as needed)
    filters = {
        "min_seq_length": min_seq_length,
        "max_seq_length": max_seq_length,
        "min_gene_count": min_gene_count,
        "max_gene_count": max_gene_count,
        "nuc_completeness": nuc_completeness,
        "host": None,  # Host filtering is done when dataset is called
        "host_taxid": host_taxid,
        "lab_passaged": lab_passaged,
        "geographic_region": geographic_region,
        # "geographic_location": geographic_location,
        "submitter_country": submitter_country,
        "min_collection_date": min_collection_date,
        "max_collection_date": max_collection_date,
        # "annotated": annotated,
        "source_database": source_database,
        # "min_release_date": min_release_date,
        "max_release_date": max_release_date,
        "min_mature_peptide_count": min_mature_peptide_count,
        "max_mature_peptide_count": max_mature_peptide_count,
        "min_protein_count": min_protein_count,
        "max_protein_count": max_protein_count,
        "max_ambiguous_chars": max_ambiguous_chars,
    }

    # Filter sequences
    filtered_sequences, filtered_metadata = filter_sequences(
        fna_file, metadata_dict, **filters
    )

    # Save filtered sequences and metadata
    if filtered_sequences:
        # Save filtered sequences to .fa file
        SeqIO.write(filtered_sequences, output_fasta_file, "fasta")

        # Save filtered metadata to .jsonl file
        with open(output_metadata_jsonl, "w") as file:
            for metadata in filtered_metadata:
                file.write(json.dumps(metadata) + "\n")

        # Save filtered metadata to CSV file
        save_metadata_to_csv(filtered_metadata, output_metadata_csv)

    # Delete temporary folder and its concent
    shutil.rmtree(temp_dir)
