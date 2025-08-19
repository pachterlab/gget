# Standard library imports for file operations, regex, JSON handling, and dates
import os
import re
import json
import logging     # For logging level checks
import shutil        # For directory operations  
import traceback     # For error traceback logging
import pandas as pd  # For data manipulation and CSV output
import requests      # For HTTP requests to NCBI API
import zipfile       # For extracting downloaded ZIP files
from datetime import datetime
from dateutil import parser  # For flexible date parsing

# Internal imports for logging, unique ID generation, and FASTA parsing
from .utils import set_up_logger, FastaIO

# Set up logger for this module
logger = set_up_logger()

from .gget_setup import UUID

# NCBI Datasets API base URL - Version 2 API endpoint
NCBI_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"


def fetch_virus_metadata(
    virus,
    accession=False,
    host=None,
    geographic_location=None,
    annotated=None,
    complete_only=None,
    min_release_date=None,
    refseq_only=False,
):
    """
    Fetch virus metadata using NCBI Datasets API.
    
    This function retrieves metadata for virus sequences from the NCBI Datasets API
    using either taxon-based or accession-based queries. It handles pagination
    automatically to retrieve all available results.
    
    Args:
        virus (str): Virus taxon name/ID or accession number
        accession (bool): Whether virus parameter is an accession number
        host (str): Host organism name filter
        geographic_location (str): Geographic location filter  
        annotated (bool): Filter for annotated genomes only
        complete_only (bool): Filter for complete genomes only  
        min_release_date (str): Minimum release date filter (YYYY-MM-DD format)
        refseq_only (bool): Limit to RefSeq genomes only
        
    Returns:
        list: List of virus metadata records from the API response
        
    Raises:
        RuntimeError: If the API request fails
    """
    
    # Choose the appropriate API endpoint based on whether we're querying by accession or taxon
    if accession:
        # For accession numbers (e.g., NC_045512.2), use the accession-specific endpoint
        url = f"{NCBI_API_BASE}/virus/accession/{virus}/dataset_report"
        logger.debug("Using accession endpoint for virus: %s", virus)
        params = {}
    else:
        # For taxon names/IDs (e.g., 'SARS-CoV-2', 'coronaviridae'), use the taxon endpoint  
        url = f"{NCBI_API_BASE}/virus/taxon/{virus}/dataset_report"
        logger.debug("Using taxon endpoint for virus: %s", virus)
        params = {}
    
    # Add API-level filters to reduce the amount of data we need to download
    # These filters are applied server-side before results are returned
    if refseq_only:
        # Limit results to RefSeq database entries only (higher quality, curated sequences)
        params['filter.refseq_only'] = 'true'
        logger.debug("Applied RefSeq-only filter")
    
    if annotated is True:
        # Only return sequences that have been annotated with gene/protein information
        params['filter.annotated_only'] = 'true'
        logger.debug("Applied annotated-only filter")
    
    if complete_only:
        # Only return complete genome sequences (not partial sequences)
        params['filter.complete_only'] = 'true'
        logger.debug("Applied complete-only filter")
        
    if host:
        # Filter by host organism name, replacing underscores with spaces for API compatibility
        params['filter.host'] = host.replace('_', ' ')
        logger.debug("Applied host filter: %s", host)
        
    if geographic_location:
        # Filter by geographic location, replacing underscores with spaces for API compatibility
        params['filter.geo_location'] = geographic_location.replace('_', ' ')
        logger.debug("Applied geographic location filter: %s", geographic_location)
        
    if min_release_date:
        # Convert date to ISO format expected by the API (YYYY-MM-DDTHH:MM:SS.sssZ)
        params['filter.released_since'] = f"{min_release_date}T00:00:00.000Z"
        logger.debug("Applied minimum release date filter: %s", min_release_date)
    
    # Note: max_release_date filtering will be done post-API call since the API doesn't support it
    
    # Set page size to maximum allowed to minimize the number of API calls needed
    # The NCBI API supports pagination for large result sets
    params['page_size'] = 1000
    logger.debug("Set page size to maximum: 1000 records per request")
    
    # Initialize variables for handling paginated results
    all_reports = []      # Will store all metadata records across all pages
    page_token = None     # Token for accessing subsequent pages
    page_count = 0        # Track number of pages processed for logging
    
    # Main pagination loop - continue until all pages are retrieved
    while True:
        page_count += 1
        logger.debug("Fetching page %d of results...", page_count)
        
        # Add pagination token if we're not on the first page
        if page_token:
            params['page_token'] = page_token
            
        try:
            # Make the HTTP GET request to the NCBI API
            logger.debug("Making API request to: %s", url)
            logger.debug("Request parameters: %s", params)
            response = requests.get(url, params=params, timeout=30)
            
            # Raise an exception if the HTTP request failed (4xx or 5xx status codes)
            response.raise_for_status()
            
            # Parse the JSON response
            data = response.json()
            logger.debug("Received response with %d bytes", len(response.content))
            
            # Extract the virus reports from the response
            reports = data.get('reports', [])
            logger.debug("Page %d contains %d virus records", page_count, len(reports))
            
            # Add this page's reports to our complete collection
            all_reports.extend(reports)
            
            # Check if there are more pages to retrieve
            next_page_token = data.get('next_page_token')
            if not next_page_token:
                logger.debug("No more pages available, pagination complete")
                break
            
            # Set up for the next page
            page_token = next_page_token
            logger.debug("Next page token received, continuing pagination...")
            
        except requests.exceptions.Timeout as e:
            # Handle timeout errors specifically
            raise RuntimeError(f"Request timed out while fetching virus metadata: {e}") from e
        except requests.exceptions.ConnectionError as e:
            # Handle connection errors (network issues, DNS failures, etc.)
            raise RuntimeError(f"Connection error while fetching virus metadata: {e}") from e
        except requests.exceptions.HTTPError as e:
            # Handle HTTP errors (404, 500, etc.)
            raise RuntimeError(f"HTTP error while fetching virus metadata: {e}") from e
        except requests.exceptions.RequestException as e:
            # Handle any other request-related errors
            raise RuntimeError(f"Failed to fetch virus metadata: {e}") from e
        except json.JSONDecodeError as e:
            # Handle JSON parsing errors
            raise RuntimeError(f"Failed to parse API response as JSON: {e}") from e
    
    # Log the final results summary
    logger.info("Successfully retrieved %d virus records from NCBI API across %d pages", 
                len(all_reports), page_count)
    
    return all_reports


def download_virus_sequences(
    virus,
    accession=False,
    host=None,
    geographic_location=None,
    annotated=None,
    complete_only=None,
    min_release_date=None,
    outdir=None,
):
    """
    Download virus genome sequences using NCBI Datasets API.
    
    This function downloads the actual sequence data (FASTA files) for viruses
    matching the specified criteria. The sequences are returned as a ZIP archive
    containing FASTA files and associated metadata.
    
    Args:
        virus (str): Virus taxon name/ID or accession number
        accession (bool): Whether virus parameter is an accession number
        host (str): Host organism name filter
        geographic_location (str): Geographic location filter
        annotated (bool): Filter for annotated genomes only
        complete_only (bool): Filter for complete genomes only
        min_release_date (str): Minimum release date filter (YYYY-MM-DD format)
        outdir (str): Output directory for downloaded files
        
    Returns:
        str: Path to the downloaded ZIP file containing sequences and metadata
        
    Raises:
        RuntimeError: If the download request fails
    """
    
    # Choose the appropriate download endpoint based on query type
    if accession:
        # For accession-based downloads, use the accession-specific endpoint
        url = f"{NCBI_API_BASE}/virus/accession/{virus}/genome/download"
        logger.debug("Using accession download endpoint for virus: %s", virus)
        
        # For accession downloads, use GET with query parameters
        params = {}
        
        # Apply the same filters as used in metadata fetching
        if host:
            # Filter by host organism, replacing underscores with spaces
            params['filter.host'] = host.replace('_', ' ')
            logger.debug("Applied host filter for download: %s", host)
        
        if geographic_location:
            # Filter by geographic location, replacing underscores with spaces
            params['filter.geo_location'] = geographic_location.replace('_', ' ')
            logger.debug("Applied geographic location filter for download: %s", geographic_location)
            
        if annotated is True:
            # Only download sequences that have annotation data
            params['filter.annotated_only'] = 'true'
            logger.debug("Applied annotated-only filter for download")
            
        if complete_only:
            # Only download complete genome sequences
            params['filter.complete_only'] = 'true'
            logger.debug("Applied complete-only filter for download")
            
        if min_release_date:
            # Filter by minimum release date in ISO format
            params['filter.released_since'] = f"{min_release_date}T00:00:00.000Z"
            logger.debug("Applied minimum release date filter for download: %s", min_release_date)
        
        # Include genomic sequences in the download package
        # This ensures we get the actual FASTA sequences, not just metadata
        params['include_annotation_type'] = 'GENOME_FASTA'
        logger.debug("Requesting genomic FASTA sequences in download")
        
        try:
            logger.info("Initiating download request for accession: %s", virus)
            # Make the download request with extended timeout for large files
            response = requests.get(url, params=params, timeout=300)
            
            # Check if the request was successful
            response.raise_for_status()
            
            # Determine output directory - use current working directory if not specified
            if not outdir:
                outdir = os.getcwd()
                logger.debug("No output directory specified, using current directory: %s", outdir)
            
            # Create a unique filename for the downloaded ZIP file
            zip_path = os.path.join(outdir, f"ncbi_virus_{virus}_{UUID}.zip")
            logger.debug("Saving download to: %s", zip_path)
            
            # Write the downloaded content to the ZIP file
            with open(zip_path, 'wb') as f:
                f.write(response.content)
                
            logger.info("Successfully downloaded virus dataset to: %s (size: %d bytes)", 
                       zip_path, len(response.content))
            return zip_path
            
        except requests.exceptions.Timeout as e:
            # Handle timeout errors - downloads can take a long time for large datasets
            raise RuntimeError(f"Download request timed out for virus sequences: {e}") from e
        except requests.exceptions.ConnectionError as e:
            # Handle connection errors during download
            raise RuntimeError(f"Connection error during virus sequence download: {e}") from e
        except requests.exceptions.HTTPError as e:
            # Handle HTTP errors during download
            raise RuntimeError(f"HTTP error during virus sequence download: {e}") from e
        except requests.exceptions.RequestException as e:
            # Handle any other request-related errors during download
            raise RuntimeError(f"Failed to download virus sequences: {e}") from e
        except IOError as e:
            # Handle file I/O errors when saving the ZIP file
            raise RuntimeError(f"Failed to save downloaded file: {e}") from e
    
    else:
        # For taxon-based downloads, use the taxon-specific endpoint
        url = f"{NCBI_API_BASE}/virus/taxon/{virus}/genome/download"
        logger.debug("Using taxon download endpoint for virus: %s", virus)
        
        # Build query parameters with the same filters
        params = {}
        
        if host:
            # Filter by host organism name
            params['filter.host'] = host.replace('_', ' ')
            logger.debug("Applied host filter for download: %s", host)
        
        if geographic_location:
            # Filter by geographic location
            params['filter.geo_location'] = geographic_location.replace('_', ' ')
            logger.debug("Applied geographic location filter for download: %s", geographic_location)
            
        if annotated is True:
            # Only download annotated sequences
            params['filter.annotated_only'] = 'true'
            logger.debug("Applied annotated-only filter for download")
            
        if complete_only:
            # Only download complete genomes
            params['filter.complete_only'] = 'true'
            logger.debug("Applied complete-only filter for download")
            
        if min_release_date:
            # Filter by minimum release date
            params['filter.released_since'] = f"{min_release_date}T00:00:00.000Z"
            logger.debug("Applied minimum release date filter for download: %s", min_release_date)
        
        # Include genomic sequences in the download package
        params['include_annotation_type'] = 'GENOME_FASTA'
        logger.debug("Requesting genomic FASTA sequences in download")
        
        try:
            logger.info("Initiating download request for taxon: %s", virus)
            # Make the download request with extended timeout
            response = requests.get(url, params=params, timeout=300)
            
            # Check if the request was successful
            response.raise_for_status()
            
            # Determine output directory
            if not outdir:
                outdir = os.getcwd()
                logger.debug("No output directory specified, using current directory: %s", outdir)
            
            # Create a unique filename for the downloaded ZIP file
            zip_path = os.path.join(outdir, f"ncbi_virus_{virus}_{UUID}.zip")
            logger.debug("Saving download to: %s", zip_path)
            
            # Write the downloaded content to the ZIP file
            with open(zip_path, 'wb') as f:
                f.write(response.content)
                
            logger.info("Successfully downloaded virus dataset to: %s (size: %d bytes)", 
                       zip_path, len(response.content))
            return zip_path
            
        except requests.exceptions.Timeout as e:
            # Handle timeout errors during download
            raise RuntimeError(f"Download request timed out for virus sequences: {e}") from e
        except requests.exceptions.ConnectionError as e:
            # Handle connection errors during download
            raise RuntimeError(f"Connection error during virus sequence download: {e}") from e
        except requests.exceptions.HTTPError as e:
            # Handle HTTP errors during download
            raise RuntimeError(f"HTTP error during virus sequence download: {e}") from e
        except requests.exceptions.RequestException as e:
            # Handle any other request-related errors during download
            raise RuntimeError(f"Failed to download virus sequences: {e}") from e
        except IOError as e:
            # Handle file I/O errors when saving the ZIP file
            raise RuntimeError(f"Failed to save downloaded file: {e}") from e


def download_sequences_by_accessions(accessions, outdir=None):
    """
    Download virus genome sequences for a specific list of accession numbers.
    
    This function downloads sequences for a pre-filtered list of accessions,
    using NCBI E-utilities API since the Datasets API virus endpoint only
    provides metadata, not actual sequence data.
    
    Args:
        accessions (list): List of accession numbers to download
        outdir (str): Output directory for downloaded files
        
    Returns:
        str: Path to the downloaded FASTA file containing sequences
        
    Raises:
        RuntimeError: If the download request fails
        ValueError: If no accessions are provided
    """
    
    if not accessions:
        raise ValueError("No accessions provided for download")
    
    logger.info("Downloading sequences for %d accessions using E-utilities API", len(accessions))
    logger.debug("Accession list: %s", accessions[:5] + ['...'] if len(accessions) > 5 else accessions)
    
    # Determine output directory - use current working directory if not specified
    if not outdir:
        outdir = os.getcwd()
        logger.debug("No output directory specified, using current directory: %s", outdir)
    
    # Create output FASTA file path
    fasta_path = os.path.join(outdir, f"virus_sequences_{UUID}.fasta")
    logger.debug("Saving sequences to: %s", fasta_path)
    
    # Use NCBI E-utilities to fetch FASTA sequences
    # This is more reliable for getting actual sequence data
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    # Build accession string (E-utils supports comma-separated IDs)
    accession_string = ",".join(accessions)
    
    params = {
        'db': 'nucleotide',
        'id': accession_string,
        'rettype': 'fasta',
        'retmode': 'text'
    }
    
    try:
        logger.info("Initiating E-utilities request for %d accessions", len(accessions))
        logger.debug("E-utilities URL: %s", base_url)
        
        # Make the request with extended timeout for large datasets
        response = requests.get(base_url, params=params, timeout=300)
        
        # Check if the request was successful
        response.raise_for_status()
        
        # Verify we got FASTA data
        if not response.text.strip().startswith('>'):
            raise RuntimeError(f"Invalid FASTA response: {response.text[:100]}")
        
        # Count sequences in response
        sequence_count = response.text.count('>')
        logger.info("Received %d sequences from E-utilities", sequence_count)
        
        # Write FASTA data to file
        with open(fasta_path, 'w', encoding='utf-8') as f:
            f.write(response.text)
            
        logger.info("Successfully saved sequences to: %s (%.2f MB)", 
                   fasta_path, len(response.text.encode('utf-8')) / 1024 / 1024)
        
        return fasta_path
        
    except requests.exceptions.RequestException as e:
        logger.error("E-utilities request failed: %s", e)
        raise RuntimeError(f"Failed to download virus sequences via E-utilities: {e}") from e
    except IOError as e:
        logger.error("Failed to save FASTA file: %s", e)
        raise RuntimeError(f"Failed to save downloaded sequences: {e}") from e


def unzip_file(zip_file_path, extract_to_path):
    """
    Unzips a ZIP file to a specified directory.
    
    This function extracts all contents from a ZIP archive to the specified
    destination directory. It creates the destination directory if it doesn't exist.

    Args:
        zip_file_path (str): Path to the ZIP file to extract
        extract_to_path (str): Directory where the contents should be extracted
        
    Raises:
        zipfile.BadZipFile: If the ZIP file is corrupted or invalid
        PermissionError: If there are insufficient permissions to create directories or files
        FileNotFoundError: If the ZIP file doesn't exist
    """
    # Create the extraction directory if it doesn't exist
    # exist_ok=True prevents errors if the directory already exists
    os.makedirs(extract_to_path, exist_ok=True)
    logger.debug("Created extraction directory: %s", extract_to_path)
    
    try:
        # Open the ZIP file in read mode
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            # Extract all contents to the specified directory
            zip_ref.extractall(extract_to_path)
            
            # Log information about the extracted contents
            file_list = zip_ref.namelist()
            logger.info("Successfully extracted %d files from %s to %s", 
                       len(file_list), zip_file_path, extract_to_path)
            logger.debug("Extracted files: %s", file_list[:10])  # Log first 10 files
            
    except zipfile.BadZipFile as e:
        raise zipfile.BadZipFile(f"Invalid or corrupted ZIP file: {zip_file_path}") from e
    except PermissionError as e:
        raise PermissionError(f"Permission denied when extracting to: {extract_to_path}") from e
    except FileNotFoundError as e:
        raise FileNotFoundError(f"ZIP file not found: {zip_file_path}") from e


def load_metadata_from_api_reports(api_reports):
    """
    Load metadata from API response reports into a dictionary.
    
    This function transforms the raw API response format into a standardized
    internal metadata format that can be used by the filtering functions.
    It maps API field names to the expected internal field names and handles
    missing or null values appropriately.
    
    Args:
        api_reports (list): List of virus metadata reports from the NCBI API
        
    Returns:
        dict: Dictionary mapping accession numbers to metadata dictionaries
              Key: accession number (str)
              Value: metadata dictionary with standardized field names
    """
    metadata_dict = {}
    processed_count = 0
    skipped_count = 0
    
    logger.debug("Processing %d API reports into metadata dictionary", len(api_reports))
    
    for report in api_reports:
        # Extract the accession number - this serves as our unique identifier
        accession = report.get("accession")
        
        if accession:
            processed_count += 1
            
            # Transform API report format to match expected internal metadata format
            # Map API fields to expected internal field names with appropriate defaults
            metadata = {
                # Basic sequence information
                "accession": accession,
                "length": report.get("length"),  # Sequence length in nucleotides
                "geneCount": report.get("gene_count"),  # Number of genes annotated
                
                # Sequence completeness status (complete, partial, etc.)
                "completeness": report.get("completeness_status", "").lower(),
                
                # Host organism information
                "host": report.get("host", {}),  # Host organism details
                "isLabHost": report.get("host", {}).get("is_lab_host", False),  # Lab-passaged flag
                "labHost": report.get("host", {}).get("is_lab_host", False),  # Alternative field name
                
                # Geographic and location information
                "location": report.get("geo_location", {}),  # Geographic location details
                
                # Submission and database information
                "submitter": report.get("submitters", [{}])[0] if report.get("submitters") else {},
                "sourceDatabase": report.get("source_database", ""),  # GenBank, RefSeq, etc.
                
                # Isolate and virus information
                "isolate": report.get("isolate", {}),  # Sample/isolate details
                "virus": report.get("virus", {}),  # Virus taxonomy and classification
                
                # Annotation status
                "isAnnotated": report.get("is_annotated", False),  # Whether sequence is annotated
                
                # Release and submission dates
                "releaseDate": report.get("release_date", ""),  # When sequence was released
                
                # Associated database records
                "sraAccessions": report.get("sra_accessions", []),  # SRA read data accessions
                "bioprojects": report.get("bioprojects", []),  # Associated BioProject IDs
                "biosample": report.get("biosample"),  # BioSample ID
                
                # Protein and peptide counts
                "proteinCount": report.get("protein_count"),  # Number of proteins
                "maturePeptideCount": report.get("mature_peptide_count"),  # Number of mature peptides
            }
            
            # Store the metadata using accession as the key
            metadata_dict[accession] = metadata
            
            logger.debug("Processed metadata for accession: %s (length: %s, host: %s)", 
                        accession, 
                        metadata.get("length"), 
                        metadata.get("host", {}).get("organismName", "Unknown"))
            
        else:
            # Skip reports without accession numbers
            skipped_count += 1
            logger.warning("Skipping API report without accession number: %s", report)
    
    logger.info("Processed %d metadata records, skipped %d records without accessions", 
                processed_count, skipped_count)
    
    return metadata_dict



def parse_date(date_str, filtername="", verbose=False):
    """
    Parse various date formats into a standardized datetime object.
    
    This function uses the dateutil parser to handle various date formats
    that might be encountered in virus metadata. It provides helpful error
    messages when date parsing fails.
    
    Args:
        date_str (str): Date string to parse (various formats accepted)
        filtername (str): Name of the filter/field for error reporting
        verbose (bool): Whether to raise detailed exceptions on parse errors
        
    Returns:
        datetime: Parsed datetime object, or None if parsing fails (when verbose=False)
        
    Raises:
        ValueError: If date parsing fails and verbose=True
        
    Note:
        Uses a default date of year 1000 for incomplete date strings to ensure
        proper comparison behavior with minimum date filters.
    """
    try:
        # Use dateutil parser for flexible date parsing
        # Default to year 1000 for partial dates to handle edge cases properly
        parsed_date = parser.parse(date_str, default=datetime(1000, 1, 1))
        logger.debug("Successfully parsed date '%s' as %s", date_str, parsed_date)
        return parsed_date
        
    except (ValueError, TypeError) as exc:
        # Handle parsing errors based on verbosity setting
        if verbose:
            # In verbose mode, raise detailed error with helpful message
            error_msg = (
                f"Invalid date detected for argument {filtername}: '{date_str}'.\n"
                "Note: Please check for errors such as incorrect day values "
                "(e.g., June 31st does not exist) or typos in the date format "
                "(should be YYYY-MM-DD)."
            )
            logger.error("Date parsing failed: %s", error_msg)
            raise ValueError(error_msg) from exc
        else:
            # In non-verbose mode, log warning and return None
            logger.warning("Failed to parse date '%s' for filter '%s': %s", 
                          date_str, filtername, exc)
            return None


# !!! TO-DO: Find ways to filter pre-download so you don't have to download a ton of seqs just to filter after
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
    annotated=None,
    source_database=None,
    # min_release_date=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
    max_ambiguous_chars=None,
    has_proteins=None,
    proteins_complete=False,
):
    """
    Filter sequences based on various metadata criteria.
    
    This function applies post-download filtering to virus sequences based on
    metadata criteria that couldn't be applied at the API level. It reads
    sequences from a FASTA file and applies multiple filters sequentially.
    
    Args:
        fna_file (str): Path to the FASTA file containing virus sequences
        metadata_dict (dict): Dictionary mapping accession numbers to metadata
        min_seq_length (int): Minimum sequence length filter
        max_seq_length (int): Maximum sequence length filter
        min_gene_count (int): Minimum number of genes filter
        max_gene_count (int): Maximum number of genes filter
        nuc_completeness (str): Completeness status filter ('complete' or 'partial')
        host (str): Host organism name filter
        host_taxid (str): Host taxonomy ID filter
        lab_passaged (bool): Lab passaging status filter
        geographic_region (str): Geographic region filter
        submitter_country (str): Submitter country filter
        min_collection_date (str): Minimum collection date filter
        max_collection_date (str): Maximum collection date filter
        annotated (bool): Annotation status filter
        source_database (str): Source database filter
        max_release_date (str): Maximum release date filter
        min_mature_peptide_count (int): Minimum mature peptide count filter
        max_mature_peptide_count (int): Maximum mature peptide count filter
        min_protein_count (int): Minimum protein count filter
        max_protein_count (int): Maximum protein count filter
        max_ambiguous_chars (int): Maximum number of 'N' characters allowed
        has_proteins (str/list): Required proteins/genes/segments
        proteins_complete (bool): Whether required proteins must be marked complete
        
    Returns:
        tuple: (filtered_sequences, filtered_metadata, protein_headers)
               - filtered_sequences: List of FastaRecord objects that passed filters
               - filtered_metadata: List of metadata dictionaries for filtered sequences
               - protein_headers: List of protein/segment information from headers
               Returns (None, None, None) if no sequences pass the filters
    """

    logger.info("Starting sequence filtering process...")
    logger.debug("Applying filters: seq_length(%s-%s), gene_count(%s-%s), completeness(%s), "
                "host(%s), host_taxid(%s), lab_passaged(%s), geo_region(%s), "
                "submitter_country(%s), collection_date(%s-%s), annotated(%s), "
                "source_db(%s), max_release_date(%s), peptide_count(%s-%s), "
                "protein_count(%s-%s), max_ambiguous(%s), has_proteins(%s), proteins_complete(%s)",
                min_seq_length, max_seq_length, min_gene_count, max_gene_count,
                nuc_completeness, host, host_taxid, lab_passaged, geographic_region,
                submitter_country, min_collection_date, max_collection_date,
                annotated, source_database, max_release_date, min_mature_peptide_count,
                max_mature_peptide_count, min_protein_count, max_protein_count,
                max_ambiguous_chars, has_proteins, proteins_complete)

    # Convert date filters to datetime objects for proper comparison
    # Use verbose=True to provide detailed error messages for date parsing failures
    min_collection_date = (
        parse_date(min_collection_date, filtername="min_collection_date", verbose=True) 
        if min_collection_date else None
    )
    max_collection_date = (
        parse_date(max_collection_date, filtername="max_collection_date", verbose=True) 
        if max_collection_date else None
    )
    max_release_date = (
        parse_date(max_release_date, filtername="max_release_date", verbose=True) 
        if max_release_date else None
    )
    
    if min_collection_date:
        logger.debug("Parsed min_collection_date: %s", min_collection_date)
    if max_collection_date:
        logger.debug("Parsed max_collection_date: %s", max_collection_date)
    if max_release_date:
        logger.debug("Parsed max_release_date: %s", max_release_date)

    # Initialize lists to store filtered results
    filtered_sequences = []    # Will store FastaRecord objects that pass filters
    filtered_metadata = []     # Will store corresponding metadata dictionaries
    protein_headers = []       # Will store protein/segment information from FASTA headers
    
    # Counters for logging filter statistics
    total_sequences = 0
    sequences_with_metadata = 0
    filter_stats = {
        'no_metadata': 0,
        'seq_length': 0,
        'gene_count': 0,
        'completeness': 0,
        'host': 0,
        'host_taxid': 0,
        'lab_passaged': 0,
        'geographic_region': 0,
        'submitter_country': 0,
        'collection_date': 0,
        'annotated': 0,
        'source_database': 0,
        'release_date': 0,
        'mature_peptide_count': 0,
        'protein_count': 0,
        'ambiguous_chars': 0,
        'has_proteins': 0,
    }

    # Read and process sequences from the FASTA file
    logger.info("Reading sequences from FASTA file: %s", fna_file)
    for record in FastaIO.parse(fna_file, "fasta"):
        total_sequences += 1
        
        # Extract accession number from the sequence ID (first part before any spaces)
        accession = record.id.split(" ")[0]
        logger.debug("Processing sequence %d: %s", total_sequences, accession)

        # Check if metadata exists for this accession number
        metadata = metadata_dict.get(accession)
        if metadata is None:
            filter_stats['no_metadata'] += 1
            logger.warning("No metadata found for sequence %s. Sequence will be dropped.", accession)
            continue
        
        sequences_with_metadata += 1

        # Apply filters sequentially - each filter can exclude the sequence
        # If any filter fails, we continue to the next sequence
        
        # FILTER 1: Sequence length filters
        if min_seq_length is not None or max_seq_length is not None:
            sequence_length = metadata.get("length")
            if sequence_length is None:
                logger.debug("Skipping %s: missing length metadata", accession)
                filter_stats['seq_length'] += 1
                continue  # Skip if length metadata is missing
                
            if min_seq_length is not None and sequence_length < min_seq_length:
                logger.debug("Skipping %s: length %d < min %d", accession, sequence_length, min_seq_length)
                filter_stats['seq_length'] += 1
                continue
                
            if max_seq_length is not None and sequence_length > max_seq_length:
                logger.debug("Skipping %s: length %d > max %d", accession, sequence_length, max_seq_length)
                filter_stats['seq_length'] += 1
                continue

        # FILTER 2: Gene count filters
        if min_gene_count is not None or max_gene_count is not None:
            gene_count = metadata.get("geneCount")
            if gene_count is None:
                logger.debug("Skipping %s: missing gene count metadata", accession)
                filter_stats['gene_count'] += 1
                continue
                
            if min_gene_count is not None and gene_count < min_gene_count:
                logger.debug("Skipping %s: gene count %d < min %d", accession, gene_count, min_gene_count)
                filter_stats['gene_count'] += 1
                continue
                
            if max_gene_count is not None and gene_count > max_gene_count:
                logger.debug("Skipping %s: gene count %d > max %d", accession, gene_count, max_gene_count)
                filter_stats['gene_count'] += 1
                continue

        # FILTER 3: Nucleotide completeness filter
        if nuc_completeness is not None:
            completeness_status = metadata.get("completeness")
            if completeness_status is None:
                logger.debug("Skipping %s: missing completeness metadata", accession)
                filter_stats['completeness'] += 1
                continue
                
            if completeness_status.lower() != nuc_completeness.lower():
                logger.debug("Skipping %s: completeness '%s' != required '%s'", 
                           accession, completeness_status, nuc_completeness)
                filter_stats['completeness'] += 1
                continue

        # FILTER 4: Host organism name filter
        if host is not None:
            # Convert host organism name to lowercase with underscores for comparison
            host_organism = "_".join(
                metadata.get("host", {}).get("organismName", "").split(" ")
            ).lower()
            
            if not host_organism:
                logger.debug("Skipping %s: missing host organism name", accession)
                filter_stats['host'] += 1
                continue
                
            if host_organism != host.lower():
                logger.debug("Skipping %s: host '%s' != required '%s'", 
                           accession, host_organism, host.lower())
                filter_stats['host'] += 1
                continue

        # FILTER 5: Host taxonomy ID filter
        if host_taxid is not None:
            host_lineage = metadata.get("host", {}).get("lineage", [])
            if not host_lineage:
                logger.debug("Skipping %s: missing host lineage metadata", accession)
                filter_stats['host_taxid'] += 1
                continue
                
            # Extract all taxonomy IDs from the host lineage
            host_lineage_taxids = {lineage["taxId"] for lineage in host_lineage}
            
            if host_taxid not in host_lineage_taxids:
                logger.debug("Skipping %s: host taxid %s not in lineage %s", 
                           accession, host_taxid, host_lineage_taxids)
                filter_stats['host_taxid'] += 1
                continue

        # FILTER 6: Lab passaging status filter
        if lab_passaged is True:
            # Only include sequences that have been lab-passaged
            from_lab = metadata.get("isLabHost")
            if not from_lab:
                logger.debug("Skipping %s: not lab-passaged (required)", accession)
                filter_stats['lab_passaged'] += 1
                continue

        if lab_passaged is False:
            # Only include sequences that have NOT been lab-passaged
            from_lab = metadata.get("isLabHost")
            if from_lab:
                logger.debug("Skipping %s: is lab-passaged (excluded)", accession)
                filter_stats['lab_passaged'] += 1
                continue

        # FILTER 7: Geographic region filter
        if geographic_region is not None:
            # Convert geographic region to lowercase with underscores for comparison
            location = "_".join(
                metadata.get("location", {}).get("geographicRegion", "").split(" ")
            ).lower()
            
            if not location:
                logger.debug("Skipping %s: missing geographic region", accession)
                filter_stats['geographic_region'] += 1
                continue
                
            if location != geographic_region.lower():
                logger.debug("Skipping %s: geographic region '%s' != required '%s'", 
                           accession, location, geographic_region.lower())
                filter_stats['geographic_region'] += 1
                continue

        # FILTER 8: Submitter country filter
        if submitter_country is not None:
            # Convert submitter country to lowercase with underscores for comparison
            submitter_country_value = "_".join(
                metadata.get("submitter", {}).get("country", "").split(" ")
            ).lower()
            
            if not submitter_country_value:
                logger.debug("Skipping %s: missing submitter country", accession)
                filter_stats['submitter_country'] += 1
                continue
                
            if submitter_country_value != submitter_country.lower():
                logger.debug("Skipping %s: submitter country '%s' != required '%s'", 
                           accession, submitter_country_value, submitter_country.lower())
                filter_stats['submitter_country'] += 1
                continue

        # FILTER 9: Collection date range filter
        if min_collection_date is not None or max_collection_date is not None:
            date_str = metadata.get("isolate", {}).get("collectionDate", "")
            
            # Parse the collection date
            date = parse_date(date_str)
            
            if date_str is None or date is None:
                logger.debug("Skipping %s: missing or invalid collection date '%s'", accession, date_str)
                filter_stats['collection_date'] += 1
                continue
                
            if min_collection_date and date < min_collection_date:
                logger.debug("Skipping %s: collection date %s < min %s", 
                           accession, date, min_collection_date)
                filter_stats['collection_date'] += 1
                continue
                
            if max_collection_date and date > max_collection_date:
                logger.debug("Skipping %s: collection date %s > max %s", 
                           accession, date, max_collection_date)
                filter_stats['collection_date'] += 1
                continue

        # FILTER 10: Annotation status filter (for annotated=False only)
        # Note: annotated=True filter is applied when API datasets call is made
        if annotated is False:
            annotated_value = metadata.get("isAnnotated")
            if annotated_value:
                logger.debug("Skipping %s: sequence is annotated (excluded)", accession)
                filter_stats['annotated'] += 1
                continue

        # FILTER 11: Source database filter
        if source_database is not None:
            source_db = metadata.get("sourceDatabase", "").lower()
            if not source_db:
                logger.debug("Skipping %s: missing source database", accession)
                filter_stats['source_database'] += 1
                continue
                
            if source_db != source_database.lower():
                logger.debug("Skipping %s: source database '%s' != required '%s'", 
                           accession, source_db, source_database.lower())
                filter_stats['source_database'] += 1
                continue

        # FILTER 12: Maximum release date filter
        # Note: minimum release date filter is applied at the API level
        if max_release_date is not None:
            release_date_str = metadata.get("releaseDate")
            
            if not release_date_str:
                logger.debug("Skipping %s: missing release date", accession)
                filter_stats['release_date'] += 1
                continue
                
            # Parse release date (remove time component if present)
            release_date_value = parse_date(release_date_str.split("T")[0])
            
            if release_date_value is None:
                logger.debug("Skipping %s: invalid release date '%s'", accession, release_date_str)
                filter_stats['release_date'] += 1
                continue
                
            if release_date_value > max_release_date:
                logger.debug("Skipping %s: release date %s > max %s", 
                           accession, release_date_value, max_release_date)
                filter_stats['release_date'] += 1
                continue

        # FILTER 13: Mature peptide count filters
        if min_mature_peptide_count is not None or max_mature_peptide_count is not None:
            mature_peptide_count = metadata.get("maturePeptideCount")
            
            if mature_peptide_count is None:
                logger.debug("Skipping %s: missing mature peptide count", accession)
                filter_stats['mature_peptide_count'] += 1
                continue
                
            if (min_mature_peptide_count is not None and 
                mature_peptide_count < min_mature_peptide_count):
                logger.debug("Skipping %s: mature peptide count %d < min %d", 
                           accession, mature_peptide_count, min_mature_peptide_count)
                filter_stats['mature_peptide_count'] += 1
                continue
                
            if (max_mature_peptide_count is not None and 
                mature_peptide_count > max_mature_peptide_count):
                logger.debug("Skipping %s: mature peptide count %d > max %d", 
                           accession, mature_peptide_count, max_mature_peptide_count)
                filter_stats['mature_peptide_count'] += 1
                continue

        # FILTER 14: Protein count filters
        if min_protein_count is not None or max_protein_count is not None:
            protein_count = metadata.get("proteinCount")
            
            if protein_count is None:
                logger.debug("Skipping %s: missing protein count", accession)
                filter_stats['protein_count'] += 1
                continue
                
            if min_protein_count is not None and protein_count < min_protein_count:
                logger.debug("Skipping %s: protein count %d < min %d", 
                           accession, protein_count, min_protein_count)
                filter_stats['protein_count'] += 1
                continue
                
            if max_protein_count is not None and protein_count > max_protein_count:
                logger.debug("Skipping %s: protein count %d > max %d", 
                           accession, protein_count, max_protein_count)
                filter_stats['protein_count'] += 1
                continue

        # FILTER 15: Ambiguous nucleotide character filter
        # Filter out sequences containing too many 'N' or 'n' characters (ambiguous nucleotides)
        if max_ambiguous_chars is not None:
            sequence_str = str(record.seq)
            # Count both uppercase and lowercase 'N' characters
            n_count = sequence_str.upper().count("N")

            if n_count > max_ambiguous_chars:
                logger.debug("Skipping %s: ambiguous chars %d > max %d", 
                           accession, n_count, max_ambiguous_chars)
                filter_stats['ambiguous_chars'] += 1
                continue

        # FILTER 16: Required proteins/genes/segments filter
        # Check if requested proteins/segments are present based on sequence header labels
        if has_proteins is not None:
            # Convert single protein to list for consistent processing
            if isinstance(has_proteins, str):
                has_proteins = [has_proteins]
                
            try:
                # Extract protein/segment information from the FASTA header
                if metadata.get("isolate", {}).get("name"):
                    # If isolate name is available, use everything after it in the header
                    prot_header = record.description.split(
                        metadata.get("isolate", {}).get("name")
                    )[-1]
                else:
                    # If sample name was not added to metadata,
                    # search the whole header for protein/segment names
                    prot_header = record.description
                    
                # Split header into segments for protein searching
                prot_parts = prot_header.split(";")

                # Check that all required proteins are present
                skip_outer_loop = False
                for protein in has_proteins:
                    # Create case-insensitive regex for each protein with flexible quote handling
                    regex = rf"(?i)\b['\",]?\(?{protein}\)?['\",]?\b"

                    if proteins_complete:
                        # Only keep sequences where proteins are marked as "complete"
                        if not any(
                            re.search(regex, part) and "complete" in part
                            for part in prot_parts
                        ):
                            logger.debug("Skipping %s: protein '%s' not marked complete", accession, protein)
                            skip_outer_loop = True
                            break
                    else:
                        # Just check for presence of the protein, regardless of completeness
                        if not any(re.search(regex, part) for part in prot_parts):
                            logger.debug("Skipping %s: protein '%s' not found", accession, protein)
                            skip_outer_loop = True
                            break
                            
                if skip_outer_loop:
                    filter_stats['has_proteins'] += 1
                    continue

            except (AttributeError, KeyError, IndexError) as e:
                logger.warning(
                    "The 'has_proteins' filter could not be applied to sequence %s due to the following error:\n%s", 
                    record.id, e
                )
                filter_stats['has_proteins'] += 1
                continue

        # If we reach this point, the sequence has passed all filters
        # Extract protein/segment information from the sequence header for output
        try:
            if metadata.get("isolate", {}).get("name"):
                # Use everything after the isolate name as protein/segment description
                prot_header = record.description.split(
                    metadata.get("isolate", {}).get("name")
                )[-1]
            else:
                # If sample name was not added to metadata,
                # use the whole header as protein/segment description
                prot_header = record.description
        except (AttributeError, KeyError, IndexError):
            # Handle cases where header parsing fails
            prot_header = pd.NA
            logger.debug("Could not extract protein header for %s", accession)

        # Add the sequence and its metadata to the results
        protein_headers.append(prot_header)
        filtered_sequences.append(record)
        filtered_metadata.append(metadata)
        
        logger.debug("Sequence %s passed all filters", accession)

    # Log comprehensive filtering statistics
    num_seqs = len(filtered_sequences)
    logger.info("Filtering complete:")
    logger.info("  Total sequences processed: %d", total_sequences)
    logger.info("  Sequences with metadata: %d", sequences_with_metadata)
    logger.info("  Sequences passing all filters: %d", num_seqs)
    
    # Log detailed filter statistics if any sequences were filtered out
    total_filtered = sum(filter_stats.values())
    if total_filtered > 0:
        logger.info("Filter statistics (sequences excluded):")
        for filter_name, count in filter_stats.items():
            if count > 0:
                logger.info("  %s: %d sequences", filter_name, count)
    
    # Validate results and return
    if num_seqs > 0:
        # Perform consistency checks
        if num_seqs != len(filtered_metadata):
            logger.warning(
                "Number of sequences (%d) and number of metadata entries (%d) do not match.", 
                num_seqs, len(filtered_metadata)
            )
        if num_seqs != len(protein_headers):
            logger.warning(
                "Number of sequences (%d) and number of protein headers (%d) do not match.", 
                num_seqs, len(protein_headers)
            )
            
        logger.debug(
            "Final counts - sequences: %d, metadata: %d, protein headers: %d", 
            num_seqs, len(filtered_metadata), len(protein_headers)
        )
        
        return filtered_sequences, filtered_metadata, protein_headers
    else:
        logger.warning("No sequences passed the provided filters.")
        return None, None, None


def save_metadata_to_csv(filtered_metadata, protein_headers, output_metadata_file):
    """
    Save filtered metadata to a CSV file with a specific column order.
    
    This function creates a comprehensive CSV file containing all relevant metadata
    for the filtered virus sequences. The output format is designed to be compatible
    with downstream analysis tools like Delphy.
    
    Args:
        filtered_metadata (list): List of metadata dictionaries for filtered sequences
        protein_headers (list): List of protein/segment information extracted from headers
        output_metadata_file (str): Path to the output CSV file
        
    Note:
        The column order is specifically designed to match requirements for
        phylogenetic analysis tools and provides a standardized format.
    """

    logger.info("Preparing metadata for CSV output...")
    logger.debug("Processing %d metadata records with %d protein headers", 
                len(filtered_metadata), len(protein_headers))

    # Define the column order for the output CSV
    # This order prioritizes the most commonly used fields and matches
    # the format expected by downstream analysis tools
    columns = [
        "accession",           # Primary identifier (lowercase for Delphy compatibility)
        "Organism Name",       # Virus species/strain name
        "GenBank/RefSeq",      # Source database (GenBank or RefSeq)
        "Submitters",          # Names of sequence submitters
        "Organization",        # Submitting organization/institution
        "Submitter Country",   # Country of submitting organization
        "Release date",        # Date when sequence was released to public databases
        "Isolate",            # Isolate/sample identifier
        "Virus Lineage",      # Taxonomic lineage of the virus
        "Length",             # Sequence length in base pairs
        "Nuc Completeness",   # Completeness status (complete/partial)
        "Proteins/Segments",  # Protein/segment information from FASTA headers
        "Geographic Region",  # Geographic region where sample was collected
        "Geographic Location",# Specific geographic location
        "Geo String",         # Full geographic information string
        "Host",               # Host organism name
        "Host Lineage",       # Taxonomic lineage of host organism
        "Lab Host",           # Whether sample was lab-passaged
        "Tissue/Specimen/Source", # Sample source/tissue type
        "Collection Date",    # Date when sample was collected
        "Sample Name",        # Sample identifier
        "Annotated",          # Whether sequence has annotation data
        "SRA Accessions",     # Associated SRA (sequencing) accessions
        "Bioprojects",        # Associated BioProject identifiers
        "Biosample",          # BioSample identifier
        "Protein count",      # Number of proteins annotated
        "Gene count",         # Number of genes annotated
    ]

    logger.debug("Using column order: %s", columns)

    # Prepare data for DataFrame creation
    data_for_df = []
    
    logger.info("Processing metadata records...")
    for i, metadata in enumerate(filtered_metadata):
        logger.debug("Processing metadata record %d/%d", i+1, len(filtered_metadata))
        
        # Extract and process geographic location information
        # The API returns nested geographic information that needs to be flattened
        location_info = metadata.get("location", {})
        
        # Filter out empty values and reverse the order for hierarchical display
        location_values = [v for v in location_info.values() if v and v != ""]
        location_values.reverse()
        
        # Create a colon-separated geographic string (e.g., "Africa:South Africa:Cape Town")
        geo_info = ":".join(location_values) if location_values else pd.NA
        
        # Extract geographic region (broadest level, e.g., "Africa", "Europe")
        try:
            geo_region = location_values[0] if location_values else pd.NA
        except IndexError:
            geo_region = pd.NA
            
        # Extract specific geographic location (e.g., "South Africa", "Germany")
        try:
            geo_loc = location_values[1] if len(location_values) > 1 else pd.NA
        except IndexError:
            geo_loc = pd.NA

        logger.debug("Geographic info for record %d: region=%s, location=%s, full=%s", 
                    i+1, geo_region, geo_loc, geo_info)

        # Build the row dictionary with all required columns
        # Use pd.NA for missing values to ensure proper CSV handling
        row = {
            # Primary identifiers
            "accession": metadata.get("accession", pd.NA),
            "Organism Name": metadata.get("virus", {}).get("organismName", pd.NA),
            
            # Database and submission information
            "GenBank/RefSeq": metadata.get("sourceDatabase", pd.NA),
            "Submitters": ", ".join(metadata.get("submitter", {}).get("names", [])),
            "Organization": metadata.get("submitter", {}).get("affiliation", pd.NA),
            "Submitter Country": metadata.get("submitter", {}).get("country", ""),
            "Release date": metadata.get("releaseDate", "").split("T")[0],  # Remove time component
            
            # Sample and isolate information
            "Isolate": metadata.get("isolate", {}).get("name", pd.NA),
            "Sample Name": metadata.get("isolate", {}).get("name", pd.NA),
            
            # Virus classification
            "Virus Lineage": metadata.get("virus", {}).get("lineage", []),
            
            # Sequence characteristics
            "Length": metadata.get("length", pd.NA),
            "Nuc Completeness": metadata.get("completeness", pd.NA),
            "Proteins/Segments": protein_headers[i] if i < len(protein_headers) else pd.NA,
            
            # Geographic information
            "Geographic Region": geo_region,
            "Geographic Location": geo_loc,
            "Geo String": geo_info,
            
            # Host information
            "Host": metadata.get("host", {}).get("organismName", pd.NA),
            "Host Lineage": metadata.get("host", {}).get("lineage", []),
            "Lab Host": metadata.get("labHost", pd.NA),
            
            # Sample source information
            "Tissue/Specimen/Source": metadata.get("isolate", {}).get("source", pd.NA),
            "Collection Date": metadata.get("isolate", {}).get("collectionDate", pd.NA),
            
            # Annotation and quality information
            "Annotated": metadata.get("isAnnotated", pd.NA),
            
            # Associated database records
            "SRA Accessions": metadata.get("sraAccessions", []),
            "Bioprojects": metadata.get("bioprojects", []),
            "Biosample": metadata.get("biosample", pd.NA),
            
            # Counts
            "Gene count": metadata.get("geneCount"),
            "Protein count": metadata.get("proteinCount"),
        }
        
        data_for_df.append(row)

    logger.info("Creating DataFrame with %d rows and %d columns", len(data_for_df), len(columns))

    # Create DataFrame with the specified column order
    df = pd.DataFrame(data_for_df, columns=columns)
    
    logger.debug("DataFrame shape: %s", df.shape)
    logger.debug("DataFrame columns: %s", list(df.columns))

    # Write DataFrame to CSV file
    try:
        df.to_csv(output_metadata_file, index=False)
        logger.info("Successfully saved metadata CSV to: %s", output_metadata_file)
        logger.debug("CSV file contains %d rows and %d columns", len(df), len(df.columns))
    except Exception as e:
        logger.error("Failed to save CSV file: %s", e)
        raise


def check_min_max(min_val, max_val, filtername, date=False):
    """
    Validate that minimum and maximum values are in the correct order.
    
    This helper function ensures that minimum values are not greater than
    maximum values for range-based filters. It handles both numeric and
    date-based comparisons.
    
    Args:
        min_val: Minimum value (can be numeric or date string)
        max_val: Maximum value (can be numeric or date string)
        filtername (str): Name of the filter for error reporting
        date (bool): Whether the values are dates that need parsing
        
    Raises:
        ValueError: If minimum value is greater than maximum value
        
    Example:
        check_min_max(100, 50, "sequence length")  # Raises ValueError
        check_min_max(100, 200, "sequence length")  # No error
    """
    # Only perform validation if both values are provided
    if min_val is not None and max_val is not None:
        logger.debug("Validating min/max values for %s: min=%s, max=%s", 
                    filtername, min_val, max_val)
        
        if date:
            # Parse date strings for comparison
            try:
                min_val = parse_date(min_val)
                max_val = parse_date(max_val)
                logger.debug("Parsed date values: min=%s, max=%s", min_val, max_val)
            except Exception as e:
                logger.error("Failed to parse dates for validation: %s", e)
                raise ValueError(f"Invalid date format in {filtername} filters") from e
        
        # Check if minimum is greater than maximum
        if min_val > max_val:
            error_msg = f"Min value ({min_val}) cannot be greater than max value ({max_val}) for {filtername}."
            logger.error("Validation failed: %s", error_msg)
            raise ValueError(error_msg)
        
        logger.debug("Min/max validation passed for %s", filtername)


def filter_metadata_only(
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
    submitter_country=None,
    min_collection_date=None,
    max_collection_date=None,
    annotated=None,
    source_database=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
    # Sequence-dependent filters are deferred
    max_ambiguous_chars=None,
    has_proteins=None,
    proteins_complete=False,
):
    """
    Filter metadata records based on metadata-only criteria.
    
    This function applies filters that can be evaluated using only metadata,
    allowing us to reduce the number of accessions before downloading sequences.
    Sequence-dependent filters (max_ambiguous_chars, has_proteins) are deferred
    to the post-download filtering step.
    
    Args:
        metadata_dict (dict): Dictionary mapping accession numbers to metadata
        (all other args): Same as filter_sequences function
        
    Returns:
        tuple: (filtered_accessions, filtered_metadata_list)
               - filtered_accessions: List of accession numbers that passed filters
               - filtered_metadata_list: List of metadata dictionaries for filtered accessions
    """
    
    logger.info("Starting metadata-only filtering process...")
    logger.debug("Applying metadata-only filters: seq_length(%s-%s), gene_count(%s-%s), "
                "completeness(%s), host(%s), host_taxid(%s), lab_passaged(%s), "
                "geo_region(%s), submitter_country(%s), collection_date(%s-%s), "
                "annotated(%s), source_db(%s), max_release_date(%s), "
                "peptide_count(%s-%s), protein_count(%s-%s)",
                min_seq_length, max_seq_length, min_gene_count, max_gene_count,
                nuc_completeness, host, host_taxid, lab_passaged, geographic_region,
                submitter_country, min_collection_date, max_collection_date,
                annotated, source_database, max_release_date, 
                min_mature_peptide_count, max_mature_peptide_count,
                min_protein_count, max_protein_count)
    
    # Convert date filters to datetime objects for proper comparison
    min_collection_date = (
        parse_date(min_collection_date, filtername="min_collection_date", verbose=True) 
        if min_collection_date else None
    )
    max_collection_date = (
        parse_date(max_collection_date, filtername="max_collection_date", verbose=True) 
        if max_collection_date else None
    )
    max_release_date = (
        parse_date(max_release_date, filtername="max_release_date", verbose=True) 
        if max_release_date else None
    )
    
    if min_collection_date:
        logger.debug("Parsed min_collection_date: %s", min_collection_date)
    if max_collection_date:
        logger.debug("Parsed max_collection_date: %s", max_collection_date)
    if max_release_date:
        logger.debug("Parsed max_release_date: %s", max_release_date)

    # Initialize lists to store filtered results
    filtered_accessions = []
    filtered_metadata_list = []
    
    # Counters for logging filter statistics
    total_sequences = len(metadata_dict)
    filter_stats = {
        'seq_length': 0,
        'gene_count': 0,
        'completeness': 0,
        'host': 0,
        'host_taxid': 0,
        'lab_passaged': 0,
        'geographic_region': 0,
        'submitter_country': 0,
        'collection_date': 0,
        'annotated': 0,
        'source_database': 0,
        'release_date': 0,
        'mature_peptide_count': 0,
        'protein_count': 0,
    }

    logger.info("Processing %d metadata records...", total_sequences)
    
    for accession, metadata in metadata_dict.items():
        logger.debug("Processing metadata for: %s", accession)
        
        # Apply filters sequentially - each filter can exclude the record
        # If any filter fails, we continue to the next record
        
        # FILTER 1: Sequence length filters
        if min_seq_length is not None or max_seq_length is not None:
            sequence_length = metadata.get("length")
            if sequence_length is None:
                logger.debug("Skipping %s: missing length metadata", accession)
                filter_stats['seq_length'] += 1
                continue
                
            if min_seq_length is not None and sequence_length < min_seq_length:
                logger.debug("Skipping %s: length %d < min %d", accession, sequence_length, min_seq_length)
                filter_stats['seq_length'] += 1
                continue
                
            if max_seq_length is not None and sequence_length > max_seq_length:
                logger.debug("Skipping %s: length %d > max %d", accession, sequence_length, max_seq_length)
                filter_stats['seq_length'] += 1
                continue

        # FILTER 2: Gene count filters
        if min_gene_count is not None or max_gene_count is not None:
            gene_count = metadata.get("geneCount")
            if gene_count is None:
                logger.debug("Skipping %s: missing gene count metadata", accession)
                filter_stats['gene_count'] += 1
                continue
                
            if min_gene_count is not None and gene_count < min_gene_count:
                logger.debug("Skipping %s: gene count %d < min %d", accession, gene_count, min_gene_count)
                filter_stats['gene_count'] += 1
                continue
                
            if max_gene_count is not None and gene_count > max_gene_count:
                logger.debug("Skipping %s: gene count %d > max %d", accession, gene_count, max_gene_count)
                filter_stats['gene_count'] += 1
                continue

        # FILTER 3: Nucleotide completeness filter
        if nuc_completeness is not None:
            completeness_status = metadata.get("completeness")
            if completeness_status is None:
                logger.debug("Skipping %s: missing completeness metadata", accession)
                filter_stats['completeness'] += 1
                continue
                
            if completeness_status.lower() != nuc_completeness.lower():
                logger.debug("Skipping %s: completeness '%s' != required '%s'", 
                           accession, completeness_status, nuc_completeness)
                filter_stats['completeness'] += 1
                continue

        # FILTER 4: Host organism name filter
        if host is not None:
            host_organism = "_".join(
                metadata.get("host", {}).get("organismName", "").split(" ")
            ).lower()
            
            if not host_organism:
                logger.debug("Skipping %s: missing host organism name", accession)
                filter_stats['host'] += 1
                continue
                
            if host_organism != host.lower():
                logger.debug("Skipping %s: host '%s' != required '%s'", 
                           accession, host_organism, host.lower())
                filter_stats['host'] += 1
                continue

        # FILTER 5: Host taxonomy ID filter
        if host_taxid is not None:
            host_lineage = metadata.get("host", {}).get("lineage", [])
            if not host_lineage:
                logger.debug("Skipping %s: missing host lineage metadata", accession)
                filter_stats['host_taxid'] += 1
                continue
                
            host_lineage_taxids = {lineage["taxId"] for lineage in host_lineage}
            
            if host_taxid not in host_lineage_taxids:
                logger.debug("Skipping %s: host taxid %s not in lineage %s", 
                           accession, host_taxid, host_lineage_taxids)
                filter_stats['host_taxid'] += 1
                continue

        # FILTER 6: Lab passaging status filter
        if lab_passaged is True:
            from_lab = metadata.get("isLabHost")
            if not from_lab:
                logger.debug("Skipping %s: not lab-passaged (required)", accession)
                filter_stats['lab_passaged'] += 1
                continue

        if lab_passaged is False:
            from_lab = metadata.get("isLabHost")
            if from_lab:
                logger.debug("Skipping %s: is lab-passaged (excluded)", accession)
                filter_stats['lab_passaged'] += 1
                continue

        # FILTER 7: Geographic region filter
        if geographic_region is not None:
            location = "_".join(
                metadata.get("location", {}).get("geographicRegion", "").split(" ")
            ).lower()
            
            if not location:
                logger.debug("Skipping %s: missing geographic region", accession)
                filter_stats['geographic_region'] += 1
                continue
                
            if location != geographic_region.lower():
                logger.debug("Skipping %s: geographic region '%s' != required '%s'", 
                           accession, location, geographic_region.lower())
                filter_stats['geographic_region'] += 1
                continue

        # FILTER 8: Submitter country filter
        if submitter_country is not None:
            submitter_country_value = "_".join(
                metadata.get("submitter", {}).get("country", "").split(" ")
            ).lower()
            
            if not submitter_country_value:
                logger.debug("Skipping %s: missing submitter country", accession)
                filter_stats['submitter_country'] += 1
                continue
                
            if submitter_country_value != submitter_country.lower():
                logger.debug("Skipping %s: submitter country '%s' != required '%s'", 
                           accession, submitter_country_value, submitter_country.lower())
                filter_stats['submitter_country'] += 1
                continue

        # FILTER 9: Collection date range filter
        if min_collection_date is not None or max_collection_date is not None:
            date_str = metadata.get("isolate", {}).get("collectionDate", "")
            
            date = parse_date(date_str)
            
            if date_str is None or date is None:
                logger.debug("Skipping %s: missing or invalid collection date '%s'", accession, date_str)
                filter_stats['collection_date'] += 1
                continue
                
            if min_collection_date and date < min_collection_date:
                logger.debug("Skipping %s: collection date %s < min %s", 
                           accession, date, min_collection_date)
                filter_stats['collection_date'] += 1
                continue
                
            if max_collection_date and date > max_collection_date:
                logger.debug("Skipping %s: collection date %s > max %s", 
                           accession, date, max_collection_date)
                filter_stats['collection_date'] += 1
                continue

        # FILTER 10: Annotation status filter (for annotated=False only)
        if annotated is False:
            annotated_value = metadata.get("isAnnotated")
            if annotated_value:
                logger.debug("Skipping %s: sequence is annotated (excluded)", accession)
                filter_stats['annotated'] += 1
                continue

        # FILTER 11: Source database filter
        if source_database is not None:
            source_db = metadata.get("sourceDatabase", "").lower()
            if not source_db:
                logger.debug("Skipping %s: missing source database", accession)
                filter_stats['source_database'] += 1
                continue
                
            if source_db != source_database.lower():
                logger.debug("Skipping %s: source database '%s' != required '%s'", 
                           accession, source_db, source_database.lower())
                filter_stats['source_database'] += 1
                continue

        # FILTER 12: Maximum release date filter
        if max_release_date is not None:
            release_date_str = metadata.get("releaseDate")
            
            if not release_date_str:
                logger.debug("Skipping %s: missing release date", accession)
                filter_stats['release_date'] += 1
                continue
                
            release_date_value = parse_date(release_date_str.split("T")[0])
            
            if release_date_value is None:
                logger.debug("Skipping %s: invalid release date '%s'", accession, release_date_str)
                filter_stats['release_date'] += 1
                continue
                
            if release_date_value > max_release_date:
                logger.debug("Skipping %s: release date %s > max %s", 
                           accession, release_date_value, max_release_date)
                filter_stats['release_date'] += 1
                continue

        # FILTER 13: Mature peptide count filters
        if min_mature_peptide_count is not None or max_mature_peptide_count is not None:
            mature_peptide_count = metadata.get("maturePeptideCount")
            
            if mature_peptide_count is None:
                logger.debug("Skipping %s: missing mature peptide count", accession)
                filter_stats['mature_peptide_count'] += 1
                continue
                
            if (min_mature_peptide_count is not None and 
                mature_peptide_count < min_mature_peptide_count):
                logger.debug("Skipping %s: mature peptide count %d < min %d", 
                           accession, mature_peptide_count, min_mature_peptide_count)
                filter_stats['mature_peptide_count'] += 1
                continue
                
            if (max_mature_peptide_count is not None and 
                mature_peptide_count > max_mature_peptide_count):
                logger.debug("Skipping %s: mature peptide count %d > max %d", 
                           accession, mature_peptide_count, max_mature_peptide_count)
                filter_stats['mature_peptide_count'] += 1
                continue

        # FILTER 14: Protein count filters
        if min_protein_count is not None or max_protein_count is not None:
            protein_count = metadata.get("proteinCount")
            
            if protein_count is None:
                logger.debug("Skipping %s: missing protein count", accession)
                filter_stats['protein_count'] += 1
                continue
                
            if min_protein_count is not None and protein_count < min_protein_count:
                logger.debug("Skipping %s: protein count %d < min %d", 
                           accession, protein_count, min_protein_count)
                filter_stats['protein_count'] += 1
                continue
                
            if max_protein_count is not None and protein_count > max_protein_count:
                logger.debug("Skipping %s: protein count %d > max %d", 
                           accession, protein_count, max_protein_count)
                filter_stats['protein_count'] += 1
                continue

        # If we reach this point, the metadata record has passed all filters
        filtered_accessions.append(accession)
        filtered_metadata_list.append(metadata)
        
        logger.debug("Metadata %s passed all filters", accession)

    # Log comprehensive filtering statistics
    num_filtered = len(filtered_accessions)
    logger.info("Metadata-only filtering complete:")
    logger.info("  Total metadata records: %d", total_sequences)
    logger.info("  Records passing filters: %d", num_filtered)
    
    # Log detailed filter statistics if any records were filtered out
    total_filtered = sum(filter_stats.values())
    if total_filtered > 0:
        logger.info("Filter statistics (records excluded):")
        for filter_name, count in filter_stats.items():
            if count > 0:
                logger.info("  %s: %d records", filter_name, count)
    
    if max_ambiguous_chars is not None or has_proteins is not None:
        logger.info("Note: Sequence-dependent filters (max_ambiguous_chars, has_proteins) will be applied after download")
    
    return filtered_accessions, filtered_metadata_list


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
    has_proteins=None,
    proteins_complete=False,
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
):
    """
    Download a virus genome dataset from the NCBI Virus database (https://www.ncbi.nlm.nih.gov/labs/virus/).

    This is the main function that orchestrates the entire virus data retrieval process,
    now optimized to download sequences only after all metadata-based filtering:
    1) Validate inputs
    2) Fetch and save API metadata (server-side filters applied)
    3) Apply metadata-only filters locally
    4) Download sequences for the filtered accession list only
    5) Apply sequence-dependent filters and save outputs
    """
    logger.info("Starting NCBI virus data retrieval process...")
    logger.info("Query parameters: virus='%s', accession=%s, outfolder='%s'", 
                virus, accession, outfolder)
    logger.debug("Applied filters: host=%s, seq_length=(%s-%s), gene_count=(%s-%s), "
                "completeness=%s, annotated=%s, lab_passaged=%s, geo_region=%s, "
                "geo_location=%s, submitter_country=%s, collection_date=(%s-%s), "
                "source_db=%s, release_date=(%s-%s), protein_count=(%s-%s), "
                "peptide_count=(%s-%s), max_ambiguous=%s, has_proteins=%s, proteins_complete=%s",
                host, min_seq_length, max_seq_length, min_gene_count, max_gene_count,
                nuc_completeness, annotated, lab_passaged, geographic_region,
                geographic_location, submitter_country, min_collection_date, max_collection_date,
                source_database, min_release_date, max_release_date, min_protein_count,
                max_protein_count, min_mature_peptide_count, max_mature_peptide_count,
                max_ambiguous_chars, has_proteins, proteins_complete)

    # SECTION 1: INPUT VALIDATION
    # Validate and normalize input arguments before proceeding
    logger.info("Validating input arguments...")
    
    # Validate nucleotide completeness argument
    if nuc_completeness is not None:
        nuc_completeness = nuc_completeness.lower()  # Normalize to lowercase
        if nuc_completeness not in ["partial", "complete"]:
            raise ValueError(
                "Argument 'nuc_completeness' must be 'partial', 'complete', or None."
            )
        logger.debug("Nucleotide completeness filter set to: %s", nuc_completeness)
    
    # Validate boolean arguments with proper type checking
    if annotated is not None and not isinstance(annotated, bool):
        raise TypeError(
            "Argument 'annotated' must be a boolean (True or False) or None."
        )
        
    if lab_passaged is not None and not isinstance(lab_passaged, bool):
        raise TypeError(
            "Argument 'lab_passaged' must be a boolean (True or False) or None."
        )
        
    if proteins_complete is not None and not isinstance(proteins_complete, bool):
        raise TypeError(
            "Argument 'proteins_complete' must be a boolean (True or False)."
        )
    
    # Convert integer virus identifiers to strings for API compatibility
    if isinstance(virus, int):
        virus = str(virus)
        logger.debug("Converted integer virus ID to string: %s", virus)

    # Validate min/max argument pairs to ensure logical consistency
    logger.debug("Validating min/max argument pairs...")
    check_min_max(
        min_seq_length,
        max_seq_length,
        "sequence length (arguments: min_seq_length and max_seq_length)",
    )
    check_min_max(
        min_gene_count,
        max_gene_count,
        "gene count (arguments: min_gene_count and max_gene_count)",
    )
    check_min_max(
        min_mature_peptide_count,
        max_mature_peptide_count,
        "mature peptide count (arguments: min_mature_peptide_count and max_mature_peptide_count)",
    )
    check_min_max(
        min_protein_count,
        max_protein_count,
        "protein count (arguments: min_protein_count and max_protein_count)",
    )
    check_min_max(
        min_release_date,
        max_release_date,
        "release data (arguments: min_release_date and max_release_date)",
        date=True,  # Enable date parsing for comparison
    )

    logger.info("Input validation completed successfully")

    # SECTION 2: OUTPUT DIRECTORY SETUP
    # Create and prepare output directory structure
    if outfolder is None:
        outfolder = os.getcwd()
        logger.info("No output folder specified, using current working directory: %s", outfolder)
    else:
        logger.info("Using specified output folder: %s", outfolder)
    
    # Ensure output folder exists
    os.makedirs(outfolder, exist_ok=True)
    logger.debug("Output folder ready: %s", outfolder)
    
    # Create temporary directory for intermediate processing
    # This will be cleaned up at the end regardless of success or failure
    temp_dir = os.path.join(outfolder, f"tmp_{UUID}")
    os.makedirs(temp_dir, exist_ok=True)
    logger.debug("Created temporary processing directory: %s", temp_dir)
    
    try:
        # SECTION 3: METADATA RETRIEVAL
        logger.info("=" * 60)
        logger.info("STEP 1: Fetching virus metadata from NCBI API")
        logger.info("=" * 60)
        logger.debug("Applying server-side filters: host=%s, geo_location=%s, annotated=%s, complete_only=%s, min_release_date=%s",
                     host, geographic_location, annotated, False, min_release_date)

        # Note: The API only supports annotated=True filter (annotated_only), not annotated=False
        # For annotated=False, we get all records and filter locally
        # Note: We don't apply complete_only at API level as it's too restrictive - filter locally instead
        api_annotated_filter = annotated if annotated is True else None
        
        api_reports = fetch_virus_metadata(
            virus,
            accession=accession,
            host=host,
            geographic_location=geographic_location,
            annotated=api_annotated_filter,  # Only pass True to API, not False
            complete_only=False,  # Don't filter at API level - too restrictive
            min_release_date=min_release_date,
        )

        if not api_reports:
            logger.warning("No virus records found matching the specified criteria.")
            logger.info("Consider relaxing your filter criteria or checking your virus identifier.")
            return

        logger.info("Successfully retrieved %d virus records from API", len(api_reports))

        # Convert API metadata to internal format
        logger.debug("Converting API metadata to internal format...")
        metadata_dict = load_metadata_from_api_reports(api_reports)
        logger.info("Processed metadata for %d sequences", len(metadata_dict))

        # Save the raw API metadata (server-side filtered) before local filtering
        virus_clean = virus.replace(' ', '_').replace('/', '_')
        output_api_metadata_jsonl = os.path.join(outfolder, f"{virus_clean}_api_metadata.jsonl")
        logger.debug("Writing API metadata to JSONL file: %s", output_api_metadata_jsonl)
        try:
            with open(output_api_metadata_jsonl, "w", encoding="utf-8") as f:
                for md in metadata_dict.values():
                    f.write(json.dumps(md) + "\n")
            logger.info("✓ Saved API metadata JSONL: %s", output_api_metadata_jsonl)
        except Exception as e:
            logger.warning("Failed to save API metadata JSONL: %s", e)

        # SECTION 4: METADATA-ONLY FILTERING
        logger.info("=" * 60)
        logger.info("STEP 2: Applying metadata-only filters")
        logger.info("=" * 60)

        filters = {
            "min_seq_length": min_seq_length,
            "max_seq_length": max_seq_length,
            "min_gene_count": min_gene_count,
            "max_gene_count": max_gene_count,
            # Apply completeness filter locally since API filter is too restrictive
            "nuc_completeness": nuc_completeness,
            # host already applied server-side when provided, so skip if it was applied at API level
            "host": None if host else None,
            "host_taxid": host_taxid,
            "lab_passaged": lab_passaged,
            "geographic_region": geographic_region,
            "submitter_country": submitter_country,
            "min_collection_date": min_collection_date,
            "max_collection_date": max_collection_date,
            # Only apply annotated filter locally if it wasn't applied at API level
            "annotated": annotated if annotated is False else None,  # API can't filter False, only True
            "source_database": source_database,
            "max_release_date": max_release_date,
            "min_mature_peptide_count": min_mature_peptide_count,
            "max_mature_peptide_count": max_mature_peptide_count,
            "min_protein_count": min_protein_count,
            "max_protein_count": max_protein_count,
            # deferred filters are still passed (function will defer them)
            "max_ambiguous_chars": max_ambiguous_chars,
            "has_proteins": has_proteins,
            "proteins_complete": proteins_complete,
        }

        filtered_accessions, filtered_metadata = filter_metadata_only(
            metadata_dict, **filters
        )

        if not filtered_accessions:
            logger.warning("No sequences passed metadata-only filters. Skipping sequence download.")
            return

        # Prepare output file paths
        output_fasta_file = os.path.join(outfolder, f"{virus_clean}_sequences.fasta")
        output_metadata_csv = os.path.join(outfolder, f"{virus_clean}_metadata.csv")
        output_metadata_jsonl = os.path.join(outfolder, f"{virus_clean}_metadata.jsonl")

        # Save filtered metadata immediately (before sequence-dependent fields)
        logger.debug("Writing filtered metadata (pre-sequence) to JSONL: %s", output_metadata_jsonl)
        try:
            with open(output_metadata_jsonl, "w", encoding="utf-8") as f:
                for md in filtered_metadata:
                    f.write(json.dumps(md) + "\n")
            logger.info("✓ Saved filtered metadata JSONL: %s", output_metadata_jsonl)
        except Exception as e:
            logger.warning("Failed to save filtered metadata JSONL: %s", e)

        # SECTION 5: DOWNLOAD SEQUENCES FOR FILTERED ACCESSIONS ONLY
        logger.info("=" * 60)
        logger.info("STEP 3: Downloading sequences for filtered accessions")
        logger.info("=" * 60)

        fna_file = download_sequences_by_accessions(filtered_accessions, outdir=temp_dir)
        if not os.path.exists(fna_file):
            raise RuntimeError(f"Download failed: FASTA file not found at {fna_file}")
        logger.info("Downloaded FASTA file: %s (%.2f MB)", fna_file, os.path.getsize(fna_file) / 1024 / 1024)

        # SECTION 6: SEQUENCE-DEPENDENT FILTERING AND SAVING
        logger.info("=" * 60)
        logger.info("STEP 4: Applying sequence-dependent filters and saving results")
        logger.info("=" * 60)

        # Restrict metadata to filtered accessions only
        filtered_metadata_dict = {acc: metadata_dict[acc] for acc in filtered_accessions}

        filtered_sequences, filtered_metadata_final, protein_headers = filter_sequences(
            fna_file,
            filtered_metadata_dict,
            **filters,
        )

        if filtered_sequences:
            logger.info("Saving %d filtered sequences and their metadata...", len(filtered_sequences))

            # Save FASTA
            FastaIO.write(filtered_sequences, output_fasta_file, "fasta")
            if os.path.exists(output_fasta_file):
                logger.info("✓ FASTA file saved: %s (%.2f MB)", output_fasta_file, os.path.getsize(output_fasta_file) / 1024 / 1024)
            else:
                logger.error("✗ Failed to create FASTA file: %s", output_fasta_file)

            # Overwrite JSONL with final filtered metadata (includes only sequences that passed all filters)
            try:
                with open(output_metadata_jsonl, "w", encoding="utf-8") as file:
                    for metadata in filtered_metadata_final:
                        file.write(json.dumps(metadata) + "\n")
                logger.info("✓ JSONL metadata file saved: %s (%.2f MB)", output_metadata_jsonl, os.path.getsize(output_metadata_jsonl) / 1024 / 1024)
            except Exception as e:
                logger.error("Failed to save JSONL metadata file: %s", e)
                raise

            # CSV
            try:
                save_metadata_to_csv(filtered_metadata_final, protein_headers, output_metadata_csv)
                if os.path.exists(output_metadata_csv):
                    logger.info("✓ CSV metadata file saved: %s (%.2f MB)", output_metadata_csv, os.path.getsize(output_metadata_csv) / 1024 / 1024)
                else:
                    logger.error("✗ Failed to create CSV file: %s", output_metadata_csv)
            except Exception as e:
                logger.error("Failed to save CSV metadata file: %s", e)
                raise

            # SECTION 10: FINAL SUMMARY
            # Provide comprehensive summary of the results
            logger.info("=" * 60)
            logger.info("PROCESS COMPLETED SUCCESSFULLY")
            logger.info("=" * 60)
            logger.info("Results summary:")
            logger.info("  Total sequences (API metadata): %d", len(metadata_dict))
            logger.info("  After metadata-only filtering: %d", len(filtered_accessions))
            logger.info("  After all filters (final): %d", len(filtered_sequences))
            logger.info("")
            logger.info("Output files saved to: %s", outfolder)
            logger.info("  📄 Sequences (FASTA): %s", os.path.basename(output_fasta_file))
            logger.info("  📊 Metadata (CSV):    %s", os.path.basename(output_metadata_csv))
            logger.info("  🔧 Metadata (JSONL):  %s", os.path.basename(output_metadata_jsonl))
            logger.info("=" * 60)
        else:
            logger.warning("=" * 60)
            logger.warning("NO SEQUENCES PASSED THE FILTERS")
            logger.warning("=" * 60)
            logger.warning("No sequences met all the specified filter criteria.")
            logger.warning("Consider:")
            logger.warning("  - Relaxing some filter parameters")
            logger.warning("  - Checking filter values for typos or incorrect formats")
            logger.warning("  - Trying a broader virus query term")
            logger.warning("  - Removing some of the more restrictive filters")
            logger.warning("=" * 60)

    except Exception as e:
        # Handle any unexpected errors during processing
        logger.error("An error occurred during virus data processing: %s", e)
        logger.error("Error type: %s", type(e).__name__)
        if logger.getEffectiveLevel() <= logging.DEBUG:
            import traceback
            logger.debug("Full traceback:\n%s", traceback.format_exc())
        raise
        
    finally:
        # SECTION 11: CLEANUP
        # Always clean up temporary files, regardless of success or failure
        logger.debug("Performing cleanup...")
        if os.path.exists(temp_dir):
            try:
                import shutil
                shutil.rmtree(temp_dir)
                logger.debug("✓ Cleaned up temporary directory: %s", temp_dir)
            except Exception as e:
                logger.warning("Failed to clean up temporary directory %s: %s", temp_dir, e)
        
        logger.info("NCBI virus data retrieval process completed.")