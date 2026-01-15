import os
import re
import json
import sys           # For accessing command line arguments
import time          # For adding delays between requests
import logging       # For logging level checks
import shutil        # For directory operations  
import subprocess    # For executing external commands
import traceback     # For error traceback logging
import platform      # For OS detection
import stat          # For file permission constants
import pandas as pd  # For data manipulation and CSV output
import requests      # For HTTP requests to NCBI API
import zipfile       # For extracting downloaded ZIP files
from datetime import datetime  # For date handling
from dateutil import parser  # For flexible date parsing
import xml.etree.ElementTree as ET # For XML parsing
import http.client
import urllib3
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

# Internal imports for logging, unique ID generation, and FASTA parsing
from .utils import set_up_logger, FastaIO
from .constants import NCBI_API_BASE, NCBI_EUTILS_BASE
from .compile import PACKAGE_PATH

# =============================================================================
# MODULE-LEVEL CONSTANTS
# =============================================================================

# API and Network Configuration
# -----------------------------------------------------------------------------
API_PAGE_SIZE = 1000  # Maximum records per API request (NCBI limit)
API_REQUEST_TIMEOUT = 30  # Timeout in seconds for API requests
API_MAX_RETRIES = 3  # Maximum retry attempts for failed API requests
API_INITIAL_RETRY_DELAY = 2.0  # Initial delay in seconds between retries
API_RETRY_BACKOFF_MULTIPLIER = 1.5  # Multiplier for exponential backoff

# E-utilities Configuration
EUTILS_TIMEOUT = 300  # Timeout in seconds for E-utilities requests
EUTILS_DEFAULT_BATCH_SIZE = 200  # Default batch size for E-utilities requests
EUTILS_INTER_BATCH_DELAY = 0.5  # Delay in seconds between batch requests
EUTILS_MIN_BATCH_SIZE_FOR_SPLIT = 50  # Minimum batch size before giving up on splitting

# GenBank Configuration
GENBANK_DEFAULT_BATCH_SIZE = 200  # Default batch size for GenBank requests
GENBANK_INTER_BATCH_DELAY = 0.5  # Delay in seconds between GenBank batch requests
GENBANK_MAX_BATCH_SIZE_WARNING = 500  # Warn user if batch size exceeds this
GENBANK_RETRY_ATTEMPTS = 5  # Number of retry attempts for GenBank requests
GENBANK_XML_CHUNK_SIZE = 10000  # Rows to process before writing to CSV

# Subprocess and Download Configuration
# -----------------------------------------------------------------------------
SUBPROCESS_VERSION_TIMEOUT = 5  # Timeout for version check commands
DOWNLOAD_OVERALL_TIMEOUT = 1800  # Maximum total download time (30 minutes)
DOWNLOAD_PROGRESS_TIMEOUT = 300  # Maximum time without progress (5 minutes)
DOWNLOAD_PROGRESS_CHECK_INTERVAL = 0.1  # Interval for checking download progress

# Chunked Download Configuration (for large datasets)
CHUNKED_DOWNLOAD_START_YEAR = 1970  # Default start year for chunked downloads
CHUNKED_DOWNLOAD_INTER_CHUNK_DELAY = 0.5  # Delay between yearly chunks

# NCBI Taxonomy IDs
# -----------------------------------------------------------------------------
NCBI_ALL_VIRUSES_TAXID = "10239"  # Taxonomy ID for all viruses
NCBI_SARS_COV2_TAXID = "2697049"  # Taxonomy ID for SARS-CoV-2
NCBI_ALPHAINFLUENZA_GENUS_TAXID = "197911"  # Alphainfluenza genus
NCBI_ALPHAINFLUENZA_SPECIES_TAXID = "2955291"  # Alphainfluenzavirus influenzae species
NCBI_INFLUENZA_A_TAXID = "11320"  # Influenza A virus

# Virus Detection Identifiers
# -----------------------------------------------------------------------------
SARS_COV2_IDENTIFIERS = {
    'sarscov2', 'sars2', '2697049', 'sarscov', 
    'severeacuterespiratorysyndromecoronavirus2',
    'covid19', 'covid', 'coronavirusdisease', 'ncov', 'hcov19'
}

ALPHAINFLUENZA_IDENTIFIERS = {
    'alphainfluenza', 'alphainfluenzavirus', 'alphainfluenzavirusinfluenzae',
    'influenzaavirus', 'influenzaa', 'flua',
    '197911',  # Alphainfluenza genus
    '2955291',  # Alphainfluenzavirus influenzae species
    '11320'  # Influenza A virus
}

# Default taxon for Alphainfluenza downloads (most comprehensive cached data)
ALPHAINFLUENZA_DEFAULT_TAXON = "Alphainfluenzavirus influenzae"

# Progress Indicator Keywords (for subprocess monitoring)
# -----------------------------------------------------------------------------
PROGRESS_INDICATORS = ['%', '=', 'downloading', 'fetching', 'MB', 'GB', 'bytes']

# Protein/Gene Keywords for Header Parsing
# -----------------------------------------------------------------------------
PROTEIN_KEYWORDS = [
    'hemagglutinin', 'neuraminidase', 'polymerase', 'nucleoprotein',
    'matrix protein', 'nonstructural protein', 'ns1', 'ns2',
    'spike', 'envelope', 'membrane', 'nucleocapsid',
    'orf', 'nsp', 'pp1a', 'pp1ab',
    'segment 1', 'segment 2', 'segment 3', 'segment 4',
    'segment 5', 'segment 6', 'segment 7', 'segment 8',
]

# Date Parsing Configuration
# -----------------------------------------------------------------------------
DATE_PARSE_DEFAULT_YEAR = 1000  # Default year for incomplete date strings

# HTTP Retry Configuration
# -----------------------------------------------------------------------------
HTTP_RETRY_STATUS_CODES = [429, 500, 502, 503, 504]  # Status codes to retry on
HTTP_MAX_LOCAL_RETRIES = 3  # Maximum local retry attempts
HTTP_INITIAL_BACKOFF = 1.0  # Initial backoff in seconds

# File Size Display Divisor
# -----------------------------------------------------------------------------
BYTES_PER_MB = 1024 * 1024  # Bytes in a megabyte for file size display

# =============================================================================
# End of Constants
# =============================================================================

# Set up logger for this module
logger = set_up_logger()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
random_suffix = os.urandom(4).hex() # random suffix for naming uniqueness

# Path to precompiled datasets binary
if platform.system() == "Windows":
    PRECOMPILED_DATASETS_PATH = os.path.join(
        PACKAGE_PATH, "bins", "Windows", "datasets.exe"
    )
else:
    PRECOMPILED_DATASETS_PATH = os.path.join(
        PACKAGE_PATH, "bins", platform.system(), "datasets"
    )

# Cache for the datasets path to avoid repeated checks
_datasets_path_cache = None


def _validate_datasets_binary(path):
    """
    Validate that a datasets binary exists and is functional.
    
    Args:
        path (str): Path to the datasets binary to validate.
        
    Returns:
        bool: True if the binary exists and runs successfully, False otherwise.
    """
    if not path:
        return False
    
    # Check if the file exists (for bundled binary) or is in PATH (for system binary)
    if not os.path.isfile(path) and not shutil.which(path):
        return False
    
    # Verify the binary actually works
    try:
        result = subprocess.run(
            [path, "--version"],
            capture_output=True,
            text=True,
            timeout=SUBPROCESS_VERSION_TIMEOUT,
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, OSError, FileNotFoundError):
        return False


def _clear_datasets_cache():
    """
    Clear the cached datasets path, forcing re-detection on next call.
    
    This is useful when the environment changes (e.g., user installs/uninstalls
    the datasets CLI) or when the cached binary becomes unavailable.
    """
    global _datasets_path_cache
    _datasets_path_cache = None
    logger.debug("Cleared datasets path cache")


def _get_datasets_path():
    """
    Get the path to the NCBI datasets CLI binary.

    This helper first checks if datasets is available in the system PATH.
    If found, it uses the system-installed version. Otherwise, it falls back
    to the precompiled binary bundled with gget.

    The result is cached after the first successful call. The cache is automatically
    invalidated if the cached binary becomes unavailable (e.g., deleted or environment
    changed), triggering re-detection.

    Returns:
        str: Path to the datasets binary ("datasets" for system PATH, or full path for bundled).

    Raises:
        RuntimeError: If no working datasets binary is available.
    """
    global _datasets_path_cache
    
    # If we have a cached path, validate it's still working
    if _datasets_path_cache is not None:
        if _validate_datasets_binary(_datasets_path_cache):
            return _datasets_path_cache
        else:
            # Cached binary is no longer valid, clear cache and re-detect
            logger.warning(
                "⚠️ Previously cached datasets binary at '%s' is no longer available. "
                "Re-detecting...", _datasets_path_cache
            )
            _clear_datasets_cache()
    
    # First, check if datasets is available in the system PATH
    datasets_path = shutil.which("datasets")
    if datasets_path:
        try:
            result = subprocess.run(
                [datasets_path, "--version"],
                capture_output=True,
                text=True,
                timeout=SUBPROCESS_VERSION_TIMEOUT,
            )
            if result.returncode == 0:
                logger.info(
                    "✅ Using system-installed NCBI datasets CLI: %s", result.stdout.strip()
                )
                _datasets_path_cache = datasets_path
                return datasets_path
        except (subprocess.TimeoutExpired, OSError):
            pass  # System binary didn't work, try bundled
    
    # Fall back to the bundled binary
    datasets_path = PRECOMPILED_DATASETS_PATH
    
    # Check if the precompiled binary exists
    if not os.path.isfile(datasets_path):
        raise RuntimeError(
            f"NCBI datasets binary not found at {datasets_path}. "
            "Please reinstall gget to restore the bundled datasets binary, "
            "or install the NCBI datasets CLI manually: "
            "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
        )
    
    # On non-Windows systems, ensure the binary is executable
    if platform.system() != "Windows":
        try:
            os.chmod(datasets_path, os.stat(datasets_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        except OSError as e:
            raise RuntimeError(
                f"Failed to make NCBI datasets binary executable: {e}"
            )
    
    # Verify the bundled binary works
    try:
        result = subprocess.run(
            [datasets_path, "--version"],
            capture_output=True,
            text=True,
            timeout=SUBPROCESS_VERSION_TIMEOUT,
        )
        if result.returncode == 0:
            logger.info(
                "✅ Using bundled NCBI datasets CLI: %s", result.stdout.strip()
            )
            _datasets_path_cache = datasets_path
            return datasets_path
    except (subprocess.TimeoutExpired, OSError) as e:
        raise RuntimeError(
            f"Failed to verify bundled NCBI datasets binary at {datasets_path}: {e}"
        )
    
    raise RuntimeError(
        f"NCBI datasets binary at {datasets_path} failed verification."
    )


def _get_modified_virus_name(virus_name):
    """
    Modify the virus name for retry attempts when the NCBI server is unreachable.
    
    This function generates an alternative virus name to try when the initial
    query fails due to server unreachability. The modification strategy is:
        1. If the name doesn't end with "virus", append " virus" to it.
        2. If the name already ends with "virus", add a space before "virus"
           (e.g., "Denguevirus" -> "Dengue virus").
    
    Args:
        virus_name (str): Original virus name that failed.
        
    Returns:
        str or None: Modified virus name to retry, or None if no modification is possible.
        
    Example:
        >>> _get_modified_virus_name("Dengue")
        'Dengue virus'
        >>> _get_modified_virus_name("Denguevirus")
        'Dengue virus'
        >>> _get_modified_virus_name("Dengue virus")
        None
    """
    if not virus_name:
        return None
    
    virus_lower = virus_name.lower().strip()
    
    # Check if the name already ends with " virus" (with space)
    if virus_lower.endswith(" virus"):
        # Already has proper " virus" suffix, no modification needed
        return None
    
    # Check if the name ends with "virus" (without space before it)
    if virus_lower.endswith("virus"):
        # Add a space before "virus"
        # Find where "virus" starts and insert a space
        idx = virus_name.lower().rfind("virus")
        if idx > 0:
            modified = virus_name[:idx] + " " + virus_name[idx:]
            logger.debug("Modified virus name by adding space before 'virus': '%s' -> '%s'", 
                        virus_name, modified)
            return modified
        return None
    
    # Name doesn't contain "virus" at the end, so append " virus"
    modified = virus_name + " virus"
    logger.debug("Modified virus name by appending ' virus': '%s' -> '%s'", 
                virus_name, modified)
    return modified


def fetch_virus_metadata(
    virus,
    accession=False,
    host=None,
    geographic_location=None,
    annotated=None,
    complete_only=None,
    min_release_date=None,
    refseq_only=None,
    failed_commands=None,
    _is_retry=False,
):
    """
    Fetch virus metadata using NCBI Datasets API.
    
    This function retrieves metadata for virus sequences from the NCBI Datasets API
    using either taxon-based or accession-based queries. It handles pagination
    automatically to retrieve all available results.
    
    When the server is unreachable, this function will automatically retry with
    a modified virus name (adding " virus" suffix or spacing) to improve query
    matching.
    
    Args:
        virus (str): Virus taxon name/ID or accession number.
        accession (bool): Whether virus parameter is an accession number.
        host (str, optional): Host organism name filter.
        geographic_location (str, optional): Geographic location filter.
        annotated (bool, optional): Filter for annotated genomes only.
        complete_only (bool, optional): Filter for complete genomes only.
        min_release_date (str, optional): Minimum release date filter (YYYY-MM-DD format).
        refseq_only (bool, optional): Limit to RefSeq genomes only.
        failed_commands (dict, optional): Dictionary to track failed operations.
        _is_retry (bool): Internal flag to prevent infinite retry loops.
        
    Returns:
        list: List of virus metadata records from the API response.
        
    Raises:
        RuntimeError: If the API request fails.
    """
    
    # Choose the appropriate API endpoint based on whether we're querying by accession or taxon
    if accession:
        # For accession numbers (e.g., NC_045512.2), use the accession-specific endpoint
        url = f"{NCBI_API_BASE}/virus/accession/{virus}/dataset_report"
        logger.debug("Using accession endpoint for virus: %s", virus)
        params = {}
    else:
        # For taxon names/IDs (e.g., 'Zika Virus', 'influenza'), use the taxon endpoint  
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
    
    # Set page size to maximum allowed to minimize the number of API calls needed
    # The NCBI API supports pagination for large result sets
    params['page_size'] = API_PAGE_SIZE
    logger.debug("Set page size to maximum: %d records per request", API_PAGE_SIZE)
    
    # Initialize variables for handling paginated results
    all_reports = []      # Will store all metadata records across all pages
    page_token = None     # Token for accessing subsequent pages
    page_count = 0        # Track number of pages processed for logging
    
    # Main pagination loop - continue until all pages are retrieved
    loop = True
    # Retry logic for handling intermittent server issues
    max_retries = API_MAX_RETRIES
    retry_delay = API_INITIAL_RETRY_DELAY
    last_exception = None
    while loop:
        page_count += 1
        logger.debug("Fetching page %d of results...", page_count)
        
        # Add pagination token if we're not on the first page
        if page_token:
            params['page_token'] = page_token
            
        for attempt in range(max_retries):
            try:
                # Make the HTTP GET request to the NCBI API  
                logger.debug("Making API request to: %s (attempt %d/%d)", url, attempt + 1, max_retries)
                logger.debug("Request parameters: %s", params)
                response = requests.get(url, params=params, timeout=API_REQUEST_TIMEOUT)
                logger.debug("Explicit URL request sent: %s", response.url)
                
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
                    loop = False
                    break
                
                # Set up for the next page
                page_token = next_page_token
                logger.debug("Next page token received, continuing pagination...")
                
                # If we got here, the request succeeded, so break out of retry loop
                break
                
            except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as e:
                last_exception = e
                # Check if this is a retryable error (5xx server errors or connection issues)
                is_retryable = False
                if isinstance(e, requests.exceptions.ConnectionError):
                    is_retryable = True
                elif isinstance(e, requests.exceptions.HTTPError) and hasattr(e, 'response') and e.response:
                    # Check if it's a 5xx server error
                    is_retryable = 500 <= e.response.status_code < 600
                
                if attempt < max_retries - 1 and is_retryable:
                    logger.warning("⚠️ Request failed (attempt %d/%d): %s. Retrying in %.1f seconds...", 
                                 attempt + 1, max_retries, e, retry_delay)
                    time.sleep(retry_delay)
                    retry_delay *= API_RETRY_BACKOFF_MULTIPLIER  # Exponential backoff
                    continue
                else:
                    # Either we've exhausted retries or this is a non-retryable error
                    # Break out of retry loop and handle the error below
                    break
            
            except (requests.exceptions.Timeout, requests.exceptions.RequestException) as e:
                last_exception = e
                # For timeout and other request exceptions, don't retry
                break
        
        # If we have an exception to handle, process it with the original error handling logic
        if last_exception:
            if isinstance(last_exception, requests.exceptions.Timeout):
                # Handle timeout errors with specific guidance for known problematic filters
                error_msg = f"Request timed out while fetching virus metadata: {last_exception}"
                
                # Track API timeout information for summary
                if failed_commands is not None:
                    # api_url = f"{NCBI_API_BASE}?virusTaxId={taxid}" if not accession else f"{NCBI_API_BASE}?accession={accession}"
                    failed_commands['api_timeout'] = {
                        'error': 'API request timeout',
                        'url': url,
                        'alternative_command': None
                    }
                
                if geographic_location:
                    error_msg += (
                        f"\n\n🔧 TIMEOUT LIKELY DUE TO GEOGRAPHIC FILTER: "
                        f"The combination of '{virus}' + geographic location '{geographic_location}' "
                        f"is known to cause server timeouts. Try removing the geographic_location parameter "
                        f"and filter the results manually after download."
                    )
                raise RuntimeError(error_msg) from last_exception
                
            elif isinstance(last_exception, requests.exceptions.ConnectionError):
                # Handle connection errors (network issues, DNS failures, etc.)
                # Try retrying with a modified virus name if this is not already a retry
                if not _is_retry and not accession:
                    modified_virus = _get_modified_virus_name(virus)
                    if modified_virus:
                        logger.warning("⚠️ Connection error with virus name '%s'. "
                                      "Retrying with modified name: '%s'", virus, modified_virus)
                        try:
                            return fetch_virus_metadata(
                                virus=modified_virus,
                                accession=accession,
                                host=host,
                                geographic_location=geographic_location,
                                annotated=annotated,
                                complete_only=complete_only,
                                min_release_date=min_release_date,
                                refseq_only=refseq_only,
                                failed_commands=failed_commands,
                                _is_retry=True,
                            )
                        except RuntimeError:
                            # If retry also fails, raise the original error
                            logger.warning("⚠️ Retry with modified virus name also failed")
                            pass
                raise RuntimeError(f"Connection error while fetching virus metadata: {last_exception}") from last_exception
                
            elif isinstance(last_exception, requests.exceptions.HTTPError):
                # Handle HTTP errors with specific guidance for known issues
                error_msg = f"HTTP error while fetching virus metadata: {last_exception}"
                
                # Check for specific server error patterns (5xx errors indicate server unreachability)
                is_server_error = False
                if hasattr(last_exception, 'response') and last_exception.response is not None:
                    is_server_error = 500 <= last_exception.response.status_code < 600
                elif "500" in str(last_exception) or "502" in str(last_exception) or "503" in str(last_exception) or "504" in str(last_exception):
                    is_server_error = True
                
                if is_server_error:
                    # Special handling for "all viruses" query
                    # If this is the first page and we're querying all viruses without date filters,
                    # the dataset is too large for NCBI to handle - need to chunk by date
                    if virus == NCBI_ALL_VIRUSES_TAXID and not accession and page_count == 1 and not min_release_date:
                        logger.warning("⚠️ NCBI API cannot handle 'all viruses' query in a single request")
                        logger.info("🔄 Automatically switching to date-chunked download strategy...")
                        logger.info("This will split the download into yearly chunks to avoid server overload")
                        # Return None to signal that chunking is needed
                        # The calling function will handle the chunking strategy
                        return None
                    
                    # Try retrying with a modified virus name if this is not already a retry
                    if not _is_retry and not accession:
                        modified_virus = _get_modified_virus_name(virus)
                        if modified_virus:
                            logger.warning("⚠️ Server error (5xx) with virus name '%s'. "
                                          "Retrying with modified name: '%s'", virus, modified_virus)
                            try:
                                return fetch_virus_metadata(
                                    virus=modified_virus,
                                    accession=accession,
                                    host=host,
                                    geographic_location=geographic_location,
                                    annotated=annotated,
                                    complete_only=complete_only,
                                    min_release_date=min_release_date,
                                    refseq_only=refseq_only,
                                    failed_commands=failed_commands,
                                    _is_retry=True,
                                )
                            except RuntimeError:
                                # If retry also fails, continue with original error handling
                                logger.warning("⚠️ Retry with modified virus name also failed")
                                pass
                    
                    error_msg += (
                        f"\n\n🔧 SERVER ERROR DETECTED: "
                        f"NCBI's API is experiencing temporary server-side issues. "
                        f"This is not a problem with your query. Try again in a few minutes, "
                        f"or consider using more specific filters to reduce the dataset size."
                    )
                raise RuntimeError(error_msg) from last_exception
                
            else:
                # Handle any other request-related errors
                raise RuntimeError(f"❌ Failed to fetch virus metadata: {last_exception}") from last_exception

    
    # Log the final results summary
    logger.info("Successfully retrieved %d virus records from NCBI API across %d pages", 
                len(all_reports), page_count)
    
    return all_reports


def fetch_virus_metadata_chunked(
    virus, 
    accession=False, 
    host=None, 
    geographic_location=None, 
    annotated=None,
    complete_only=False,
    min_release_date=None,
    refseq_only=False,
    failed_commands=None
):
    """
    Fetch virus metadata using a chunked date-range strategy for very large datasets.
    
    This function is used as a fallback when the standard fetch_virus_metadata fails
    due to dataset size limitations. It breaks down the request into yearly chunks
    starting from a reasonable start date (2000-01-01 or user's min_release_date) to the present.
    
    Args:
        Same as fetch_virus_metadata.
        
    Returns:
        list: Combined list of virus metadata records from all date chunks.
        
    Raises:
        RuntimeError: If any chunk fails to download.
    """
    
    logger.info("=" * 80)
    logger.info("📦 CHUNKED DOWNLOAD MODE ACTIVATED")
    logger.info("=" * 80)
    logger.info("The 'all viruses' dataset is too large for NCBI to handle in one request.")
    logger.info("Splitting download into yearly chunks to ensure successful completion.")
    logger.info("This may take a while, but ensures all data is retrieved.")
    logger.info("=" * 80)
    
    # Define date range for chunking
    # If user specified min_release_date, use it; otherwise start from default year
    if min_release_date:
        # Extract year from user's min_release_date
        start_year = int(min_release_date.split('-')[0])
        logger.info(f"Starting from user-specified year: {start_year}")
    else:
        # Start from default year as most valuable viral sequence data is from then onwards
        start_year = CHUNKED_DOWNLOAD_START_YEAR
        logger.info("Starting from year %d (default for 'all viruses' downloads)", CHUNKED_DOWNLOAD_START_YEAR)
    
    current_date = datetime.now()
    current_year = current_date.year
    
    all_reports = []
    total_chunks = current_year - start_year + 1
    
    logger.info(f"Will process {total_chunks} year(s) from {start_year} to {current_year}")
    logger.info("=" * 80)
    
    for year in range(start_year, current_year + 1):
        chunk_start = f"{year}-01-01"
        chunk_end = f"{year}-12-31"
        
        # For the current year, use today's date as the end
        if year == current_year:
            chunk_end = current_date.strftime("%Y-%m-%d")
        
        chunk_num = year - start_year + 1
        logger.info("")
        logger.info("=" * 80)
        logger.info(f"📥 Chunk {chunk_num}/{total_chunks}: Fetching data for year {year}")
        logger.info(f"   Date range: {chunk_start} to {chunk_end}")
        logger.info("=" * 80)
        
        try:
            # Fetch metadata for this date chunk
            chunk_reports = fetch_virus_metadata(
                virus=virus,
                accession=accession,
                host=host,
                geographic_location=geographic_location,
                annotated=annotated,
                complete_only=complete_only,
                min_release_date=chunk_start,
                refseq_only=refseq_only,
                failed_commands=failed_commands
            )
            
            # If we got None, it means even this chunk is too large (shouldn't happen for yearly chunks)
            if chunk_reports is None:
                logger.error(f"❌ Chunk for year {year} is too large even for NCBI to handle")
                logger.error("This is unexpected and may indicate an API issue")
                raise RuntimeError(f"Year {year} chunk failed - dataset too large even when split by year")
            
            chunk_count = len(chunk_reports)
            all_reports.extend(chunk_reports)
            
            logger.info(f"✅ Chunk {chunk_num}/{total_chunks} complete: Retrieved {chunk_count:,} records")
            logger.info(f"   Running total: {len(all_reports):,} records")
            
            # Add a small delay between chunks to be respectful to NCBI servers
            if year < current_year:
                time.sleep(CHUNKED_DOWNLOAD_INTER_CHUNK_DELAY)
                
        except Exception as e:
            logger.error(f"❌ Failed to fetch chunk for year {year}: {e}")
            raise RuntimeError(f"Chunked download failed at year {year}") from e
    
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"✅ CHUNKED DOWNLOAD COMPLETE")
    logger.info(f"   Total records retrieved: {len(all_reports):,}")
    logger.info(f"   Total chunks processed: {total_chunks}")
    logger.info("=" * 80)
    
    return all_reports


def is_sars_cov2_query(virus, accession=False):
    """
    Check if the query is for SARS-CoV-2 to determine if optimized downloads should be used.
    
    NCBI provides optimized cached data packages for SARS-CoV-2 that are faster and more
    reliable than the general API endpoints. This function detects SARS-CoV-2 queries.
    
    Args:
        virus (str): Virus taxon name/ID or accession number.
        accession (bool): Whether virus parameter is an accession number.
        
    Returns:
        bool: True if this is a SARS-CoV-2 query that should use cached downloads.
    """
    if accession:
        # When in accession mode, let the user explicitly set is_sars_cov2=True
        # rather than trying to detect it
        return False
    
    # Check for common SARS-CoV-2 identifiers in taxon names
    virus_lower = virus.lower().replace('-', '').replace('_', '').replace(' ', '')
    
    # Check if the query matches any SARS-CoV-2 identifier
    for identifier in SARS_COV2_IDENTIFIERS:
        if identifier in virus_lower:
            logger.debug("Detected SARS-CoV-2 query: %s matches %s", virus, identifier)
            return True
    
    logger.debug("Not a SARS-CoV-2 query: %s", virus)
    return False


def is_alphainfluenza_query(virus, accession=False):
    """
    Check if the query is for Alphainfluenza to determine if optimized downloads should be used.
    
    NCBI provides optimized cached data packages for Alphainfluenza that are faster and more
    reliable than the general API endpoints. This function detects Alphainfluenza queries.
    
    Cached packages are available for:
        - Alphainfluenza (genus, taxid: 197911)
        - Alphainfluenzavirus influenzae (species, taxid: 2955291)
        - Influenza A virus (no-rank, taxid: 11320)
    
    Args:
        virus (str): Virus taxon name/ID or accession number.
        accession (bool): Whether virus parameter is an accession number.
        
    Returns:
        bool: True if this is an Alphainfluenza query that should use cached downloads.
    """
    if accession:
        # When in accession mode, let the user explicitly set is_alphainfluenza=True
        # rather than trying to detect it
        return False
    
    # Check for common Alphainfluenza identifiers in taxon names
    virus_lower = virus.lower().replace('-', '').replace('_', '').replace(' ', '')
    
    # Check if the query matches any Alphainfluenza identifier
    for identifier in ALPHAINFLUENZA_IDENTIFIERS:
        if identifier in virus_lower:
            logger.debug("Detected Alphainfluenza query: %s matches %s", virus, identifier)
            return True
    
    logger.debug("Not an Alphainfluenza query: %s", virus)
    return False


def _process_cached_download(zip_file, virus_type="virus"):
    """
    Process a cached download ZIP file and extract sequences with metadata.
    
    This helper function extracts sequences from a cached ZIP download and loads the
    rich metadata from data_report.jsonl (if available). The metadata is essential 
    for post-download filtering operations.
    
    NCBI cached downloads typically include:
        - genomic.fna: FASTA sequences
        - data_report.jsonl: Rich metadata with virus genome information
        - dataset_catalog.json: List of files in the package
    
    Args:
        zip_file (str): Path to the downloaded ZIP file.
        virus_type (str): Type of virus for logging messages.
        
    Returns:
        tuple: (sequences, metadata_dict, success)
            - sequences: List of all sequence records from the cached download.
            - metadata_dict: Dictionary mapping accessions to metadata (rich metadata
              from data_report.jsonl if available, or basic metadata from FASTA headers).
            - success: Boolean indicating if processing was successful.
            
    Raises:
        RuntimeError: If no valid sequences are found in the cached data.
    """
    if not zip_file or not os.path.exists(zip_file):
        return None, None, False
    
    # Extract directory path from zip file name
    extract_dir = os.path.splitext(zip_file)[0]
    unzip_file(zip_file, extract_dir)
    
    if not os.path.exists(extract_dir):
        logger.warning("Extraction directory not found: %s", extract_dir)
        return None, None, False
    
    logger.info("🔬 PROCESSING CACHED DATA...")
    logger.info("Extracted cached data to: %s", extract_dir)
    
    # Find and load metadata from data_report.jsonl (rich metadata from NCBI)
    metadata_files = []
    fasta_files = []
    for root, dirs, files in os.walk(extract_dir):
        for file in files:
            if file == 'data_report.jsonl':
                metadata_files.append(os.path.join(root, file))
            elif file.endswith(('.fasta', '.fa', '.fna')):
                fasta_files.append(os.path.join(root, file))
    
    # Load rich metadata from data_report.jsonl if available
    cached_metadata_dict = {}
    if metadata_files:
        logger.info("Found %d metadata file(s) in cached download", len(metadata_files))
        for metadata_file in metadata_files:
            try:
                with open(metadata_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        if line.strip():
                            report = json.loads(line)
                            # Extract accession from the report
                            accession = report.get('accession', '')
                            if not accession:
                                continue
                            
                            # Transform the NCBI report format to our internal metadata format
                            # This mirrors the logic in load_metadata_from_api_reports
                            metadata = {
                                'accession': accession,
                                'length': report.get('length'),
                                'source': 'cached_data_report',
                            }
                            
                            # Extract virus info
                            virus_info = report.get('virus', {})
                            metadata['virus_name'] = virus_info.get('sci_name')
                            metadata['virus_tax_id'] = virus_info.get('tax_id')
                            metadata['virus_pangolin_classification'] = virus_info.get('pangolin_classification')
                            
                            # Extract host info
                            host_info = report.get('host', {})
                            metadata['host_name'] = host_info.get('sci_name')
                            metadata['host_common_name'] = host_info.get('common_name')
                            metadata['host_tax_id'] = host_info.get('tax_id')
                            
                            # Extract biosample info
                            biosample_info = report.get('biosample_info', {})
                            metadata['biosample_accession'] = biosample_info.get('accession')
                            
                            # Extract isolate info
                            isolate_info = report.get('isolate', {})
                            metadata['isolate_name'] = isolate_info.get('name')
                            metadata['collection_date'] = isolate_info.get('collection_date')
                            
                            # Extract location info
                            location_info = report.get('location', {})
                            metadata['geo_location'] = location_info.get('geo_location')
                            metadata['usa_state'] = location_info.get('usa_state')
                            
                            # Extract nucleotide completeness
                            completeness_info = report.get('nucleotide_completeness', {})
                            metadata['nuc_completeness'] = completeness_info.get('value')
                            
                            # Extract other fields
                            metadata['release_date'] = report.get('release_date')
                            metadata['update_date'] = report.get('update_date')
                            metadata['is_annotated'] = report.get('annotation', {}).get('is_annotated', False)
                            metadata['is_refseq'] = report.get('is_refseq', False)
                            metadata['is_lab_host'] = report.get('is_lab_host', False)
                            
                            # Gene and protein counts
                            metadata['gene_count'] = report.get('gene_count')
                            metadata['protein_count'] = report.get('protein_count')
                            metadata['mature_peptide_count'] = report.get('mature_peptide_count')
                            
                            cached_metadata_dict[accession] = metadata
                            
                logger.info("Loaded %d metadata records from %s", len(cached_metadata_dict), metadata_file)
            except Exception as e:
                logger.warning("❌ Failed to load metadata file %s: %s", metadata_file, e)
                continue
    else:
        logger.warning("No data_report.jsonl found in cached download. Post-download filters may be limited.")
    
    # Process all available FASTA files
    all_cached_sequences = []
    for fasta_file in fasta_files:
        try:
            sequences = list(FastaIO.parse(fasta_file, "fasta"))
            all_cached_sequences.extend(sequences)
            logger.info("Loaded %d sequences from %s", len(sequences), fasta_file)
        except Exception as e:
            logger.warning("❌ Failed to load FASTA file %s: %s", fasta_file, e)
            continue
    
    if not all_cached_sequences:
        logger.warning("No valid sequences found in cached data.")
        raise RuntimeError("No valid sequences found in cached data")
    
    # If no rich metadata was loaded, create minimal metadata from FASTA headers
    if not cached_metadata_dict:
        logger.info("Creating basic metadata from FASTA headers (no data_report.jsonl available)")
        for seq in all_cached_sequences:
            accession = seq.id.split('.')[0]  # Remove version if present
            cached_metadata_dict[accession] = {
                'accession': accession,
                'description': seq.description,
                'length': len(seq.seq),
                'source': 'cached_fasta_header'
            }
        logger.info("Created basic metadata for %d sequences", len(cached_metadata_dict))
    
    if not all_cached_sequences:
        logger.warning("No valid sequences found in cached data")
        raise RuntimeError("No valid sequences found in cached data")
    
    logger.info("🎉 CACHED DATA LOADING SUCCESSFUL!")
    logger.info("Loaded %d sequences from cached %s data", len(all_cached_sequences), virus_type)
    if metadata_files:
        logger.info("Rich metadata available from data_report.jsonl for post-download filtering")
    
    return all_cached_sequences, cached_metadata_dict, True


def _monitor_subprocess_with_progress(process, cmd, timeout=None, progress_timeout=None):
    """
    Monitor a subprocess with progress tracking and timeout handling.
    
    This helper function monitors a running subprocess, checking for progress
    indicators in stderr output. It implements a two-tier timeout strategy:
        - Overall timeout: Maximum total execution time.
        - Progress timeout: Maximum time without seeing progress.
    
    Args:
        process: subprocess.Popen instance to monitor.
        cmd (list): Command that was executed (for error reporting).
        timeout (int): Maximum total execution time in seconds. Defaults to DOWNLOAD_OVERALL_TIMEOUT.
        progress_timeout (int): Maximum time without progress in seconds. Defaults to DOWNLOAD_PROGRESS_TIMEOUT.
        
    Returns:
        subprocess.CompletedProcess: Result of the completed process.
        
    Raises:
        subprocess.TimeoutExpired: If timeout conditions are met.
    """
    # Apply default timeouts if not specified
    if timeout is None:
        timeout = DOWNLOAD_OVERALL_TIMEOUT
    if progress_timeout is None:
        progress_timeout = DOWNLOAD_PROGRESS_TIMEOUT
    
    start_time = time.time()
    last_progress = start_time
    
    while True:
        # Check if process has finished
        retcode = process.poll()
        if retcode is not None:
            break
            
        # Read stderr without blocking
        stderr = process.stderr.readline()
        if stderr:
            # Log the stderr for debugging
            logger.debug("Progress output: %s", stderr.strip())
            
            # If we see any progress indicator, update the last_progress time
            if any(indicator.lower() in stderr.lower() for indicator in PROGRESS_INDICATORS):
                last_progress = time.time()
                logger.debug("Progress detected, updating last_progress time")
        
        # Check timeout conditions:
        # 1. Less than total timeout, continue
        # 2. If more than total timeout but progress in last progress_timeout, continue
        # 3. Otherwise, timeout
        current_time = time.time()
        total_time = current_time - start_time
        time_since_progress = current_time - last_progress
        
        if total_time > timeout and time_since_progress > progress_timeout:
            process.kill()
            raise subprocess.TimeoutExpired(cmd, timeout)
        
        time.sleep(DOWNLOAD_PROGRESS_CHECK_INTERVAL)  # Prevent CPU spin
    
    stdout, stderr = process.communicate()
    return subprocess.CompletedProcess(
        args=cmd,
        returncode=retcode,
        stdout=stdout,
        stderr=stderr
    )


def _download_optimized_cached(
    virus_type,
    strategies,
    zip_path,
    outdir,
    use_accession=False,
    accession=None,
    requested_filters=None
):
    """
    Execute optimized cached download strategies with fallback.
    
    This is a generic implementation of the hierarchical fallback download pattern
    used for both SARS-CoV-2 and Alphainfluenza. It tries each strategy in order
    until one succeeds, with comprehensive error handling and logging.
    
    Args:
        virus_type (str): Type of virus for error messages ('SARS-CoV-2', 'Alphainfluenza', etc.).
        strategies (list): List of tuples (strategy_name, cmd, applied_filters).
        zip_path (str): Path where ZIP file should be saved.
        outdir (str): Output directory for download.
        use_accession (bool): Whether using accession-based download.
        accession (str, optional): Accession number if using accession-based download.
        requested_filters (dict, optional): Dictionary of originally requested filters.
        
    Returns:
        str: Path to the successfully downloaded ZIP file.
        
    Raises:
        RuntimeError: If all strategies fail or datasets CLI is not available.
        
    Example:
        >>> strategies = [
        ...     ("Strategy 1 (specific)", ["datasets", "download", ...], ["complete-only"]),
        ...     ("Strategy 2 (general)", ["datasets", "download", ...], [])
        ... ]
        >>> zip_file = _download_optimized_cached(
        ...     "SARS-CoV-2", strategies, "/path/to/output.zip", "/output/dir"
        ... )
    """
    
    # Get the path to the datasets CLI binary (uses precompiled binary bundled with gget)
    datasets_path = _get_datasets_path()
    
    last_error = None
    
    for strategy_name, cmd, applied_filters in strategies:
        # Replace "datasets" with the actual path to the binary
        if cmd and cmd[0] == "datasets":
            cmd = [datasets_path] + cmd[1:]
        
        logger.info("🔄 Trying %s...", strategy_name)
        
        if applied_filters:
            logger.debug("Applied filters: %s", ", ".join(applied_filters))
        else:
            logger.debug("No specific filters applied")
        
        logger.debug("Command: %s", " ".join(cmd))
        
        try:
            # Log the exact command being executed
            cmd_str = " ".join(cmd)
            logger.info("📋 Executing command: %s", cmd_str)

            # Start subprocess for progress monitoring
            # Note: We don't use cwd=outdir because the command already includes full paths
            try:
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
            except FileNotFoundError as fnf_error:
                # Datasets binary not found - this shouldn't happen if bundled correctly
                error_msg = (
                    f"❌ datasets binary not found at expected path.\n\n"
                    f"Expected path: {datasets_path}\n\n"
                    "🔧 SOLUTION:\n"
                    "• Please reinstall gget to restore the bundled datasets binary.\n"
                )
                logger.error(error_msg)
                raise RuntimeError(error_msg) from fnf_error
            
            # Monitor progress with timeout handling using helper function
            result = _monitor_subprocess_with_progress(process, cmd)
            
            # Check if the command was successful
            if result.returncode == 0 and os.path.exists(zip_path):
                file_size = os.path.getsize(zip_path)
                logger.info("✅ %s successful: %s (%.2f MB)", 
                           strategy_name, os.path.basename(zip_path), file_size / 1024 / 1024)
                
                # Log any important output from the datasets CLI
                if result.stdout:
                    logger.debug("datasets CLI output: %s", result.stdout.strip())
                
                # Check which filters from the original request weren't applied in this strategy
                if requested_filters:
                    requested_filter_list = []
                    for key, value in requested_filters.items():
                        if value is not None and value is not False:
                            if isinstance(value, bool):
                                requested_filter_list.append(key)
                            else:
                                requested_filter_list.append(f"{key}={value}")
                    
                    missing_filters = [f for f in requested_filter_list if f not in applied_filters]
                    if missing_filters:
                        logger.warning("⚠️ Some requested filters were not applied in successful strategy:")
                        logger.warning("   Filters applied: %s", ", ".join(applied_filters) if applied_filters else "none")
                        logger.warning("   Filters missing: %s", ", ".join(missing_filters))
                        logger.warning("   These filters will need to be applied through post-processing")
                
                return zip_path
            else:
                # Strategy failed, prepare error message
                error_msg = f"❌ {strategy_name} failed with return code {result.returncode}"
                if result.stderr:
                    error_msg += f": {result.stderr.strip()}"
                logger.warning("❌ %s", error_msg)
                last_error = error_msg
                
                # If this was an accession download that failed, provide specific guidance
                if use_accession:
                    error_msg = (
                        f"❌ Failed to download {virus_type} sequence for accession '{accession}'. "
                        f"Please verify that this is a valid {virus_type} accession number. "
                        f"If you're not sure, try without the is_{virus_type.lower().replace('-', '_').replace(' ', '_')} flag."
                    )
                    raise RuntimeError(error_msg)
                
                # Clean up failed download file if it exists
                if os.path.exists(zip_path):
                    try:
                        os.remove(zip_path)
                    except OSError:
                        pass
                continue # Try next strategy
                
        except subprocess.TimeoutExpired:
            error_msg = f"{strategy_name} timed out after 30 minutes"
            logger.warning("❌ %s", error_msg)
            last_error = error_msg
            continue
            
        except subprocess.CalledProcessError as e:
            error_msg = f"{strategy_name} execution failed: {e}"
            logger.warning("❌ %s", error_msg)
            last_error = error_msg
            continue
            
        except Exception as e:
            error_msg = f"{strategy_name} unexpected error: {e}"
            logger.warning("❌ %s", error_msg)
            last_error = error_msg
            continue
    
    # All strategies failed
    logger.warning("🚨 All cached download strategies failed. Last error: %s", last_error)
    
    # Provide helpful guidance based on virus type
    example_taxon = "SARS-CoV-2" if "sars" in virus_type.lower() else virus_type
    guidance_messages = [
        "🔧 TROUBLESHOOTING SUGGESTIONS:",
        "1. Check your internet connection",
        "2. Try running the command manually to see detailed error messages:",
        f"   {datasets_path} download virus genome taxon \"{example_taxon}\" --filename test.zip",
        "3. NCBI servers may be temporarily unavailable - try again later",
        f"4. Consider using the general API method by removing {virus_type} specific terms from your query"
    ]
    
    for msg in guidance_messages:
        logger.info(msg)
    
    # Raise error with the last failure details
    raise RuntimeError(
        f"❌ All {virus_type} cached download strategies failed. "
        f"Last error: {last_error}. "
        f"Consider using the general API method instead."
    )


def download_sars_cov2_optimized(
    host=None,
    complete_only=None,
    annotated=None,
    outdir=None,
    lineage=None,
    accession=None,
    use_accession=False,
):
    """
    Download SARS-CoV-2 sequences using NCBI's optimized cached data packages.
    
    NCBI provides pre-computed, highly compressed cached packages for SARS-CoV-2
    that offer faster and more reliable downloads than the general API endpoints.
    This function uses the datasets CLI to download these optimized packages with
    hierarchical fallback from specific to general cached files.
    
    Download strategies (in order of precedence):
        1. If use_accession=True: Direct accession download using accession endpoint.
        2. If use_accession=False:
           a. Specific lineage + complete + host filters using taxon endpoint.
           b. Complete genomes only using taxon endpoint.
           c. All SARS-CoV-2 genomes using taxon endpoint (default fallback).
    
    Args:
        host (str, optional): Host organism filter (optimized for 'human').
        complete_only (bool, optional): Whether to download only complete genomes.
        annotated (bool, optional): Whether to download only annotated genomes.
        outdir (str, optional): Output directory for downloaded files.
        lineage (str, optional): SARS-CoV-2 lineage filter (e.g., 'B.1.1.7', 'P.1').
        accession (str, optional): Specific SARS-CoV-2 accession or taxon ID.
        use_accession (bool): Whether to use accession endpoint. Defaults to False.
        
    Returns:
        str: Path to the downloaded ZIP file containing sequences and metadata.
        
    Raises:
        RuntimeError: If the datasets CLI is not available or download fails.
    """
    
    # Determine filter specificity for logging
    filter_count = sum(1 for param in [host, complete_only, annotated, lineage] if param is not None)
    if filter_count > 0:
        logger.info("Attempting SARS-CoV-2 cached download with %d specific filters", filter_count)
    else:
        logger.info("Attempting general SARS-CoV-2 cached download (no specific filters)")
    
    # Determine output directory
    if not outdir:
        outdir = os.getcwd()
        logger.debug("No output directory specified, using current directory: %s", outdir)
    
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    logger.debug("Output directory ready: %s", outdir)
    
    # Create descriptive filename with timestamp and random suffix
    zip_filename = f"sars_cov_2_{timestamp}_{random_suffix}.zip"
    zip_path = os.path.join(outdir, zip_filename)
    
    # Define which filters are available for this download
    logger.debug("Available filters for SARS-CoV-2 download:")
    if complete_only:
        logger.debug("- complete-only filter")
    if lineage:
        logger.debug("- lineage filter: %s", lineage)
    if host:
        logger.debug("- host filter: %s", host)
    if annotated:
        logger.debug("- annotated filter")
    
    # Define fallback strategies in order of preference
    strategies = []
    
    if use_accession:
        # Strategy 1: Direct accession download
        cmd1 = ["datasets", "download", "virus", "genome", "accession", accession]
        cmd1.extend(["--filename", zip_path])
        strategies.append(("Strategy 1 (direct accession)", cmd1, [f"accession={accession}"]))
    elif lineage or complete_only or host or annotated:
        # Strategy 1: Try with specific filters using taxon endpoint
        cmd1 = ["datasets", "download", "virus", "genome", "taxon", "SARS-CoV-2"]
        filters1 = []
        
        if complete_only:
            cmd1.append("--complete-only")
            filters1.append("complete-only")
        
        if lineage:
            cmd1.extend(["--lineage", lineage])
            filters1.append(f"lineage={lineage}")

        if host:
            cmd1.extend(["--host", host])
            filters1.append(f"host={host}")
        
        if annotated:
            cmd1.append("--annotated")
            filters1.append("annotated")

        cmd1.extend(["--filename", zip_path])
        strategies.append(("Strategy 1 (specific filters)", cmd1, filters1))
    
    # Strategy 2: Try complete-only and host if it was requested (without lineage)
    if complete_only and host and lineage:  # Only add this if we had lineage in strategy 1
        cmd2 = ["datasets", "download", "virus", "genome", "taxon", "SARS-CoV-2", "--complete-only", "--host", host, "--filename", zip_path]
        strategies.append(("Strategy 2 (complete-only and host)", cmd2, ["complete-only", f"host={host}"]))

    # Strategy 3: Try complete-only if it was requested
    if complete_only and (host or lineage):  
        cmd3 = ["datasets", "download", "virus", "genome", "taxon", "SARS-CoV-2", "--complete-only", "--filename", zip_path]
        strategies.append(("Strategy 3 (complete-only)", cmd3, ["complete-only"]))

    # Strategy 4: Try host if it was requested 
    if host and (complete_only or lineage):  
        cmd4 = ["datasets", "download", "virus", "genome", "taxon", "SARS-CoV-2", "--host", host, "--filename", zip_path]
        strategies.append(("Strategy 4 (host)", cmd4, [f"host={host}"]))

    # Strategy 5: Try lineage if it was requested 
    if lineage and (complete_only or host):  
        cmd5 = ["datasets", "download", "virus", "genome", "taxon", "SARS-CoV-2", "--lineage", lineage, "--filename", zip_path]
        strategies.append(("Strategy 5 (lineage)", cmd5, [f"lineage={lineage}"]))

    # Strategy 6: General SARS-CoV-2 package (no filters)
    cmd6 = ["datasets", "download", "virus", "genome", "taxon", "SARS-CoV-2", "--filename", zip_path]
    strategies.append(("Strategy 6 (general package)", cmd6, []))

    # Use the common download function with all strategies
    requested_filters_dict = {
        'complete-only': complete_only,
        'lineage': lineage,
        'host': host,
        'annotated': annotated
    }
    
    return _download_optimized_cached(
        virus_type="SARS-CoV-2",
        strategies=strategies,
        zip_path=zip_path,
        outdir=outdir,
        use_accession=use_accession,
        accession=accession,
        requested_filters=requested_filters_dict
    )


def download_alphainfluenza_optimized(
    host=None,
    complete_only=None,
    annotated=None,
    outdir=None,
    accession=None,
    use_accession=False,
):
    """
    Download Alphainfluenza sequences using NCBI's optimized cached data packages.
    
    NCBI provides pre-computed, highly compressed cached packages for Alphainfluenza
    that offer faster and more reliable downloads than the general API endpoints.
    This function uses the datasets CLI to download these optimized packages with
    hierarchical fallback from specific to general cached files.
    
    Cached packages are available for the following Alphainfluenza taxonomic nodes:
        1. Alphainfluenza (genus, taxid: 197911)
        2. Alphainfluenzavirus influenzae (species, taxid: 2955291)
        3. Influenza A virus (no-rank, taxid: 11320)
    
    For each taxon, filtered sets are available:
        1. All genomes
        2. Human host only
        3. Human host only & complete
        4. Complete only
    
    Args:
        host (str, optional): Host organism filter (optimized for 'human').
        complete_only (bool, optional): Whether to download only complete genomes.
        annotated (bool, optional): Whether to download only annotated genomes.
        outdir (str, optional): Output directory for downloaded files.
        accession (str, optional): Specific Alphainfluenza accession or taxon ID.
        use_accession (bool): Whether to use accession endpoint. Defaults to False.
        
    Returns:
        str: Path to the downloaded ZIP file containing sequences and metadata.
        
    Raises:
        RuntimeError: If the datasets CLI is not available or download fails.
    """
    
    # Determine filter specificity for logging
    filter_count = sum(1 for param in [host, complete_only, annotated] if param is not None)
    if filter_count > 0:
        logger.info("Attempting Alphainfluenza cached download with %d specific filters", filter_count)
    else:
        logger.info("Attempting general Alphainfluenza cached download (no specific filters)")
    
    # Determine output directory
    if not outdir:
        outdir = os.getcwd()
        logger.debug("No output directory specified, using current directory: %s", outdir)
    
    # Ensure output directory exists before passing path to datasets CLI
    os.makedirs(outdir, exist_ok=True)
    logger.debug("Output directory ready: %s", outdir)
    
    # Create descriptive filename with timestamp and random suffix
    zip_filename = f"alphainfluenza_{timestamp}_{random_suffix}.zip"
    zip_path = os.path.join(outdir, zip_filename)
    
    # Ensure the parent directory exists (in case outdir has subdirectories)
    os.makedirs(os.path.dirname(zip_path), exist_ok=True)
    
    # Define which filters are available for this download
    logger.debug("Available filters for Alphainfluenza download:")
    if complete_only:
        logger.debug("- complete-only filter")
    if host:
        logger.debug("- host filter: %s", host)
    if annotated:
        logger.debug("- annotated filter")
    
    # Define fallback strategies in order of preference
    strategies = []
    
    # Default taxon to use (most specific: Alphainfluenzavirus influenzae species)
    # This taxon ID has the most comprehensive cached data
    default_taxon = ALPHAINFLUENZA_DEFAULT_TAXON
    
    if use_accession:
        # Strategy 1: Direct accession download
        cmd1 = ["datasets", "download", "virus", "genome", "accession", accession]
        cmd1.extend(["--filename", zip_path])
        strategies.append(("Strategy 1 (direct accession)", cmd1, [f"accession={accession}"]))
    elif complete_only or host or annotated:
        # Strategy 1: Try with specific filters using taxon endpoint
        cmd1 = ["datasets", "download", "virus", "genome", "taxon", default_taxon]
        filters1 = []
        
        if complete_only:
            cmd1.append("--complete-only")
            filters1.append("complete-only")

        if host:
            cmd1.extend(["--host", host])
            filters1.append(f"host={host}")
        
        if annotated:
            cmd1.append("--annotated")
            filters1.append("annotated")

        cmd1.extend(["--filename", zip_path])
        strategies.append(("Strategy 1 (specific filters)", cmd1, filters1))
    
    # Strategy 2: Try complete-only and host if both were requested
    if complete_only and host:
        cmd2 = ["datasets", "download", "virus", "genome", "taxon", default_taxon, "--complete-only", "--host", host, "--filename", zip_path]
        strategies.append(("Strategy 2 (complete-only and host)", cmd2, ["complete-only", f"host={host}"]))

    # Strategy 3: Try complete-only if it was requested
    if complete_only and (host or annotated):  
        cmd3 = ["datasets", "download", "virus", "genome", "taxon", default_taxon, "--complete-only", "--filename", zip_path]
        strategies.append(("Strategy 3 (complete-only)", cmd3, ["complete-only"]))

    # Strategy 4: Try host if it was requested 
    if host and (complete_only or annotated):  
        cmd4 = ["datasets", "download", "virus", "genome", "taxon", default_taxon, "--host", host, "--filename", zip_path]
        strategies.append(("Strategy 4 (host)", cmd4, [f"host={host}"]))

    # Strategy 5: General Alphainfluenza package (no filters)
    cmd5 = ["datasets", "download", "virus", "genome", "taxon", default_taxon, "--filename", zip_path]
    strategies.append(("Strategy 5 (general package)", cmd5, []))

    # Use the common download function with all strategies
    requested_filters_dict = {
        'complete-only': complete_only,
        'host': host,
        'annotated': annotated
    }
    
    return _download_optimized_cached(
        virus_type="Alphainfluenza",
        strategies=strategies,
        zip_path=zip_path,
        outdir=outdir,
        use_accession=use_accession,
        accession=accession,
        requested_filters=requested_filters_dict
    )


def download_sequences_by_accessions(accessions, outdir=None, batch_size=200, failed_commands=None):
    """
    Download virus genome sequences for a specific list of accession numbers.
    
    This function downloads sequences for a pre-filtered list of accessions,
    using NCBI E-utilities API with batching to avoid URL length limitations.
    Large requests are automatically split into smaller batches.
    
    Args:
        accessions (list): List of accession numbers to download.
        outdir (str, optional): Output directory for downloaded files.
        batch_size (int): Maximum number of accessions per batch. Defaults to 200.
        failed_commands (dict, optional): Dictionary to track failed operations.
        
    Returns:
        str: Path to the downloaded FASTA file containing sequences.
        
    Raises:
        RuntimeError: If the download request fails.
        ValueError: If no accessions are provided.
    """
    
    if not accessions:
        raise ValueError("No accessions provided for download")
    
    logger.info("Downloading sequences for %d accessions using E-utilities API", len(accessions))
    logger.debug("Accession list: %s", accessions[:5] + ['...'] if len(accessions) > 5 else accessions)
    
    # Determine output directory - use current working directory if not specified
    if not outdir:
        outdir = os.getcwd()
        logger.debug("No output directory specified, using current directory: %s", outdir)
    
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    logger.debug("Ensured output directory exists: %s", outdir)
    
    # Create output FASTA file path
    fasta_path = os.path.join(outdir, f"virus_sequences_{timestamp}_{random_suffix}.fasta")
    logger.debug("Saving sequences to: %s", fasta_path)
    
    # If we have many accessions, split into batches to avoid URL length limits
    if len(accessions) > batch_size:
        logger.info("Large request detected (%d accessions). Using batch processing with batch size %d", 
                   len(accessions), batch_size)
        return _download_sequences_batched(accessions, NCBI_EUTILS_BASE, fasta_path, batch_size, failed_commands)
    
    # For smaller requests, use single request
    return _download_sequences_single_batch(accessions, NCBI_EUTILS_BASE, fasta_path, failed_commands)


def _download_sequences_single_batch(accessions, NCBI_EUTILS_BASE, fasta_path, failed_commands=None):
    """
    Download sequences in a single E-utilities request.
    
    This function handles downloading virus sequences for a list of accessions
    using a single HTTP request to NCBI E-utilities. It's optimized for
    smaller batches (< 200 accessions) to avoid URL length limitations.
    
    Args:
        accessions (list): List of accession numbers to download
        NCBI_EUTILS_BASE (str): Base URL for NCBI E-utilities API
        fasta_path (str): Path where FASTA file should be saved
        
    Returns:
        str: Path to the saved FASTA file
        
    Raises:
        RuntimeError: If the download fails or response is invalid
        
    Note:
        - Validates FASTA format before saving
        - Includes extended timeout for large datasets
        - Automatically retries with batching if URL is too long
        
    Example:
        >>> accessions = ['NC_045512.2', 'MN908947.3']
        >>> path = _download_sequences_single_batch(accessions, BASE_URL, 'output.fasta')
    """
    
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
        logger.debug("E-utilities URL: %s", NCBI_EUTILS_BASE)
        
        # Make the request with extended timeout for large datasets
        response = requests.get(NCBI_EUTILS_BASE, params=params, timeout=EUTILS_TIMEOUT)
        
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
        logger.error("❌ E-utilities request failed: %s", e)
        
        # Check for specific URL length error
        if "414" in str(e) or "Request-URI Too Long" in str(e):
            logger.info("URL too long error detected. Retrying with batch processing...")
            # Retry with smaller batches (half of default)
            return _download_sequences_batched(accessions, NCBI_EUTILS_BASE, fasta_path, batch_size=EUTILS_DEFAULT_BATCH_SIZE // 2, failed_commands=failed_commands)
        raise RuntimeError(f"❌ Failed to download virus sequences via E-utilities: {e}") from e
    except IOError as e:
        logger.error("❌ Failed to save FASTA file: %s", e)
        raise RuntimeError(f"❌ Failed to save downloaded sequences: {e}") from e


def _download_sequences_batched(accessions, NCBI_EUTILS_BASE, fasta_path, batch_size, failed_commands=None):
    """
    Download sequences using multiple batched E-utilities requests with incremental file writing.
    
    This function handles large sequence downloads by splitting them into smaller
    batches and writing results incrementally to avoid memory issues. It includes
    robust error handling with automatic retries for failed batches.
    
    Key features:
    - Batched requests to avoid URL length limits
    - Incremental file writing to manage memory
    - Automatic retry with smaller batch sizes for failures
    - Progress tracking and detailed logging
    - Graceful handling of partial failures
    
    Args:
        accessions (list): List of accession numbers to download
        NCBI_EUTILS_BASE (str): Base URL for NCBI E-utilities API
        fasta_path (str): Path where FASTA file should be saved
        batch_size (int): Number of accessions per batch
        
    Returns:
        str: Path to the saved FASTA file containing all downloaded sequences
        
    Raises:
        RuntimeError: If all batches fail or no sequences are downloaded
        
    Note:
        - Respects NCBI rate limits with 0.5s delays between batches
        - Automatically reduces batch size for URL length errors
        - Continues processing even if some batches fail
        - Writes sequences immediately to reduce memory usage
        
    Example:
        >>> large_accession_list = ['NC_045512.2', 'MN908947.3', ...]  # 1000+ accessions
        >>> path = _download_sequences_batched(large_accession_list, BASE_URL, 'out.fasta', 200)
    """
    
    # Split accessions into batches
    batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
    logger.info("Downloading %d accessions in %d batches of size %d", 
               len(accessions), len(batches), batch_size)
    
    total_downloaded = 0
    batch_failed_count = 0
    
    # Open file once and write batches incrementally to avoid storing all data in memory
    try:
        with open(fasta_path, 'w', encoding='utf-8') as f:
            for batch_num, batch_accessions in enumerate(batches, 1):
                logger.info("Processing batch %d/%d (%d accessions)", 
                           batch_num, len(batches), len(batch_accessions))
                
                # Build accession string for this batch
                accession_string = ",".join(batch_accessions)
                params = {
                    'db': 'nucleotide',
                    'id': accession_string,
                    'rettype': 'fasta',
                    'retmode': 'text'
                }
                
                try:
                    # Make the request with timeout
                    response = requests.get(NCBI_EUTILS_BASE, params=params, timeout=EUTILS_TIMEOUT)
                    response.raise_for_status()
                    
                    # Verify we got FASTA data
                    if not response.text.strip().startswith('>'):
                        logger.warning("Batch %d returned invalid FASTA data: %s", 
                                     batch_num, response.text[:100])
                        batch_failed_count += 1
                        continue
                    
                    # Count sequences in this batch
                    batch_sequence_count = response.text.count('>')
                    total_downloaded += batch_sequence_count
                    
                    # Write sequences immediately to file (incremental write)
                    f.write(response.text)
                    if not response.text.endswith('\n'):
                        f.write('\n')  # Ensure proper line endings between batches
                    f.flush()  # Force write to disk immediately
                    
                    logger.info("Batch %d: Downloaded and wrote %d sequences (%.2f MB)", 
                               batch_num, batch_sequence_count, 
                               len(response.text.encode('utf-8')) / BYTES_PER_MB)
                    
                    # Add small delay between requests to be respectful to NCBI servers
                    if batch_num < len(batches):  # Don't delay after the last batch
                        time.sleep(EUTILS_INTER_BATCH_DELAY)
                        
                except requests.exceptions.RequestException as e:
                    logger.error("Batch %d failed: %s", batch_num, e)
                    batch_failed_count += 1
                    
                    # Track failed batch for summary (if tracking is enabled)
                    if failed_commands is not None:
                        retry_url = f"{NCBI_EUTILS_BASE}?db=nucleotide&id={','.join(batch_accessions)}&rettype=fasta&retmode=text"
                        failed_commands['sequence_batches'].append({
                            'batch_num': batch_num,
                            'accessions': batch_accessions,
                            'error': str(e),
                            'retry_url': retry_url
                        })
                    
                    # Check for URL length error even in batch mode
                    if "414" in str(e) and batch_size > EUTILS_MIN_BATCH_SIZE_FOR_SPLIT:
                        logger.warning("⚠️ URL still too long for batch size %d. Retrying batch %d with smaller size...", 
                                     batch_size, batch_num)
                        # Recursively retry this batch with smaller size by splitting it further
                        temp_batch_path = f"temp_batch_{batch_num}.fasta"
                        try:
                            _download_sequences_batched(
                                batch_accessions, NCBI_EUTILS_BASE, temp_batch_path, batch_size // 2, failed_commands
                            )
                            # Read the temporary file and append to main file
                            with open(temp_batch_path, 'r', encoding='utf-8') as temp_f:
                                batch_content = temp_f.read()
                                f.write(batch_content)
                                if not batch_content.endswith('\n'):
                                    f.write('\n')
                                f.flush()
                                # Count sequences in this recovered batch
                                recovered_count = batch_content.count('>')
                                total_downloaded += recovered_count
                                batch_failed_count -= 1  # This batch succeeded after retry
                                logger.info("Recovered batch %d with smaller size: %d sequences", 
                                           batch_num, recovered_count)
                            os.remove(temp_batch_path)  # Clean up temp file
                        except Exception as file_error:
                            logger.warning("❌ Failed to recover batch %d with smaller size: %s", batch_num, file_error)
                            continue
                    else:
                        logger.warning("❌ Batch %d failed and will be skipped: %s", batch_num, e)
                        continue
        
        # Check if we downloaded anything
        if total_downloaded == 0:
            raise RuntimeError("❌ All batches failed. No sequences were downloaded.")
        
        if batch_failed_count > 0:
            logger.warning("⚠️ %d out of %d batches failed. Successfully downloaded %d sequences.", 
                          batch_failed_count, len(batches), total_downloaded)
        
        file_size = os.path.getsize(fasta_path)
        logger.info("Successfully saved %d sequences to: %s (%.2f MB)", 
                   total_downloaded, fasta_path, file_size / BYTES_PER_MB)
        return fasta_path
        
    except IOError as e:
        logger.error("❌ Failed to write FASTA file: %s", e)
        raise RuntimeError(f"❌ Failed to save downloaded sequences: {e}") from e


def unzip_file(zip_file_path, extract_to_path):
    """
    Unzip a ZIP file to a specified directory.
    
    Args:
        zip_file_path (str): Path to the ZIP file to extract.
        extract_to_path (str): Directory where contents will be extracted.
        
    Raises:
        zipfile.BadZipFile: If the ZIP file is invalid or corrupted.
        PermissionError: If there are permission issues with the target directory.
        FileNotFoundError: If the ZIP file does not exist.
    """
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
                "accession": accession,
                "length": report.get("length"),  # Sequence length in nucleotides
                "geneCount": report.get("gene_count"),  # Number of genes annotated
                "completeness": report.get("completeness", "").lower(), # Completeness status (e.g., complete, partial)
                "host": report.get("host", {}),  # Host organism details
                "isLabHost": report.get("host", {}).get("is_lab_host", False),  # Lab-passaged flag
                "labHost": report.get("host", {}).get("is_lab_host", False),  # Alternative field name
                "location": report.get("location", {}).get("geographic_location", None),  # Geographic location details
                "region": report.get("location", {}).get("geographic_region", None),  # Broad region
                "submitter": report.get("submitters", [{}])[0] if report.get("submitters") else {},
                "sourceDatabase": report.get("source_database", ""),  # GenBank, RefSeq, etc.
                "isolate": report.get("isolate", {}),  # Sample/isolate details
                "virus": report.get("virus", {}),  # Virus taxonomy and classification
                "isAnnotated": report.get("is_annotated", False),  # Whether sequence is annotated
                "releaseDate": report.get("release_date", ""),  # When sequence was released
                "sraAccessions": report.get("sra_accessions", []),  # SRA read data accessions 
                "bioprojects": report.get("bioprojects", []),  # Associated BioProject IDs 
                "biosample": report.get("biosample"),  # BioSample ID 
                "proteinCount": report.get("protein_count"),  # Number of proteins
                "maturePeptideCount": report.get("mature_peptide_count"),  # Number of mature peptides
            }
            
            # Store the metadata using accession as the key
            metadata_dict[accession] = metadata
            logger.debug("Processed metadata for accession: %s (length: %s, host: %s)", 
                        accession, 
                        metadata.get("length"), 
                        metadata.get("host", {}).get("organism_name", "Unknown"))
            
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
        # Default to a very early year for partial dates to handle edge cases properly
        parsed_date = parser.parse(date_str, default=datetime(DATE_PARSE_DEFAULT_YEAR, 1, 1))
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
            logger.error("❌ Date parsing failed: %s", error_msg)
            raise ValueError(error_msg) from exc
        else:
            # In non-verbose mode, log warning and return None
            logger.warning("⚠️ Failed to parse date '%s' for filter '%s': %s", 
                          date_str, filtername, exc)
            return None


def _check_protein_requirements(record, metadata, has_proteins, proteins_complete):
    """
    Check if a sequence meets protein/gene requirements.
    
    This function validates whether a virus sequence contains required proteins
    or genes, and optionally whether those proteins are complete. It checks both
    the FASTA header (for gene/protein names) and metadata (for counts).
    
    Args:
        record: FastaRecord object containing the sequence and description
        metadata (dict): Metadata dictionary for this accession
        has_proteins (str/list/None): Required protein(s)/gene(s) to check for
                                      Can be a single string or list of strings
        proteins_complete (bool): Whether to require complete protein annotations
        
    Returns:
        bool: True if protein requirements are met, False otherwise
        
    Example:
        >>> # Check for spike protein
        >>> _check_protein_requirements(record, metadata, "spike", False)
        True
        >>> # Check for multiple proteins
        >>> _check_protein_requirements(record, metadata, ["HA", "NA"], False)
        True
    """
    
    # If proteins_complete is True, check if sequence has protein annotations
    if proteins_complete:
        protein_count = metadata.get("proteinCount", 0)
        gene_count = metadata.get("geneCount", 0)
        
        # Consider a sequence complete if it has at least one protein or gene annotated
        if protein_count is None or protein_count == 0:
            if gene_count is None or gene_count == 0:
                logger.debug("Sequence %s has no protein/gene annotations (protein_count=%s, gene_count=%s)",
                           record.id, protein_count, gene_count)
                return False
    
    # If specific proteins are required, check for them in the header
    if has_proteins is not None:
        # Convert single string to list for uniform processing
        required_proteins = [has_proteins] if isinstance(has_proteins, str) else has_proteins
        
        # Get the full FASTA header description
        header = record.description.lower() if hasattr(record, 'description') else str(record.id).lower()
        
        # Check if all required proteins are mentioned in the header
        for protein in required_proteins:
            protein_lower = protein.lower()
            if protein_lower not in header:
                logger.debug("Sequence %s missing required protein: %s (header: %s)",
                           record.id, protein, record.description[:100])
                return False
        
        logger.debug("Sequence %s has all required proteins: %s",
                   record.id, required_proteins)
    
    return True


def _extract_protein_info_from_header(header):
    """
    Extract protein/segment information from FASTA header.
    
    This function parses FASTA headers to extract useful protein or segment
    information, which is particularly important for segmented viruses like
    influenza that have multiple genomic segments encoding different proteins.
    
    Common patterns extracted:
    - Segment numbers (e.g., "segment 4", "segment 8")
    - Gene/protein names (e.g., "hemagglutinin", "HA", "neuraminidase", "NA")
    - CDS/protein annotations
    
    Args:
        header (str): FASTA header/description line
        
    Returns:
        str: Extracted protein/segment information, or empty string if none found
        
    Example:
        >>> _extract_protein_info_from_header("A/H1N1 segment 4 (HA)")
        "segment 4 (HA)"
        >>> _extract_protein_info_from_header("spike glycoprotein gene")
        "spike glycoprotein"
    """
    
    if not header:
        return ""
    
    header_lower = header.lower()
    
    # Try to find any of these keywords in the header
    found_info = []
    for keyword in PROTEIN_KEYWORDS:
        if keyword in header_lower:
            # Find the position and extract surrounding context
            pos = header_lower.find(keyword)
            # Extract a reasonable chunk around the keyword
            start = max(0, pos - 10)
            end = min(len(header), pos + len(keyword) + 20)
            chunk = header[start:end].strip()
            
            # Clean up the chunk
            chunk = chunk.split('|')[0].strip()  # Remove accession info after |
            chunk = chunk.split(',')[0].strip()   # Remove comma-separated info
            
            if chunk and chunk not in found_info:
                found_info.append(chunk)
    
    if found_info:
        return "; ".join(found_info)
    
    # If no specific keywords found, try to extract gene/protein from common patterns
    # Pattern: "gene for X" or "X gene" or "X protein"
    import re
    gene_pattern = r'(?:gene for |protein )?([a-zA-Z0-9\-]+(?:\s+[a-zA-Z0-9\-]+){0,2})(?:\s+gene|\s+protein)'
    match = re.search(gene_pattern, header_lower)
    if match:
        return match.group(1).strip()
    
    return ""


def filter_sequences(
    fna_file,
    metadata_dict,
    max_ambiguous_chars=None,
    has_proteins=None,
    proteins_complete=False,
):
    """
    Apply sequence-dependent filters to downloaded sequences.
    
    This function only applies filters that require the actual sequence data:
    - Ambiguous character counting
    - Protein/feature analysis if required
    
    Note: All metadata-only filters should have been applied by filter_metadata_only
    before downloading sequences. The metadata-related parameters are kept for
    backwards compatibility but are ignored.
    
    Args:
        fna_file (str): Path to FASTA file containing sequences
        metadata_dict (dict): Dictionary mapping accession numbers to metadata
        max_ambiguous_chars (int): Maximum number of ambiguous nucleotides allowed
        has_proteins (str/list): Required proteins/genes filter
        proteins_complete (bool): Whether proteins must be complete
        
    Returns:
        tuple: (filtered_sequences, filtered_metadata, protein_headers)
    """
    logger.info("Applying sequence-dependent filters...")
    logger.debug("Sequence filters: max_ambiguous=%s, proteins=%s, complete=%s",
                max_ambiguous_chars, has_proteins, proteins_complete)
    
    # Initialize lists to store filtered results
    filtered_sequences = []    # Will store FastaRecord objects that pass filters
    filtered_metadata = []     # Will store corresponding metadata dictionaries
    protein_headers = []       # Will store protein/segment information from FASTA headers
    
    # Counters for logging filter statistics
    total_sequences = 0
    filter_stats = {
        'seq_length': 0,
        'ambiguous_chars': 0,
        'proteins': 0
    }

    # Read and process sequences from the FASTA file
    logger.info("Reading sequences from FASTA file: %s", fna_file)
    for record in FastaIO.parse(fna_file, "fasta"):
        total_sequences += 1
        record_passes = True
            
        # Count ambiguous characters (N's)
        if max_ambiguous_chars is not None:
            ambiguous_count = record.seq.upper().count('N')
            if ambiguous_count > max_ambiguous_chars:
                filter_stats['ambiguous_chars'] += 1
                record_passes = False
                continue
        
        # Get metadata for this record to check protein information
        record_metadata = metadata_dict.get(record.id, {})
        
        # Check protein requirements if specified
        if has_proteins is not None or proteins_complete:
            protein_check_passed = _check_protein_requirements(
                record, 
                record_metadata, 
                has_proteins, 
                proteins_complete
            )
            
            if not protein_check_passed:
                filter_stats['proteins'] += 1
                record_passes = False
                logger.debug("Sequence %s failed protein requirements", record.id)
                continue
        
        # If sequence passed all filters, keep it and its metadata
        if record_passes:
            filtered_sequences.append(record)
            filtered_metadata.append(record_metadata)
            
            # Extract protein/segment information from FASTA header for CSV output
            # This is useful for segmented viruses like influenza
            protein_info = _extract_protein_info_from_header(record.description)
            protein_headers.append(protein_info)

    # Log filtering results
    logger.info("Sequence filter results:")
    logger.info("- Total sequences processed: %d", total_sequences)
    logger.info("- ⚠️ Failed length filter: %d", filter_stats['seq_length'])
    logger.info("- ⚠️ Failed ambiguous char filter: %d", filter_stats['ambiguous_chars'])
    logger.info("- ⚠️ Failed protein requirements: %d", filter_stats['proteins'])
    logger.info("- Sequences passing all filters: %d", len(filtered_sequences))
    
    return filtered_sequences, filtered_metadata, protein_headers


def save_command_summary(
    outfolder,
    command_line,
    total_api_records,
    total_after_metadata_filter,
    total_final_sequences,
    output_files,
    filtered_metadata,
    success=True,
    error_message=None,
    failed_commands=None,
    genbank_error=None
):
    """
    Save a summary file documenting the command execution and results.
    
    This function creates a comprehensive summary text file that includes:
    - The exact command line that was run
    - Statistics about sequences, hosts, locations, etc.
    - List of output files generated
    - Information about what filters were applied
    - Any errors or warnings encountered
    - Failed download commands (batches, URLs) for user retry
    
    Args:
        outfolder (str): Output directory where summary file will be saved
        command_line (str): The command line that was executed
        total_api_records (int): Number of records returned from API
        total_after_metadata_filter (int): Number after metadata filtering
        total_final_sequences (int): Final number of sequences after all filters
        output_files (dict): Dictionary of output file paths
        filtered_metadata (list): List of metadata dictionaries for statistics
        success (bool): Whether the command completed successfully
        error_message (str): Error message if command failed
        failed_commands (dict): Dictionary containing failed operations with retry commands/URLs
            Expected keys: 'api_timeout', 'sequence_batches', 'genbank_batches'
        genbank_error (str): Error message if GenBank metadata retrieval failed
    """
    
    summary_file = os.path.join(outfolder, "command_summary.txt")
    
    try:
        with open(summary_file, 'w', encoding='utf-8') as f:
            # Header
            f.write("=" * 80 + "\n")
            f.write("GGET VIRUS COMMAND SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            # Timestamp
            f.write(f"Execution Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Output Folder: {outfolder}\n\n")
            
            # Command line
            f.write("-" * 80 + "\n")
            f.write("COMMAND LINE\n")
            f.write("-" * 80 + "\n")
            f.write(f"{command_line}\n\n")
            
            # Execution status
            f.write("-" * 80 + "\n")
            f.write("EXECUTION STATUS\n")
            f.write("-" * 80 + "\n")
            if success:
                if genbank_error:
                    f.write("Command completed with warnings\n")
                    f.write(f"⚠️ GenBank metadata retrieval failed: {genbank_error}\n\n")
                else:
                    f.write("✅ Command completed successfully\n\n")
            else:
                f.write("✗ Command failed\n")
                if error_message:
                    f.write(f"Error: {error_message}\n\n")
            
            # Statistics
            f.write("-" * 80 + "\n")
            f.write("SEQUENCE STATISTICS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Total records from API: {total_api_records}\n")
            f.write(f"After metadata filtering: {total_after_metadata_filter}\n")
            f.write(f"Final sequences (after all filters): {total_final_sequences}\n\n")
            
            # Detailed statistics from metadata
            if filtered_metadata and len(filtered_metadata) > 0:
                f.write("-" * 80 + "\n")
                f.write("DETAILED STATISTICS\n")
                f.write("-" * 80 + "\n")
                
                # Unique hosts
                hosts = set()
                for meta in filtered_metadata:
                    host_name = meta.get('host', {}).get('organism_name') if isinstance(meta.get('host'), dict) else None
                    if host_name:
                        hosts.add(host_name)
                f.write(f"Unique hosts: {len(hosts)}\n")
                if hosts and len(hosts) <= 20:
                    for host in sorted(hosts):
                        f.write(f"  - {host}\n")
                elif hosts:
                    f.write(f"  (Top 20 of {len(hosts)} hosts)\n")
                    for host in sorted(hosts)[:20]:
                        f.write(f"  - {host}\n")
                f.write("\n")
                
                # Unique geographic locations
                locations = set()
                for meta in filtered_metadata:
                    location = meta.get('location', {}).get('geographic_location') if isinstance(meta.get('location'), dict) else None
                    if location:
                        locations.add(location)
                f.write(f"Unique geographic locations: {len(locations)}\n")
                if locations and len(locations) <= 20:
                    for loc in sorted(locations):
                        f.write(f"  - {loc}\n")
                elif locations:
                    f.write(f"  (Top 20 of {len(locations)} locations)\n")
                    for loc in sorted(locations)[:20]:
                        f.write(f"  - {loc}\n")
                f.write("\n")
                
                # Sequence length statistics
                lengths = [meta.get('length') for meta in filtered_metadata if meta.get('length')]
                if lengths:
                    f.write(f"Sequence length range: {min(lengths)} - {max(lengths)} bp\n")
                    f.write(f"Average sequence length: {sum(lengths) / len(lengths):.0f} bp\n\n")
                
                # Completeness breakdown
                completeness_counts = {}
                for meta in filtered_metadata:
                    comp = meta.get('completeness', 'unknown')
                    completeness_counts[comp] = completeness_counts.get(comp, 0) + 1
                f.write("Completeness breakdown:\n")
                for comp, count in sorted(completeness_counts.items()):
                    f.write(f"  - {comp}: {count}\n")
                f.write("\n")
                
                # Source database breakdown
                source_counts = {}
                for meta in filtered_metadata:
                    source = meta.get('sourceDatabase', 'unknown')
                    source_counts[source] = source_counts.get(source, 0) + 1
                f.write("Source database breakdown:\n")
                for source, count in sorted(source_counts.items()):
                    f.write(f"  - {source}: {count}\n")
                f.write("\n")
                
                # Submitter countries
                countries = set()
                for meta in filtered_metadata:
                    submitter = meta.get('submitter', {})
                    if isinstance(submitter, dict):
                        country = submitter.get('country')
                        if country:
                            countries.add(country)
                f.write(f"Unique submitter countries: {len(countries)}\n")
                if countries and len(countries) <= 20:
                    for country in sorted(countries):
                        f.write(f"  - {country}\n")
                elif countries:
                    f.write(f"  (Top 20 of {len(countries)} countries)\n")
                    for country in sorted(countries)[:20]:
                        f.write(f"  - {country}\n")
                f.write("\n")
            
            # Output files
            f.write("-" * 80 + "\n")
            f.write("OUTPUT FILES\n")
            f.write("-" * 80 + "\n")
            if output_files:
                for file_type, file_path in output_files.items():
                    if file_path and os.path.exists(file_path):
                        file_size = os.path.getsize(file_path) / 1024 / 1024  # MB
                        f.write(f"{file_type}: {os.path.basename(file_path)} ({file_size:.2f} MB)\n")
                    elif file_path:
                        f.write(f"{file_type}: {os.path.basename(file_path)} (not created)\n")
            else:
                f.write("No output files generated\n")
            f.write("\n")
            
            # Footer
            f.write("=" * 80 + "\n")
            f.write("END OF SUMMARY\n")
            f.write("=" * 80 + "\n")
        
        logger.info("✅ Command summary saved: %s", summary_file)
        return summary_file
        
    except Exception as e:
        logger.error("Failed to save command summary: %s", e)
        logger.error("Traceback: %s", traceback.format_exc())
        return None


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
            "Geographic Region": metadata.get("region", pd.NA),
            "Geographic Location": metadata.get("location", pd.NA),
            
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
        logger.error("❌ Failed to save CSV file: %s", e)
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
                logger.error("❌ Failed to parse dates for validation: %s", e)
                raise ValueError(f"Invalid date format in {filtername} filters") from e
        
        # Check if minimum is greater than maximum
        if min_val > max_val:
            error_msg = f"Min value ({min_val}) cannot be greater than max value ({max_val}) for {filtername}."
            logger.error("❌ Validation failed: %s", error_msg)
            raise ValueError(error_msg)
        
        logger.debug("Min/max validation passed for %s", filtername)


def fetch_genbank_metadata(accessions, genbank_full_xml_path, genbank_full_csv_path, batch_size=200, delay=0.5, failed_log_path=None):
    """
    Fetch detailed GenBank metadata for a list of accession numbers using NCBI E-utilities.
    
    This function provides an optimized alternative to Biopython's approach for retrieving
    GenBank information. It uses NCBI's E-utilities API directly with HTTP requests to
    fetch GenBank records in XML format, then parses them using Python's built-in XML
    library to extract comprehensive metadata.
    
    Key advantages over Biopython approach:
    - No external dependencies beyond standard library
    - Batch processing for improved performance
    - Structured output suitable for CSV/analysis
    - Respectful rate limiting for NCBI servers
    - Comprehensive error handling and logging
    
    Extracted metadata includes:
    - Basic sequence information (organism, length, definition)
    - Collection metadata (date, location, host, strain)
    - Publication references (authors, titles, journals, PubMed IDs)
    - Database information (create/update dates, comments)
    - Taxonomic classification
    - Assembly information
    
    Args:
        accessions (list): List of accession numbers to fetch GenBank data for
        batch_size (int): Maximum number of accessions per API request (default: 200)
                         Recommended range: 50-500 depending on server load
        delay (float): Delay in seconds between batch requests (default: 0.5)
                      Helps avoid overloading NCBI servers
        
    Returns:
        dict: Dictionary mapping accession numbers to GenBank metadata dictionaries
              Key: accession number (str) 
              Value: dictionary with structure:
                     {
                         'accession': str,
                         'genbank_data': {
                             'organism': str,
                             'sequence_length': int,
                             'collection_date': str,
                             'country': str,
                             'host': str,
                             'references': [{'title': str, 'authors': str, ...}],
                             ... (20+ additional fields)
                         }
                     }
              
    Raises:
        RuntimeError: If API requests fail or XML parsing encounters errors
        ValueError: If no accessions are provided
        
    Example:
        >>> accessions = ['NC_045512.2', 'MN908947.3']
        >>> genbank_data = fetch_genbank_metadata(accessions)
        >>> print(genbank_data['NC_045512.2']['genbank_data']['organism'])
        'Severe acute respiratory syndrome coronavirus 2'
    """
    if failed_log_path is None:
        failed_log_path = os.path.join(os.path.dirname(genbank_full_xml_path), "genbank_failed_batches.log")
    if os.path.exists(failed_log_path):
        os.remove(failed_log_path)

    if not accessions:
        raise ValueError("No accessions provided for GenBank metadata retrieval")
    
    logger.info("Fetching GenBank metadata for %d accessions using E-utilities", len(accessions))
    logger.debug("First 5 accessions: %s", accessions[:5])
    
    # Initialize tracking variables
    all_metadata = {}
    all_xml_text = ""
    failed_batches = []
    
    # Split accessions into batches to avoid URL length limits and server overload
    if len(accessions) > batch_size:
        batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
        logger.info("Processing %d accessions in %d batches of size %d", 
                   len(accessions), len(batches), batch_size)
    else:
        batches = [accessions]
        logger.info("Processing %d accessions in 1 batch", len(accessions))
    
    # Process each batch
    for batch_num, batch_accessions in enumerate(batches, 1):
        logger.info("Processing GenBank batch %d/%d (%d accessions)", 
                   batch_num, len(batches), len(batch_accessions))
        
        try:
            # Fetch GenBank XML data using E-utilities efetch
            batch_metadata, batch_xml_text = _fetch_genbank_batch(batch_accessions, failed_log_path=failed_log_path)
            
            if batch_metadata:
                all_metadata.update(batch_metadata)

            if batch_xml_text:
                # Clean XML before concatenating using helper function
                cleaned_xml = _clean_xml_declarations(batch_xml_text)
                
                # Add to the global XML string
                all_xml_text += cleaned_xml + "\n"
                
                logger.info("Batch %d: Successfully retrieved metadata for %d accessions", 
                           batch_num, len(batch_metadata))
            else:
                # Batch failed, add to failed_batches for retry
                logger.warning("Batch %d returned no data, will retry later", batch_num)
                failed_batches.append(batch_accessions)
            
            # Add delay between requests to be respectful to NCBI servers
            if batch_num < len(batches) and delay > 0:
                logger.debug("Adding %.1f second delay before next batch", delay)
                time.sleep(delay)  # Uses the delay parameter passed to the function
                
        except Exception as e:
            logger.error("⚠️ Batch %d failed: %s", batch_num, e)
            failed_batches.append(batch_accessions)
            logger.info("Added batch %d to retry list", batch_num)
            # Continue with remaining batches rather than failing completely
            continue
    
    # Retry failed batches at the end
    if failed_batches:
        logger.info("Retrying %d failed batches", len(failed_batches))
        retry_success = []
        for batch_accessions in failed_batches:
            try:
                meta, xml = _fetch_genbank_batch(batch_accessions, failed_log_path=failed_log_path)
                if meta:
                    all_metadata.update(meta)
                    retry_success.append(batch_accessions)
                if xml:
                    # Clean XML using helper function
                    cleaned_xml = _clean_xml_declarations(xml)
                    all_xml_text += cleaned_xml + "\n"
                    logger.info("Successfully retried batch with %d accessions", len(batch_accessions))
            except Exception as e:
                logger.warning("Final retry failed for batch: %s", batch_accessions)
        
        if retry_success:
            logger.info("Successfully recovered %d/%d failed batches on retry", len(retry_success), len(failed_batches))

    logger.info("GenBank metadata retrieval complete: %d/%d accessions processed", 
                len(all_metadata), len(accessions))

    if all_xml_text:
        final_xml = "<AllGBSets>\n" + all_xml_text + "</AllGBSets>"
        logger.debug("Saving full GenBank XML and CSV data to: %s and %s", genbank_full_xml_path, genbank_full_csv_path)
        _save_genbank_xml_and_csv(final_xml, genbank_full_xml_path, genbank_full_csv_path)
    else:
        logger.warning("No GenBank XML content retrieved to save")

    if not all_metadata:
        logger.warning("No GenBank metadata was successfully retrieved")
    
    missing_accessions = set(accessions) - set(all_metadata.keys())
    if missing_accessions:
        logger.info("❌ The following accessions could not be downloaded:")
        for acc in missing_accessions:
            efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={acc}&rettype=gb&retmode=xml"
            logger.info(f"  {acc}: {efetch_url}")
        logger.info(f"A log of failed batches and efetch URLs is saved at: {failed_log_path}")
    
    # Return both metadata and the path to the failed batches log for summary tracking
    return all_metadata, failed_log_path if os.path.exists(failed_log_path) else None


def _fetch_genbank_batch(accessions, failed_log_path=None):
    """
    Fetch GenBank metadata for a single batch of accessions.
    
    This function handles retrieval of GenBank XML data for a batch of accessions
    using NCBI E-utilities. It includes comprehensive error handling with retries,
    exponential backoff, and automatic batch splitting for problematic requests.
    
    Error handling features:
    - Retry logic with exponential backoff for transient errors
    - Automatic batch splitting for URL length errors
    - Logging of failed accessions with download URLs
    - Graceful degradation for partial failures
    
    Args:
        accessions (list): List of accession numbers for this batch
        failed_log_path (str, optional): Path to log file for failed batches
        
    Returns:
        tuple: (metadata_dict, xml_text)
            - metadata_dict: Dictionary mapping accessions to parsed metadata
            - xml_text: Raw XML content from the response
            
    Raises:
        RuntimeError: If the E-utilities request fails or XML parsing fails
        
    Note:
        - Implements retry logic for network errors (max 3 attempts)
        - Uses exponential backoff between retries (1s, 2s, 4s)
        - Automatically splits batches that are too large
        - Logs failed batches with individual download URLs for debugging
        
    Example:
        >>> batch = ['NC_045512.2', 'MN908947.3', 'MT020781.1']
        >>> metadata, xml = _fetch_genbank_batch(batch, 'failed_batches.log')
        >>> print(f"Retrieved metadata for {len(metadata)} accessions")
    """
    
    # Build E-utilities efetch URL for GenBank XML format
    # base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    accession_string = ",".join(accessions)
    
    params = {
        'db': 'nucleotide',           # Nucleotide database (includes virus sequences)
        'id': accession_string,       # Comma-separated accession numbers
        'rettype': 'gb',              # GenBank format
        'retmode': 'xml',             # XML output for structured parsing
        'complexity': 0,              # Requesting the entire blob of data content to be returned
    }
    
    # Create a requests.Session with urllib3 Retry/HTTPAdapter for robust retries
    session = requests.Session()
    try:
        retry_strategy = Retry(
            total=GENBANK_RETRY_ATTEMPTS,
            backoff_factor=HTTP_INITIAL_BACKOFF,
            status_forcelist=HTTP_RETRY_STATUS_CODES,
            allowed_methods=frozenset(['GET', 'POST'])
        )
    except TypeError:
        # Fallback for older urllib3 versions that use method_whitelist
        retry_strategy = Retry(
            total=GENBANK_RETRY_ATTEMPTS,
            backoff_factor=HTTP_INITIAL_BACKOFF,
            status_forcelist=HTTP_RETRY_STATUS_CODES,
            method_whitelist=frozenset(['GET', 'POST'])
        )

    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    headers = {'Connection': 'close', 'User-Agent': 'gget/1.0'}
    # Local retry loop for transient chunk/connection errors with exponential backoff
    max_attempts = HTTP_MAX_LOCAL_RETRIES
    attempt = 0
    backoff = HTTP_INITIAL_BACKOFF
    efetch_url = None  # Initialize to track the URL for logging

    while attempt < max_attempts:
        try:
            logger.debug("Making E-utilities request for %d accessions (attempt %d)", len(accessions), attempt + 1)
            logger.debug("Request URL: %s", NCBI_EUTILS_BASE)
            logger.debug("Request parameters: %s", {k: (v[:50] + '...' if isinstance(v, str) and len(v) > 50 else v) for k, v in params.items()})

            response = session.get(NCBI_EUTILS_BASE, params=params, timeout=EUTILS_TIMEOUT, headers=headers)
            efetch_url = response.url  # Capture the full URL for logging
            logger.debug("explicit URL requested: %s", response.url)
            response.raise_for_status()

            # Verify we got XML data
            if not response.text.strip().startswith('<?xml') and not response.text.strip().startswith('<'):
                raise RuntimeError(f"Invalid XML response: {response.text[:200]}")

            logger.debug("Received XML response: %d characters", len(response.text))

            # Parse the GenBank XML and extract metadata
            metadata_dict = _parse_genbank_xml(response.text)

            logger.debug("✅ Successfully parsed metadata for %d records", len(metadata_dict))
            return metadata_dict, response.text

        except (requests.exceptions.ChunkedEncodingError,
                urllib3.exceptions.ProtocolError,
                http.client.IncompleteRead,
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout) as e:
            attempt += 1
            logger.warning("Transient network error during E-utilities request (attempt %d/%d): %s", attempt, max_attempts, e)
            if attempt < max_attempts:
                logger.debug("Backing off for %.1f seconds before retrying", backoff)
                time.sleep(backoff)
                backoff *= 2
                continue
            else:
                logger.warning("Maximum local retry attempts reached for this batch (%d). Will attempt to split batch if possible.", max_attempts)
                break

        except requests.exceptions.RequestException as e:
            # Non-recoverable requests error (HTTP error status etc.)
            logger.error("❌ E-utilities request failed: %s", e)
            if hasattr(e, 'response') and e.response is not None:
                efetch_url = e.response.url
            break  # Exit retry loop

        except Exception as e:
            logger.error("❌ GenBank XML parsing failed: %s", e)
            break  # Exit retry loop

    # If local retries failed, try splitting into smaller sub-batches
    if len(accessions) > 1:
        logger.info("Attempting to split failing batch of %d accessions into smaller sub-batches", len(accessions))
        mid = len(accessions) // 2
        left_acc = accessions[:mid]
        right_acc = accessions[mid:]

        try:
            left_meta, left_xml = _fetch_genbank_batch(left_acc, failed_log_path=failed_log_path)
            right_meta, right_xml = _fetch_genbank_batch(right_acc, failed_log_path=failed_log_path)

            combined_meta = {}
            combined_meta.update(left_meta if left_meta else {})
            combined_meta.update(right_meta if right_meta else {})

            combined_xml = (left_xml or "") + "\n" + (right_xml or "")
            return combined_meta, combined_xml
        except Exception as split_error:
            logger.error("Failed to split batch: %s", split_error)

    # Cannot split further, log failed batch
    if failed_log_path:
        # Build the URL for manual download
        if efetch_url is None:
            # Construct URL if we never got a response
            base_url = f"{NCBI_EUTILS_BASE}?db=nucleotide&id={accession_string}&rettype=gb&retmode=xml"
            efetch_url = base_url
        
        with open(failed_log_path, "a") as flog:
            flog.write(f"FAILED_BATCH: {accessions}\n")
            flog.write(f"URL: {efetch_url}\n")
            flog.write(f"Individual URLs:\n")
            for acc in accessions:
                individual_url = f"{NCBI_EUTILS_BASE}?db=nucleotide&id={acc}&rettype=gb&retmode=xml"
                flog.write(f"  {acc}: {individual_url}\n")
            flog.write("\n")
        logger.error(f"Failed to fetch GenBank metadata for batch: {accessions[:5]}... Logged to {failed_log_path}")
    
    return {}, None


def _clean_xml_declarations(xml_text):
    """
    Remove XML declarations and DOCTYPE declarations from XML text.
    
    When concatenating multiple XML documents, we need to remove the XML
    declaration (<?xml...?>) and DOCTYPE declarations from each document
    to create valid combined XML.
    
    Args:
        xml_text (str): Raw XML text that may contain declarations
        
    Returns:
        str: Cleaned XML text without declarations
        
    Example:
        >>> xml = '<?xml version="1.0"?>\n<!DOCTYPE GBSet>\n<GBSet>...</GBSet>'
        >>> _clean_xml_declarations(xml)
        '<GBSet>...</GBSet>'
    """
    cleaned_lines = []
    for line in xml_text.splitlines():
        # Skip XML declarations and DOCTYPE declarations
        if line.strip().startswith("<?xml") or line.strip().startswith("<!DOCTYPE"):
            continue
        cleaned_lines.append(line)
    return "\n".join(cleaned_lines)


def _local_name(tag):
    """
    Return the local name of an XML tag (strip namespace if present).
    
    XML tags may include namespace prefixes (e.g., '{http://namespace}TagName').
    This helper function extracts just the tag name without the namespace.
    
    Args:
        tag (str): XML tag string, potentially with namespace
        
    Returns:
        str: Tag name without namespace prefix
        
    Example:
        >>> _local_name('{http://www.ncbi.nlm.nih.gov}GBSeq')
        'GBSeq'
        >>> _local_name('GBSeq')
        'GBSeq'
    """
    return tag.split('}')[-1] if '}' in tag else tag


def _genbank_xml_to_csv(xml_path, csv_path, chunk_size=None):
    """
    Convert GenBank XML to structured CSV with streaming and dynamic qualifier columns.
    
    This function parses a GenBank XML file and extracts feature data into a flat CSV
    format suitable for analysis. It uses streaming parsing to handle large files
    efficiently and dynamically discovers all qualifier columns present in the data.
    
    The output CSV includes:
    - Basic sequence information (accession, sequence)
    - Feature information (key, location, intervals)
    - All qualifiers as separate columns (discovered dynamically)
    
    Args:
        xml_path (str): Path to input GenBank XML file
        csv_path (str): Path to output CSV file
        chunk_size (int): Number of rows to process before writing to disk (default: 10000)
                         Larger values use more memory but may be faster
    
    Note:
        - Uses streaming to handle files larger than available memory
        - Dynamically discovers qualifier columns from the XML data
        - The sequence is only added to the last row for each accession to save space
        - Writes incrementally in chunks to manage memory usage
        
    Example:
        >>> _genbank_xml_to_csv('genbank_data.xml', 'output.csv', chunk_size=5000)
    """
    # Apply default chunk size if not specified
    if chunk_size is None:
        chunk_size = GENBANK_XML_CHUNK_SIZE
    
    qualifier_names = set()
    rows = []
    header_written = False

    csv_file = open(csv_path, "w", newline='', encoding='utf-8')
    writer = None

    # Stream-parse XML
    for event, elem in ET.iterparse(xml_path, events=("end",)):
        if _local_name(elem.tag) == "GBSeq":
            accession = elem.findtext(".//GBSeq_accession-version", "").strip()
            sequence = elem.findtext(".//GBSeq_sequence", "").strip()
            features = elem.findall(".//GBFeature")

            seq_rows = []  # rows for this GBSeq

            for feature in features:
                feature_key = feature.findtext("GBFeature_key", "").strip()
                feature_location = feature.findtext("GBFeature_location", "").strip()
                feature_operator = feature.findtext("GBFeature_operator", "").strip()

                intervals = feature.findall(".//GBInterval")
                quals = feature.findall(".//GBQualifier")

                qualifiers = {}
                for q in quals:
                    name = q.findtext("GBQualifier_name", "").strip()
                    value = q.findtext("GBQualifier_value", "").strip()
                    if name:
                        qualifiers[name] = value
                        qualifier_names.add(name)

                if not intervals:
                    intervals = [None]

                for interval in intervals:
                    interval_from = (interval.findtext("GBInterval_from") if interval is not None else "") or ""
                    interval_to = (interval.findtext("GBInterval_to") if interval is not None else "") or ""
                    interval_acc = (interval.findtext("GBInterval_accession") if interval is not None else "") or ""

                    interval_from = interval_from.strip()
                    interval_to = interval_to.strip()
                    interval_acc = interval_acc.strip()

                    order_flag = feature_operator if feature_operator.startswith("order(") else ""

                    row = {
                        "accession": accession,
                        "Feature_key": feature_key,
                        "Feature_location": feature_location,
                        "Interval_from": interval_from,
                        "Interval_to": interval_to,
                        "Interval_accession": interval_acc,
                        "order": order_flag,
                        "sequence": ""  # leave blank for now
                    }

                    for qn in qualifier_names:
                        row[qn] = qualifiers.get(qn, "")

                    seq_rows.append(row)

            # Add sequence only to the last row for this GBSeq
            if seq_rows:
                seq_rows[-1]["sequence"] = sequence

            rows.extend(seq_rows)

            # Write chunk
            if len(rows) >= chunk_size:
                df = pd.DataFrame(rows)
                if not header_written:
                    df.to_csv(csv_file, index=False)
                    header_written = True
                else:
                    df.to_csv(csv_file, index=False, header=False)
                rows.clear()

            elem.clear()

    # Write remaining rows
    if rows:
        df = pd.DataFrame(rows)
        if not header_written:
            df.to_csv(csv_file, index=False)
        else:
            df.to_csv(csv_file, index=False, header=False)

    csv_file.close()
    logger.debug(f"✅ Finished writing {csv_path}")


def _save_genbank_xml_and_csv(xml_content, xml_file_name, csv_file_name):
    """
    Save GenBank XML content and flattened CSV representation.
    
    This function takes raw GenBank XML data and saves it in two formats:
    1. Raw XML file - preserves all original data
    2. Flattened CSV file - extracts key fields into tabular format
    
    The CSV conversion makes the data more accessible for analysis tools
    that work with tabular data, while the XML preserves all information
    in its original structure.
    
    Args:
        xml_content (str): Raw XML content from E-utilities efetch
        xml_file_name (str): Path where XML file should be saved
        csv_file_name (str): Path where CSV file should be saved
        
    Raises:
        RuntimeError: If XML parsing fails or file writing encounters errors
        ET.ParseError: If XML content is invalid or malformed
        
    Example:
        >>> xml = fetch_genbank_xml(['NC_045512.2'])
        >>> _save_genbank_xml_and_csv(xml, 'data.xml', 'data.csv')
    """
    try:
        root = ET.fromstring(xml_content)
        logger.debug("✅ XML parsing successful, root element: %s", root.tag)

        logger.debug("Saving GenBank XML content to file: %s", xml_file_name)
        with open(xml_file_name, "w", encoding="utf-8") as f:
            f.write(xml_content)
        logger.debug("✅ Saved GenBank XML content to: %s", xml_file_name)

        logger.debug("Saving the loss-less genbank information into a CSV if the flag is set.")
        logger.debug("Reformatting to extract important data from GenBank XML to save as CSV: %s", csv_file_name)
        _genbank_xml_to_csv(xml_file_name, csv_file_name)
        logger.debug("✅ Saved GenBank data in a loss-less csv file to: %s", csv_file_name)

    except ET.ParseError as e:
        logger.error("❌ XML parsing failed: %s", e)
        logger.debug("XML content preview: %s", xml_content[:500])
        raise RuntimeError(f"Invalid XML format in GenBank response: {e}") from e


def _parse_genbank_xml(xml_content):
    """
    Parse GenBank XML response and extract high-level metadata fields.
    This function processes the GenBank XML format returned by E-utilities efetch
    and extracts key metadata including collection dates, geographic information,
    host details, publication references, and sequence features.
    
    Args:
        xml_content (str): Raw XML content from E-utilities efetch
        
    Returns:
        dict: Dictionary mapping accession numbers to metadata dictionaries
        
    Note:
        Uses xml.etree.ElementTree for parsing to avoid external dependencies.
        The GenBank XML schema is documented by NCBI and contains structured
        information about sequence records.
    """
    
    # Parse the XML content
    try:
        root = ET.fromstring(xml_content)
        logger.debug("✅ XML parsing successful, root element: %s", root.tag)

    except ET.ParseError as e:
        logger.error("❌ XML parsing failed: %s", e)
        logger.debug("XML content preview: %s", xml_content[:500])
        raise RuntimeError(f"Invalid XML format in GenBank response: {e}") from e
    
    metadata_dict = {}
    
    # Process each GenBank sequence record in the XML
    for gbseq in root.findall('.//GBSeq'):
        try:
            # Extract accession number as the primary key
            accession_elem = gbseq.find('GBSeq_accession-version')
            if accession_elem is None:
                accession_elem = gbseq.find('GBSeq_primary-accession')
            
            if accession_elem is None:
                logger.warning("Skipping GenBank record without accession number")
                continue
                
            accession = accession_elem.text
            logger.debug("Processing GenBank record: %s", accession)
            
            # Initialize metadata dictionary for this record
            metadata = {
                'accession': accession,
                'genbank_data': {}  # Store GenBank-specific fields
            }
            
            # Extract basic sequence information
            length_elem = gbseq.find('GBSeq_length')
            metadata['genbank_data']['sequence_length'] = int(length_elem.text) if length_elem is not None else None
            
            organism_elem = gbseq.find('GBSeq_organism')
            metadata['genbank_data']['organism'] = organism_elem.text if organism_elem is not None else ""
            
            definition_elem = gbseq.find('GBSeq_definition')
            metadata['genbank_data']['definition'] = definition_elem.text if definition_elem is not None else ""
            
            # Extract taxonomy information
            taxonomy_elem = gbseq.find('GBSeq_taxonomy')
            metadata['genbank_data']['taxonomy'] = taxonomy_elem.text if taxonomy_elem is not None else ""
            
            # Extract creation and update dates
            create_date_elem = gbseq.find('GBSeq_create-date')
            metadata['genbank_data']['create_date'] = create_date_elem.text if create_date_elem is not None else ""
            
            update_date_elem = gbseq.find('GBSeq_update-date')
            metadata['genbank_data']['update_date'] = update_date_elem.text if update_date_elem is not None else ""
            
            # Extract references (publications)
            references = []
            for ref in gbseq.findall('.//GBReference'):
                ref_data = {}
                
                title_elem = ref.find('GBReference_title')
                ref_data['title'] = title_elem.text if title_elem is not None else ""
                
                authors_elem = ref.find('GBReference_authors')
                if authors_elem is not None:
                    authors = [a.text for a in authors_elem.findall('GBAuthor') if a.text]
                    ref_data['authors'] = ', '.join(authors)
                else:
                    ref_data['authors'] = ""
                # ref_data['authors'] = authors_elem.text if authors_elem is not None else ""
                
                journal_elem = ref.find('GBReference_journal')
                ref_data['journal'] = journal_elem.text if journal_elem is not None else ""
                
                pubmed_elem = ref.find('GBReference_pubmed')
                ref_data['pubmed_id'] = pubmed_elem.text if pubmed_elem is not None else ""
                
                if any(ref_data.values()):  # Only add if we got some reference data
                    references.append(ref_data)
            
            metadata['genbank_data']['references'] = references
            
            # Extract features (collection_date, geographic location, host, etc.)
            features_data = {}
            
            for feature in gbseq.findall('.//GBFeature'):
                feature_key_elem = feature.find('GBFeature_key')
                if feature_key_elem is None:
                    continue
                    
                feature_key = feature_key_elem.text
                
                # Extract qualifiers for this feature
                feature_qualifiers = {}
                for qual in feature.findall('.//GBQualifier'):
                    qual_name_elem = qual.find('GBQualifier_name')
                    qual_value_elem = qual.find('GBQualifier_value')
                    
                    if qual_name_elem is not None and qual_value_elem is not None:
                        qual_name = qual_name_elem.text
                        qual_value = qual_value_elem.text
                        feature_qualifiers[qual_name] = qual_value
                
                if feature_qualifiers:
                    features_data[feature_key] = feature_qualifiers
            
            # Extract specific fields of interest from source feature
            source_feature = features_data.get('source', {})
            
            metadata['genbank_data']['collection_date'] = source_feature.get('collection_date', '')
            metadata['genbank_data']['geographic_location'] = source_feature.get('geo_loc_name', '')
            metadata['genbank_data']['host'] = source_feature.get('host', '')
            metadata['genbank_data']['isolation_source'] = source_feature.get('isolation_source', '')
            metadata['genbank_data']['strain'] = source_feature.get('strain', '')
            metadata['genbank_data']['isolate'] = source_feature.get('isolate', '')
            # metadata['genbank_data']['collected_by'] = source_feature.get('collected_by', '')
            # metadata['genbank_data']['specimen_voucher'] = source_feature.get('specimen_voucher', '')
            
            # Store all features for potential future use
            metadata['genbank_data']['all_features'] = features_data
            
            # Extract comment field (often contains additional metadata)
            comment_elem = gbseq.find('GBSeq_comment')
            comment_text = comment_elem.text if comment_elem is not None else ""
            metadata['genbank_data']['comment'] = comment_text
            
            # Parse assembly name from comment if present (used in some studies)
            assembly_name = ""
            if comment_text:
                assembly_match = re.search(r'Assembly Name :: (\S+)', comment_text)
                if assembly_match:
                    assembly_name = assembly_match.group(1)
            metadata['genbank_data']['assembly_name'] = assembly_name
            
            # Store the metadata for this accession
            metadata_dict[accession] = metadata
            
            logger.debug("Extracted GenBank metadata for %s: organism=%s, collection_date=%s, geographic_location=%s", 
                        accession, 
                        metadata['genbank_data']['organism'],
                        metadata['genbank_data']['collection_date'],
                        metadata['genbank_data']['geographic_location'])
            
        except Exception as e:
            logger.warning("❌ Failed to parse GenBank record %s: %s", 
                          accession if 'accession' in locals() else 'unknown', e)
            continue
    
    logger.info("✅ Successfully parsed GenBank metadata for %d records", len(metadata_dict))
    return metadata_dict


def save_genbank_metadata_to_csv(genbank_metadata, output_file, virus_metadata=None):
    """
    Save GenBank metadata to a human-readable CSV file.
    
    This function creates a comprehensive CSV file containing GenBank-specific metadata
    that complements the standard virus metadata. The output includes collection dates,
    geographic information, host details, publication references, and other fields
    extracted from GenBank records.
    
    Args:
        genbank_metadata (dict): Dictionary mapping accessions to GenBank metadata
        output_file (str): Path to the output CSV file
        virus_metadata (list, optional): List of virus metadata dictionaries to merge
        
    Note:
        The CSV format is designed to be easily readable and compatible with
        downstream analysis tools. Complex nested data (like references) is
        flattened into separate columns or JSON-encoded strings.
    """
    
    logger.info("Preparing GenBank metadata for CSV output...")
    logger.debug("Processing %d GenBank records", len(genbank_metadata))
    
    # Define the column order for the GenBank metadata CSV
    # Prioritize the most commonly used and important fields
    columns = [
        # Basic identifiers
        "accession",
        "organism",
        "definition",
        "sequence_length",
        
        # Collection and geographic information
        "collection_date",
        "geographic_location",
        "strain",
        "isolate",
        
        # Host and source information
        "host",
        "isolation_source",
        
        # Database and version information
        "create_date",
        "update_date",
        "assembly_name",
        
        # Taxonomic information
        "taxonomy",
        
        # Publication information
        "authors",
        "title",
        "journal",
        "pubmed_id",
        "reference_count",
        
        # Additional metadata
        "comment",
    ]
    
    logger.debug("Using column order: %s", columns)
    
    # Prepare data for DataFrame creation
    data_for_df = []
    
    for accession, metadata in genbank_metadata.items():
        logger.debug("Processing GenBank metadata for accession: %s", accession)
        
        genbank_data = metadata.get('genbank_data', {})
        
        # Extract publication information (use first reference if available)
        references = genbank_data.get('references', [])
        first_ref = references[0] if references else {}
        
        # Build the row dictionary
        row = {
            "accession": accession,
            "organism": genbank_data.get('organism', ''),
            "definition": genbank_data.get('definition', ''),
            "sequence_length": genbank_data.get('sequence_length', ''),
            
            "collection_date": genbank_data.get('collection_date', ''),
            # "country": genbank_data.get('country', ''),
            "geographic_location": genbank_data.get('geographic_location', ''),
            "strain": genbank_data.get('strain', ''),
            "isolate": genbank_data.get('isolate', ''),
            
            "host": genbank_data.get('host', ''),
            "isolation_source": genbank_data.get('isolation_source', ''),
            # "collected_by": genbank_data.get('collected_by', ''),
            # "specimen_voucher": genbank_data.get('specimen_voucher', ''),
            
            "create_date": genbank_data.get('create_date', ''),
            "update_date": genbank_data.get('update_date', ''),
            "assembly_name": genbank_data.get('assembly_name', ''),
            
            "taxonomy": genbank_data.get('taxonomy', ''),
            
            "authors": first_ref.get('authors', ''),
            "title": first_ref.get('title', ''),
            "journal": first_ref.get('journal', ''),
            "pubmed_id": first_ref.get('pubmed_id', ''),
            "reference_count": len(references),
            
            "comment": genbank_data.get('comment', ''),
        }
        
        data_for_df.append(row)
    
    logger.info("Creating DataFrame with %d rows and %d columns", len(data_for_df), len(columns))
    
    # Create DataFrame with the specified column order
    df = pd.DataFrame(data_for_df, columns=columns)
    
    # If virus metadata is provided, try to merge it
    if virus_metadata:
        logger.info("Merging with virus metadata (%d records)", len(virus_metadata))
        try:
            # Create virus metadata DataFrame
            virus_df = pd.DataFrame(virus_metadata)
            if 'accession' in virus_df.columns:
                # Merge on accession number
                df = pd.merge(df, virus_df, on='accession', how='outer', suffixes=('_genbank', '_virus'))
                logger.info("✅ Successfully merged GenBank and virus metadata")
            else:
                logger.warning("Cannot merge: virus metadata missing 'accession' column")
        except Exception as e:
            logger.warning("❌ Failed to merge GenBank and virus metadata: %s", e)

    # Write DataFrame to CSV file
    try:
        df.to_csv(output_file, index=False, encoding='utf-8')
        logger.info("✅ Merged GenBank metadata successfully saved to: %s", output_file)
        logger.info("CSV file contains %d rows and %d columns", len(df), len(df.columns)) 
    except Exception as e:
        logger.error("❌ Failed to save GenBank metadata CSV: %s", e)
        raise RuntimeError(f"❌ Failed to save GenBank metadata to {output_file}: {e}") from e
    

def filter_metadata_only(
    metadata_dict,
    min_seq_length=None,
    max_seq_length=None,
    min_gene_count=None,
    max_gene_count=None,
    nuc_completeness=None,
    lab_passaged=None,
    submitter_country=None,
    min_collection_date=None,
    max_collection_date=None,
    source_database=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
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
                "completeness(%s), lab_passaged(%s), "
                "submitter_country(%s), collection_date(%s-%s), source_db(%s), max_release_date(%s), "
                "peptide_count(%s-%s), protein_count(%s-%s)",
                min_seq_length, max_seq_length, min_gene_count, max_gene_count,
                nuc_completeness, lab_passaged,
                submitter_country, min_collection_date, max_collection_date, source_database, max_release_date, 
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
        'lab_passaged': 0,
        'submitter_country': 0,
        'collection_date': 0,
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
        if nuc_completeness is not None and not nuc_completeness.lower() == "complete":
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

        # FILTER 4: Lab passaging status filter
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

        # FILTER 5: Submitter country filter
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

        # FILTER 6: Collection date range filter
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

        # FILTER 7: Source database filter
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

        # FILTER 8: Maximum release date filter
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

        # FILTER 9: Mature peptide count filters
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

        # FILTER 10: Protein count filters
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
    logger.info("✅ Metadata-only filtering process complete.")
    logger.info("=================================")
    logger.info("Metadata-only filtering complete:")
    logger.info("  Total metadata records: %d", total_sequences)
    logger.info("  Records passing filters: %d", num_filtered)

    if num_filtered == 0:
        logger.info("No records passed the metadata-only filters.")

    # Log detailed filter statistics if any records were filtered out
    total_filtered = sum(filter_stats.values())
    if total_filtered > 0:
        logger.info("Filter statistics (records excluded):")
        for filter_name, count in filter_stats.items():
            if count > 0:
                logger.info("  %s: %d records", filter_name, count)
    
    return filtered_accessions, filtered_metadata_list


def virus(
    virus,
    is_accession=False,
    outfolder=None,
    host=None,
    min_seq_length=None,
    max_seq_length=None,
    min_gene_count=None,
    max_gene_count=None,
    nuc_completeness=None,
    has_proteins=None,
    proteins_complete=False,
    lab_passaged=None,
    geographic_location=None,
    submitter_country=None,
    min_collection_date=None,
    max_collection_date=None,
    annotated=None,
    refseq_only=False,
    keep_temp=False,
    source_database=None,
    min_release_date=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
    max_ambiguous_chars=None,
    is_sars_cov2=False,
    is_alphainfluenza=False,
    lineage=None,
    genbank_metadata=False,
    genbank_batch_size=200,
    download_all_accessions=False,
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
    6) Optionally fetch detailed GenBank metadata and save to CSV
    
    Args:
        virus (str): Virus taxon name/ID or accession number
        is_accession (bool): Whether virus parameter is an accession number
        outfolder (str): Output directory for files
        host (str): Host organism name filter
        min_seq_length (int): Minimum sequence length filter
        max_seq_length (int): Maximum sequence length filter
        min_gene_count (int): Minimum gene count filter
        max_gene_count (int): Maximum gene count filter
        nuc_completeness (str): Nucleotide completeness filter ('complete' or 'partial')
        has_proteins (str/list): Required proteins/genes filter
        proteins_complete (bool): Whether proteins must be complete
        lab_passaged (bool): Lab passaging status filter
        geographic_location (str): Geographic location filter
        submitter_country (str): Submitter country filter
        min_collection_date (str): Minimum collection date filter (YYYY-MM-DD)
        max_collection_date (str): Maximum collection date filter (YYYY-MM-DD)
        annotated (bool): Annotation status filter
        source_database (str): Source database filter
        min_release_date (str): Minimum release date filter (YYYY-MM-DD)
        max_release_date (str): Maximum release date filter (YYYY-MM-DD)
        min_mature_peptide_count (int): Minimum mature peptide count filter
        max_mature_peptide_count (int): Maximum mature peptide count filter
        min_protein_count (int): Minimum protein count filter
        max_protein_count (int): Maximum protein count filter
        max_ambiguous_chars (int): Maximum ambiguous nucleotide character filter
        is_sars_cov2 (bool): Flag to indicate if the accession is for SARS-CoV-2, enabling optimized download method (default: False)
        is_alphainfluenza (bool): Flag to indicate if the query is for Alphainfluenza, enabling optimized download method (default: False)
        lineage (str): Virus lineage filter (SARS-CoV-2 specific)
        genbank_metadata (bool): Whether to fetch detailed GenBank metadata (default: False)
        genbank_batch_size (int): Batch size for GenBank API requests (default: 200)
        keep_temp (bool): Flag to indicate if all output files should be saved, including intermediate files (default: False)
        refseq_only (bool): Whether to restrict to RefSeq sequences only

    Returns:
        None: Files are saved to the output directory
        
    Note:
        When genbank_metadata=True, an additional CSV file with detailed GenBank
        metadata will be saved alongside the standard output files. This includes
        collection dates, geographic information, host details, publication
        references, and other fields extracted from GenBank records.
    """
    logger.info("Starting NCBI virus data retrieval process...")
    
    # Capture the command line for summary
    command_line = " ".join(sys.argv) if len(sys.argv) > 0 else "virus (called programmatically)"
    
    # Initialize variables for tracking results
    total_api_records = 0
    total_after_metadata_filter = 0
    total_final_sequences = 0
    output_files_dict = {}
    final_metadata_for_summary = []
    
    # Initialize failed commands tracker
    failed_commands = {
        'api_timeout': None,
        'sequence_batches': [],
        'genbank_batches': []
    }
    
    # Track if GenBank metadata was successfully retrieved
    genbank_success = False
    genbank_error_msg = None

    if download_all_accessions:
        logger.info("ATTENTION: Download all accessions mode is active.")
        logger.info("This will download all virus accessions from NCBI, which can be a very large dataset and take a long time.")
        virus = NCBI_ALL_VIRUSES_TAXID  # NCBI taxonomy ID for all Viruses
        is_accession = False
        logger.info("Overriding virus query to fetch all viruses using taxon ID: %s", virus)

    logger.info("Query parameters: virus='%s', is_accession=%s, outfolder='%s'", 
                virus, is_accession, outfolder)
    logger.debug("Applied filters: host=%s, seq_length=(%s-%s), gene_count=(%s-%s), completeness=%s, annotated=%s, refseq_only=%s, keep_temp=%s, lab_passaged=%s, geo_location=%s, submitter_country=%s, collection_date=(%s-%s), source_db=%s, release_date=(%s-%s), protein_count=(%s-%s), peptide_count=(%s-%s), max_ambiguous=%s, has_proteins=%s, proteins_complete=%s, genbank_metadata=%s, genbank_batch_size=%s",
    host, min_seq_length, max_seq_length, min_gene_count, max_gene_count, nuc_completeness, annotated, refseq_only, keep_temp, lab_passaged, geographic_location, submitter_country, min_collection_date, max_collection_date,source_database, min_release_date, max_release_date, min_protein_count, max_protein_count, min_mature_peptide_count, max_mature_peptide_count, max_ambiguous_chars, has_proteins, proteins_complete, genbank_metadata, genbank_batch_size)

    # SECTION 1: INPUT VALIDATION AND OUTPUT DIRECTORY SETUP
    # Validate and normalize input arguments before proceeding
    logger.info("=" * 60)
    logger.info("STEP 1: VALIDATING INPUT ARGUMENTS AND OUTPUT DIRECTORY SETUP...")
    logger.info("=" * 60)
    
    # Validate virus parameter
    if virus is None or (isinstance(virus, str) and virus.strip() == ""):
        raise ValueError(
            "Argument 'virus' must be a non-empty string (virus name, taxon ID, or accession number)."
        )
    
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

    if refseq_only is not None and not isinstance(refseq_only, bool):
        raise TypeError(
            "Argument 'refseq_only' must be a boolean (True or False)."
        )
    
    if keep_temp is not None and not isinstance(keep_temp, bool):
        raise TypeError(
            "Argument 'keep_temp' must be a boolean (True or False)."
        )

    if is_accession is not None and not isinstance(is_accession, bool):
        accession = is_accession
        raise TypeError(
            "Argument 'is_accession' must be a boolean (True or False)."
        )
    
    # Validate GenBank metadata parameters
    if genbank_metadata is not None and not isinstance(genbank_metadata, bool):
        raise TypeError(
            "Argument 'genbank_metadata' must be a boolean (True or False)."
        )
    
    if genbank_batch_size is not None:
        if not isinstance(genbank_batch_size, int) or genbank_batch_size <= 0:
            raise ValueError(
                "Argument 'genbank_batch_size' must be a positive integer."
            )
        if genbank_batch_size > GENBANK_MAX_BATCH_SIZE_WARNING:
            logger.warning("Large genbank_batch_size (%d) may cause API timeouts. Consider using smaller batches.", genbank_batch_size)
    
    # Log GenBank metadata configuration
    if genbank_metadata:
        logger.info("GenBank metadata retrieval enabled (batch_size=%d)", genbank_batch_size)
    else:
        logger.debug("GenBank metadata retrieval disabled")

    
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
        min_collection_date,
        max_collection_date,
        "collection date (arguments: min_collection_date and max_collection_date)",
        date=True,  # Enable date parsing for comparison
    )
    check_min_max(
        min_release_date,
        max_release_date,
        "release data (arguments: min_release_date and max_release_date)",
        date=True,  # Enable date parsing for comparison
    )

    logger.info("Input validation completed successfully")

    virus_clean = virus.replace(' ', '_').replace('/', '_')
    # Create and prepare output directory structure
    if outfolder is None:
        currentfolder = os.getcwd()
        outfolder = f"{currentfolder}/{virus_clean}_{timestamp}"
        logger.info("No output folder specified, creating a subdirectory in current directory: %s", outfolder)
    else:
        logger.info("Using specified output folder: %s", outfolder)
    
    # Ensure output folder exists
    os.makedirs(outfolder, exist_ok=True)
    logger.debug("Output folder ready: %s", outfolder)
    
    # SECTION 2: SARS-CoV-2 CACHED DATA PROCESSING
    # For SARS-CoV-2 queries, use cached data packages with hierarchical fallback
    logger.info("=" * 60)
    logger.info("STEP 2: CHECKING FOR SARS-CoV-2 QUERY...")
    logger.info("=" * 60)

    # Initialize variables to track cached download results
    cached_sequences = None
    cached_metadata_dict = None
    used_cached_download = False

    if is_sars_cov2 or is_sars_cov2_query(virus, is_accession):
        logger.info("DETECTED SARS-CoV-2 QUERY - USING CACHED DATA PACKAGES")
        logger.info("SARS-CoV-2 queries will use NCBI's optimized cached data packages")
        logger.info("with hierarchical fallback from specific to general cached files.")
        
        # Use the download_sars_cov2_optimized function which handles fallback strategies internally
        params = {
            'host': host,
            'complete_only': (nuc_completeness == "complete"),
            'annotated': annotated,
            'outdir': outfolder,
            'lineage': lineage,
            'accession': virus,
            'use_accession': is_accession
        }
            
        zip_file = download_sars_cov2_optimized(**params)
        
        try:
            cached_sequences, cached_metadata_dict, used_cached_download = _process_cached_download(
                zip_file, virus_type="SARS-CoV-2"
            )
            if used_cached_download:
                logger.info("Cached download completed. Server-side filters (host, complete_only, annotated, lineage) applied.")
                logger.info("All other filters will be applied in the unified filtering pipeline.")
        except Exception as cache_error:
            logger.warning("❌ Cached data processing failed: %s", cache_error)
            logger.info("Proceeding with regular API workflow...")
    
    # SECTION 2b: ALPHAINFLUENZA CACHED DATA PROCESSING
    # For Alphainfluenza queries, use cached data packages with hierarchical fallback
    logger.info("=" * 60)
    logger.info("STEP 2b: CHECKING FOR ALPHAINFLUENZA QUERY...")
    logger.info("=" * 60)

    if is_alphainfluenza or is_alphainfluenza_query(virus, is_accession):
        logger.info("DETECTED ALPHAINFLUENZA QUERY - USING CACHED DATA PACKAGES")
        logger.info("Alphainfluenza queries will use NCBI's optimized cached data packages")
        logger.info("with hierarchical fallback from specific to general cached files.")
        
        # Use the download_alphainfluenza_optimized function which handles fallback strategies internally
        params = {
            'host': host,
            'complete_only': (nuc_completeness == "complete"),
            'annotated': annotated,
            'outdir': outfolder,
            'accession': virus,
            'use_accession': is_accession
        }
            
        zip_file = download_alphainfluenza_optimized(**params)
        
        try:
            cached_sequences, cached_metadata_dict, used_cached_download = _process_cached_download(
                zip_file, virus_type="Alphainfluenza"
            )
            if used_cached_download:
                logger.info("Cached download completed. Server-side filters (host, complete_only, annotated) applied.")
                logger.info("All other filters will be applied in the unified filtering pipeline.")
        except Exception as cache_error:
            logger.warning("❌ Cached data processing failed: %s", cache_error)
            logger.info("Proceeding with regular API workflow...")
    
    # Create temporary directory for intermediate processing
    # This will be cleaned up at the end regardless of success or failure
    temp_dir = os.path.join(outfolder, f"tmp_{timestamp}_{random_suffix}")
    os.makedirs(temp_dir, exist_ok=True)
    logger.debug("Created temporary processing directory: %s", temp_dir)

    # File names which will be referenced later
    genbank_csv_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata.csv")
    genbank_full_xml_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata_full.xml")
    genbank_full_csv_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata_full.csv")
    output_api_metadata_jsonl = os.path.join(outfolder, f"{virus_clean}_api_metadata.jsonl")

    try:
    # SECTION 3: METADATA RETRIEVAL WHILE APPLYING SERVER-SIDE FILTERS
        # logger.info("=" * 60)
        # logger.info("STEP 3: Fetching virus metadata from NCBI API")
        # logger.info("=" * 60)

        # Check if we're using cached download data
        if used_cached_download and cached_metadata_dict:
            logger.info("=" * 60)
            logger.info("STEP 3: Saving metadata from cached download")
            logger.info("=" * 60)
            logger.info("Using metadata from cached download (skipping API metadata fetch)")
            metadata_dict = cached_metadata_dict
            total_api_records = len(metadata_dict)
            logger.info("Loaded %d metadata records from cached download", total_api_records)
            
            # Save the cached metadata as API metadata for consistency
            logger.debug("Writing cached metadata to JSONL file: %s", output_api_metadata_jsonl)
            try:
                with open(output_api_metadata_jsonl, "w", encoding="utf-8") as f:
                    for md in metadata_dict.values():
                        f.write(json.dumps(md) + "\n")
                logger.info("✅ Saved cached metadata JSONL: %s", output_api_metadata_jsonl)
            except Exception as e:
                logger.warning("❌ Failed to save cached metadata JSONL: %s", e)
        else:
            # Regular API metadata fetch
            logger.info("=" * 60)
            logger.info("STEP 3: Fetching virus metadata from NCBI API")
            logger.info("=" * 60)
            api_annotated_filter = annotated if annotated is True else None
            api_complete_filter = True if nuc_completeness=="complete" else False

            logger.debug("Applying server-side filters: host=%s, geo_location=%s, annotated=%s, complete_only=%s, min_release_date=%s, refseq_only=%s", host, geographic_location, annotated, api_complete_filter, min_release_date, refseq_only)
            
            api_reports = fetch_virus_metadata(
                virus,
                accession=is_accession,
                host=host,
                geographic_location=geographic_location,
                annotated=api_annotated_filter,
                complete_only=api_complete_filter,
                min_release_date=min_release_date,
                refseq_only=refseq_only,
                failed_commands=failed_commands,
            )

            # If fetch_virus_metadata returns None, it means the dataset is too large
            # and we need to use the chunked download strategy
            if api_reports is None:
                logger.info("Standard download failed due to dataset size - switching to chunked download")
                api_reports = fetch_virus_metadata_chunked(
                    virus,
                    accession=is_accession,
                    host=host,
                    geographic_location=geographic_location,
                    annotated=api_annotated_filter,
                    complete_only=api_complete_filter,
                    min_release_date=min_release_date,
                    refseq_only=refseq_only,
                    failed_commands=failed_commands,
                )

            if not api_reports:
                logger.warning("No virus records found matching the specified criteria.")
                logger.info("Consider relaxing your filter criteria or checking your virus identifier.")
                return

            logger.info("Successfully retrieved %d virus records from API", len(api_reports))
            total_api_records = len(api_reports)

            # Convert API metadata to internal format
            logger.debug("Converting API metadata to internal format...")
            metadata_dict = load_metadata_from_api_reports(api_reports)
            logger.info("Processed metadata for %d sequences", len(metadata_dict))

            # Save the raw API metadata (server-side filtered) before local filtering
            logger.debug("Writing API metadata to JSONL file: %s", output_api_metadata_jsonl)
            try:
                with open(output_api_metadata_jsonl, "w", encoding="utf-8") as f:
                    for md in metadata_dict.values():
                        f.write(json.dumps(md) + "\n")
                logger.info("✅ Saved API metadata JSONL: %s", output_api_metadata_jsonl)
            except Exception as e:
                logger.warning("❌ Failed to save API metadata JSONL: %s", e)

    # SECTION 4: METADATA-ONLY FILTERING
        logger.info("=" * 60)
        logger.info("STEP 4: Applying metadata-only filters")
        logger.info("=" * 60)

        filters = {
            "min_seq_length": min_seq_length,
            "max_seq_length": max_seq_length,
            "min_gene_count": min_gene_count,
            "max_gene_count": max_gene_count,
            "nuc_completeness": nuc_completeness, #only for partial cases
            "lab_passaged": lab_passaged,
            "submitter_country": submitter_country,
            "min_collection_date": min_collection_date,
            "max_collection_date": max_collection_date,
            "source_database": source_database,
            "max_release_date": max_release_date,
            "min_mature_peptide_count": min_mature_peptide_count,
            "max_mature_peptide_count": max_mature_peptide_count,
            "min_protein_count": min_protein_count,
            "max_protein_count": max_protein_count,
        }

        all_metadata_filters_none_except_nuc = all(
            v is None for k, v in filters.items() if k != "nuc_completeness"
        )

        # Prepare output file paths (defined early for use in cleanup even if filters return early)
        output_fasta_file = os.path.join(outfolder, f"{virus_clean}_sequences.fasta")
        output_metadata_csv = os.path.join(outfolder, f"{virus_clean}_metadata.csv")
        output_metadata_jsonl = os.path.join(outfolder, f"{virus_clean}_metadata.jsonl")

        if all_metadata_filters_none_except_nuc and filters["nuc_completeness"]!="partial":
            logger.info("No metadata-only filters specified, skipping this step.")
            filtered_accessions = list(metadata_dict.keys())
            filtered_metadata = list(metadata_dict.values())
            logger.info("All %d sequences will proceed to sequence download", len(filtered_accessions))
        else:
            filtered_accessions, filtered_metadata = filter_metadata_only(metadata_dict, **filters)
            if not filtered_accessions:
                logger.warning("No sequences passed metadata-only filters. Skipping sequence download.")
                total_after_metadata_filter = 0
                total_final_sequences = 0
                # Save command summary even if no sequences passed filters
                save_command_summary(
                    outfolder=outfolder,
                    command_line=command_line,
                    total_api_records=total_api_records,
                    total_after_metadata_filter=0,
                    total_final_sequences=0,
                    output_files=output_files_dict,
                    filtered_metadata=[],
                    success=True,
                    error_message="No sequences passed the metadata filters",
                    failed_commands=failed_commands
                )
                return
        
        total_after_metadata_filter = len(filtered_accessions)

        # Save filtered metadata immediately (before sequence-dependent fields)
        logger.debug("Writing filtered metadata (pre-sequence) to JSONL: %s", output_metadata_jsonl)
        try:
            with open(output_metadata_jsonl, "w", encoding="utf-8") as f:
                for md in filtered_metadata:
                    f.write(json.dumps(md) + "\n")
            logger.info("✅ Saved filtered metadata JSONL: %s", output_metadata_jsonl)
        except Exception as e:
            logger.warning("❌ Failed to save filtered metadata JSONL: %s", e)

    # SECTION 5: DOWNLOAD SEQUENCES FOR FILTERED ACCESSIONS ONLY
        logger.info("=" * 60)
        logger.info("STEP 5: Downloading sequences for filtered accessions")
        logger.info("=" * 60)

        # Check if we're using cached sequences
        if used_cached_download and cached_sequences:
            logger.info("Using sequences from cached download (skipping sequence download)")
            
            # Filter cached sequences based on filtered_accessions
            # Create a mapping of accession to sequence
            cached_seq_dict = {seq.id.split('.')[0]: seq for seq in cached_sequences}
            
            # Get only the sequences that passed metadata filtering
            filtered_cached_seqs = []
            for acc in filtered_accessions:
                if acc in cached_seq_dict:
                    filtered_cached_seqs.append(cached_seq_dict[acc])
                else:
                    logger.debug("Accession %s not found in cached sequences", acc)
            
            logger.info("After metadata filtering: %d sequences from cache", len(filtered_cached_seqs))
            
            # Save to temporary FASTA file for consistency with regular pipeline
            fna_file = os.path.join(temp_dir, f"{virus_clean}_cached_sequences.fasta")
            FastaIO.write(filtered_cached_seqs, fna_file, "fasta")
            logger.info("Saved filtered cached sequences to temporary file: %s", fna_file)
        else:
            # Regular sequence download
            fna_file = download_sequences_by_accessions(filtered_accessions, outdir=temp_dir, failed_commands=failed_commands)
            if not os.path.exists(fna_file):
                raise RuntimeError(f"❌ Download failed: FASTA file not found at {fna_file}")
            logger.info("Downloaded FASTA file: %s (%.2f MB)", fna_file, os.path.getsize(fna_file) / 1024 / 1024)

    # SECTION 6: SEQUENCE-DEPENDENT FILTERING
        logger.info("=" * 60)
        logger.info("STEP 6: Applying sequence-dependent filters and saving results")
        logger.info("=" * 60)

        filters_seq={
            "max_ambiguous_chars": max_ambiguous_chars,
            "has_proteins": has_proteins,
            "proteins_complete": proteins_complete,
        }

        if filters_seq["max_ambiguous_chars"] is None and filters_seq["has_proteins"] is None and not filters_seq["proteins_complete"]:
            logger.info("No sequence-dependent filters specified, skipping this step.")
            filtered_sequences = list(FastaIO.parse(fna_file, "fasta"))
            filtered_metadata_final = filtered_metadata  # No change to metadata
            protein_headers = []
            logger.info("All %d downloaded sequences will be saved", len(filtered_sequences))
        else:
            # Restrict metadata to filtered accessions only
            filtered_metadata_dict = {acc: metadata_dict[acc] for acc in filtered_accessions}

            filtered_sequences, filtered_metadata_final, protein_headers = filter_sequences(
                fna_file,
                filtered_metadata_dict,
                **filters_seq,
            )

    # SECTION 7: SAVING FINAL OUTPUT FILES
        logger.info("=" * 60)
        logger.info("STEP 7: Saving final output files")
        logger.info("=" * 60)

        if filtered_sequences:
            logger.info("Saving %d filtered sequences and their metadata...", len(filtered_sequences))

            # Save FASTA
            FastaIO.write(filtered_sequences, output_fasta_file, "fasta")
            if os.path.exists(output_fasta_file):
                logger.info("✅ FASTA file saved: %s (%.2f MB)", output_fasta_file, os.path.getsize(output_fasta_file) / 1024 / 1024)
                output_files_dict['FASTA Sequences'] = output_fasta_file
            else:
                logger.error("❌ Failed to create FASTA file: %s", output_fasta_file)
            
            # Track final counts and metadata for summary
            total_final_sequences = len(filtered_sequences)
            final_metadata_for_summary = filtered_metadata_final

            # Overwrite JSONL with final filtered metadata (includes only sequences that passed all filters)
            try:
                with open(output_metadata_jsonl, "w", encoding="utf-8") as file:
                    for metadata in filtered_metadata_final:
                        file.write(json.dumps(metadata) + "\n")
                logger.info("✅ JSONL metadata file saved: %s (%.2f MB)", output_metadata_jsonl, os.path.getsize(output_metadata_jsonl) / 1024 / 1024)
                output_files_dict['JSONL Metadata'] = output_metadata_jsonl
            except Exception as e:
                logger.error("❌ Failed to save JSONL metadata file: %s", e)
                raise

            # CSV
            try:
                save_metadata_to_csv(filtered_metadata_final, protein_headers, output_metadata_csv)
                if os.path.exists(output_metadata_csv):
                    logger.info("✅ CSV metadata file saved: %s (%.2f MB)", output_metadata_csv, os.path.getsize(output_metadata_csv) / 1024 / 1024)
                    output_files_dict['CSV Metadata'] = output_metadata_csv
                else:
                    logger.error("❌ Failed to create CSV file: %s", output_metadata_csv)
            except Exception as e:
                logger.error("❌ Failed to save CSV metadata file: %s", e)
                raise

        # SECTION 8: GENBANK METADATA RETRIEVAL (OPTIONAL)
        logger.info("=" * 60)
        logger.info("STEP 8: Fetching detailed GenBank metadata")
        logger.info("=" * 60)
        logger.info("GenBank metadata retrieval requested - fetching detailed information...")
            
        # genbank_csv_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata.csv")
        if genbank_metadata:
            try:
                # Extract accession numbers from filtered sequences
                final_accessions = []
                for seq_record in filtered_sequences:
                    acc = seq_record.id.split()[0] if hasattr(seq_record, 'id') else str(seq_record)
                    if acc:
                        final_accessions.append(acc)
                
                if final_accessions:
                    logger.info("Fetching GenBank metadata for %d sequences...", len(final_accessions))
                    
                    # genbank_full_xml_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata_full_xml_data.xml")
                    # genbank_full_csv_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata_full_csv_data.csv")

                    # Fetch GenBank metadata
                    genbank_data, genbank_failed_log = fetch_genbank_metadata(
                        accessions=list(set(final_accessions)),  # Remove duplicates
                        genbank_full_xml_path=genbank_full_xml_path, genbank_full_csv_path=genbank_full_csv_path,
                        batch_size=genbank_batch_size,
                        delay=0.5  # Be respectful to NCBI servers
                    )
                    
                    # Parse GenBank failed batches log if it exists
                    if genbank_failed_log and os.path.exists(genbank_failed_log):
                        try:
                            with open(genbank_failed_log, 'r') as flog:
                                content = flog.read()
                                # Parse the log file to extract failed batch information
                                batch_pattern = r'FAILED_BATCH: \[([^\]]+)\][\s\S]*?URL: ([^\n]+)'
                                matches = re.findall(batch_pattern, content)
                                for accessions_str, url in matches:
                                    # Clean up accessions string
                                    batch_accessions = [acc.strip().strip("'").strip('"') for acc in accessions_str.split(',')]
                                    failed_commands['genbank_batches'].append({
                                        'accessions': batch_accessions,
                                        'retry_url': url.strip()
                                    })
                        except Exception as parse_error:
                            logger.debug("Could not parse GenBank failed batches log: %s", parse_error)
                    
                    if genbank_data:
                        # Save GenBank metadata to CSV
                        save_genbank_metadata_to_csv(
                            genbank_metadata=genbank_data,
                            output_file=genbank_csv_path,
                            virus_metadata=filtered_metadata_final
                        )
                        logger.info("✅ GenBank metadata CSV saved: %s (%.2f MB)", 
                                    genbank_csv_path, os.path.getsize(genbank_csv_path) / 1024 / 1024)
                        output_files_dict['GenBank CSV Metadata'] = genbank_csv_path
                        if os.path.exists(genbank_full_xml_path):
                            output_files_dict['GenBank Full XML'] = genbank_full_xml_path
                        if os.path.exists(genbank_full_csv_path):
                            output_files_dict['GenBank Full CSV'] = genbank_full_csv_path
                        genbank_success = True
                    else:
                        logger.warning("No GenBank metadata was retrieved")
                        genbank_error_msg = "No GenBank metadata was retrieved"
                else:
                    logger.warning("No accession numbers found for GenBank metadata lookup")
                    genbank_error_msg = "No accession numbers found for GenBank metadata lookup"
                    
            except Exception as genbank_error:
                logger.error("❌ GenBank metadata retrieval failed: %s", genbank_error)
                logger.warning("Continuing without GenBank metadata - standard output files are still available")
                genbank_error_msg = str(genbank_error)
                # Don't raise the error - continue with the rest of the process
            
            logger.info("GenBank metadata processing completed")
        else:
            logger.info("GenBank metadata retrieval not requested, skipping this step")

        # SECTION 9: FINAL SUMMARY
        # Provide comprehensive summary of the results
        if filtered_sequences:
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
            if not os.path.exists(genbank_csv_path):
                logger.info("  📊 Metadata (CSV):    %s", os.path.basename(output_metadata_csv))
            logger.info("  🔧 Metadata (JSONL):  %s", os.path.basename(output_metadata_jsonl))

            # Check if GenBank metadata CSV was created
            if genbank_metadata:
                if genbank_success and os.path.exists(genbank_csv_path):
                    logger.info("  📊 Metadata (including Genbank information):  %s", os.path.basename(genbank_csv_path))
                    if os.path.exists(genbank_full_xml_path):
                        logger.info("  🧬 GenBank-only full XML:      %s", os.path.basename(genbank_full_xml_path))
                    if os.path.exists(genbank_full_csv_path):
                        logger.info("  🧬 GenBank-only full CSV (readable XML format):      %s", os.path.basename(genbank_full_csv_path))
                else:
                    logger.warning("")
                    logger.warning("⚠️  GenBank metadata was requested but NOT saved due to errors:")
                    logger.warning("    %s", genbank_error_msg)
                    logger.warning("    Standard metadata files are still available.")

            logger.info("=" * 60)
            
            # Save command summary
            save_command_summary(
                outfolder=outfolder,
                command_line=command_line,
                total_api_records=total_api_records,
                total_after_metadata_filter=total_after_metadata_filter,
                total_final_sequences=total_final_sequences,
                output_files=output_files_dict,
                filtered_metadata=final_metadata_for_summary,
                success=True,
                error_message=None,
                failed_commands=failed_commands,
                genbank_error=genbank_error_msg if genbank_metadata and not genbank_success else None
            )
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
            
            # Save command summary even when no sequences pass
            save_command_summary(
                outfolder=outfolder,
                command_line=command_line,
                total_api_records=total_api_records,
                total_after_metadata_filter=total_after_metadata_filter,
                total_final_sequences=0,
                output_files=output_files_dict,
                filtered_metadata=[],
                success=True,
                error_message="No sequences passed all filters",
                failed_commands=failed_commands
            )

    except Exception as e:
        # Handle any unexpected errors during processing
        error_msg = str(e)
        
        # Check if this is a server-side issue that we can provide guidance for
        if any(indicator in error_msg.lower() for indicator in ['timeout', '500 server error', 'internal server error']):
            logger.error("=" * 80)
            logger.error("❌ SERVER-SIDE ERROR DETECTED")
            logger.error("=" * 80)
            logger.error("The NCBI API is experiencing server-side issues.")
            logger.error("This is not a problem with your query or parameters.")
            logger.error("")
            logger.error("Error details: %s", e)
            logger.error("")
            
            # Provide alternative commands based on the problematic parameters
            if geographic_location:
                logger.error("🔧 SUGGESTED SOLUTION:")
                logger.error("The geographic location filter appears to be causing server issues.")
                logger.error("Try running without the geographic filter and filter manually afterward:")
                logger.error("")
                
                # Build alternative command
                cmd_parts = [f"gget.virus('{virus}'"]
                
                # Add all non-problematic filters
                if host:
                    cmd_parts.append(f"host='{host}'")
                if min_seq_length:
                    cmd_parts.append(f"min_seq_length={min_seq_length}")
                if max_seq_length:
                    cmd_parts.append(f"max_seq_length={max_seq_length}")
                if min_gene_count:
                    cmd_parts.append(f"min_gene_count={min_gene_count}")
                if max_gene_count:
                    cmd_parts.append(f"max_gene_count={max_gene_count}")
                if nuc_completeness:
                    cmd_parts.append(f"nuc_completeness='{nuc_completeness}'")
                if annotated is not None:
                    cmd_parts.append(f"annotated={annotated}")
                if lab_passaged is not None:
                    cmd_parts.append(f"lab_passaged={lab_passaged}")
                if min_collection_date:
                    cmd_parts.append(f"min_collection_date='{min_collection_date}'")
                if max_collection_date:
                    cmd_parts.append(f"max_collection_date='{max_collection_date}'")
                if source_database:
                    cmd_parts.append(f"source_database='{source_database}'")
                if min_release_date:
                    cmd_parts.append(f"min_release_date='{min_release_date}'")
                if max_release_date:
                    cmd_parts.append(f"max_release_date='{max_release_date}'")
                if max_ambiguous_chars is not None:
                    cmd_parts.append(f"max_ambiguous_chars={max_ambiguous_chars}")
                if has_proteins:
                    if isinstance(has_proteins, list):
                        cmd_parts.append(f"has_proteins={has_proteins}")
                    else:
                        cmd_parts.append(f"has_proteins='{has_proteins}'")
                
                cmd_parts.append(f"outfolder='{virus_clean}_data'")
                
                alternative_cmd = ", ".join(cmd_parts) + ")"
                logger.error("📋 ALTERNATIVE COMMAND:")
                logger.error("  %s", alternative_cmd)
                logger.error("")
                logger.error("After download completes, filter the output CSV file by the")
                logger.error("'Geographic Location' column to get sequences from '%s'.", geographic_location)
            
            elif any(x in virus.lower() for x in ['sars-cov-2', 'covid', 'influenza']) and not host:
                logger.error("🔧 SUGGESTED SOLUTION:")
                logger.error("Large datasets like '%s' may cause server timeouts.", virus)
                logger.error("Try adding a host filter to reduce the dataset size:")
                logger.error("")
                
                # Build alternative command with host filter
                cmd_parts = [f"gget.virus('{virus}'", "host='human'"]
                
                # Add existing filters
                if min_seq_length:
                    cmd_parts.append(f"min_seq_length={min_seq_length}")
                if max_seq_length:
                    cmd_parts.append(f"max_seq_length={max_seq_length}")
                if nuc_completeness:
                    cmd_parts.append(f"nuc_completeness='{nuc_completeness}'")
                if annotated is not None:
                    cmd_parts.append(f"annotated={annotated}")
                
                cmd_parts.append(f"outfolder='{virus_clean}_data'")
                
                alternative_cmd = ", ".join(cmd_parts) + ")"
                logger.error("📋 ALTERNATIVE COMMAND:")
                logger.error("  %s", alternative_cmd)
            
            else:
                logger.error("🔧 SUGGESTED SOLUTIONS:")
                logger.error("1. Wait a few minutes and try again (server issues are often temporary)")
                logger.error("2. Try using more specific filters to reduce dataset size")
                logger.error("3. Use host='human' filter if studying human pathogens")
                logger.error("4. Add date range filters to limit the time period")
            
            logger.error("=" * 80)
        else:
            # For non-server errors, show the original error message
            logger.error("An error occurred during virus data processing: %s", e)
            logger.error("Error type: %s", type(e).__name__)
            if logger.getEffectiveLevel() <= logging.DEBUG:
                logger.debug("Full traceback:\n%s", traceback.format_exc())
        
        # Save command summary with error information
        save_command_summary(
            outfolder=outfolder if outfolder else os.getcwd(),
            command_line=command_line,
            total_api_records=total_api_records,
            total_after_metadata_filter=total_after_metadata_filter,
            total_final_sequences=total_final_sequences,
            output_files=output_files_dict,
            filtered_metadata=final_metadata_for_summary,
            success=False,
            error_message=str(e),
            failed_commands=failed_commands
        )
        
        raise
        
    # SECTION 10: CLEANUP
    finally:
        # Always clean up temporary files, regardless of success or failure
        logger.debug("Performing cleanup...")
        if os.path.exists(temp_dir) and keep_temp is False:
            try:
                shutil.rmtree(temp_dir)
                if os.path.exists(output_api_metadata_jsonl):
                    os.remove(output_api_metadata_jsonl)
                # Only remove metadata CSV if GenBank was successfully retrieved (as it's superseded by GenBank CSV)
                if genbank_metadata and genbank_success and os.path.exists(output_metadata_csv):
                    os.remove(output_metadata_csv)
                logger.debug("✅ Cleaned up temporary directory: %s", temp_dir)
            except Exception as e:
                logger.warning("❌ Failed to clean up temporary directory %s: %s", temp_dir, e)
        elif keep_temp and os.path.exists(output_api_metadata_jsonl):
            logger.debug("Preserving temporary directory as per user request: %s", temp_dir)
            shutil.move(output_api_metadata_jsonl, os.path.join(temp_dir, os.path.basename(output_api_metadata_jsonl)))
            if genbank_metadata and genbank_success and os.path.exists(genbank_csv_path):
                shutil.move(output_metadata_csv, os.path.join(temp_dir, os.path.basename(output_metadata_csv)))
                
                
        
        logger.info("NCBI virus data retrieval process completed.")


if __name__ == "__main__":
    # Main module entry point
    pass
