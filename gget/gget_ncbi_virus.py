import os
import re
import json
import time          # For adding delays between requests
import logging       # For logging level checks
import shutil        # For directory operations  
import subprocess    # For executing external commands
import traceback     # For error traceback logging
import pandas as pd  # For data manipulation and CSV output
import requests      # For HTTP requests to NCBI API
import zipfile       # For extracting downloaded ZIP files
from datetime import datetime
from dateutil import parser  # For flexible date parsing

# Internal imports for logging, unique ID generation, and FASTA parsing
from .utils import set_up_logger, FastaIO
from .constants import NCBI_API_BASE, NCBI_EUTILS_BASE

# Set up logger for this module
logger = set_up_logger()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
random_suffix = os.urandom(4).hex() # random suffix for naming uniqueness

def fetch_virus_metadata(
    virus,
    accession=False,
    host=None,
    geographic_location=None,
    # geographic_region=None,
    annotated=None,
    complete_only=None,
    min_release_date=None,
    refseq_only=None,
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
    loop = True
    # Retry logic for handling intermittent server issues
    max_retries = 3
    retry_delay = 2.0  # Start with 2 second delay
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
                response = requests.get(url, params=params, timeout=30)
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
                    logger.warning("Request failed (attempt %d/%d): %s. Retrying in %.1f seconds...", 
                                 attempt + 1, max_retries, e, retry_delay)
                    time.sleep(retry_delay)
                    retry_delay *= 1.5  # Exponential backoff
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
                raise RuntimeError(f"Connection error while fetching virus metadata: {last_exception}") from last_exception
                
            elif isinstance(last_exception, requests.exceptions.HTTPError):
                # Handle HTTP errors with specific guidance for known issues
                error_msg = f"HTTP error while fetching virus metadata: {last_exception}"
                
                # Check for specific server error patterns
                if "500" in str(last_exception):
                    error_msg += (
                        f"\n\n🔧 SERVER ERROR DETECTED: "
                        f"NCBI's API is experiencing temporary server-side issues. "
                        f"This is not a problem with your query. Try again in a few minutes, "
                        f"or consider using more specific filters to reduce the dataset size."
                    )
                raise RuntimeError(error_msg) from last_exception
                
            else:
                # Handle any other request-related errors
                raise RuntimeError(f"Failed to fetch virus metadata: {last_exception}") from last_exception

    
    # Log the final results summary
    logger.info("Successfully retrieved %d virus records from NCBI API across %d pages", 
                len(all_reports), page_count)
    
    return all_reports


def is_sars_cov2_query(virus, accession=False):
    """
    Check if the query is for SARS-CoV-2 to determine if optimized cached downloads should be used.
    
    NCBI provides optimized cached data packages for SARS-CoV-2 that are faster and more reliable
    than the general API endpoints. This function detects SARS-CoV-2 queries so we can use
    the optimized download method.
    
    Args:
        virus (str): Virus taxon name/ID or accession number
        accession (bool): Whether virus parameter is an accession number
        
    Returns:
        bool: True if this is a SARS-CoV-2 query that should use cached downloads
    """
    if accession:
        # When in accession mode, let the user explicitly set is_sars_cov2=True
        # rather than trying to detect it
        return False
    
    # Check for common SARS-CoV-2 identifiers in taxon names
    virus_lower = virus.lower().replace('-', '').replace('_', '').replace(' ', '')
    sars_cov2_identifiers = {
        'sarscov2', 'sars2', '2697049', 'sarscov', 'severeacuterespiratorysyndromecoronavirus2',
        'covid19', 'covid', 'coronavirusdisease', 'ncov', 'hcov19'
    }
    
    # Check if the query matches any SARS-CoV-2 identifier
    for identifier in sars_cov2_identifiers:
        if identifier in virus_lower:
            logger.debug("Detected SARS-CoV-2 query: %s matches %s", virus, identifier)
            return True
    
    logger.debug("Not a SARS-CoV-2 query: %s", virus)
    return False


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
    1. If use_accession=True: Direct accession download using accession endpoint
    2. If use_accession=False:
       a. Specific lineage + complete + host filters using taxon endpoint
       b. Complete genomes only using taxon endpoint
       c. All SARS-CoV-2 genomes using taxon endpoint (default fallback)
    
    Args:
        host (str): Host organism filter (optimized for 'human', others handled in post-processing)
        complete_only (bool): Whether to download only complete genomes
        annotated (bool): Whether to download only annotated genomes  
        outdir (str): Output directory for downloaded files
        lineage (str): SARS-CoV-2 lineage filter (e.g., 'B.1.1.7', 'P.1')
        accession (str): Specific SARS-CoV-2 accession or taxon ID
        use_accession (bool): Whether to use the accession endpoint (True) or taxon endpoint (False)
        
    Returns:
        str: Path to the downloaded ZIP file containing sequences and metadata
        
    Raises:
        RuntimeError: If the datasets CLI is not available or download fails
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

    # Try each strategy in order
    last_error = None
    
    for strategy_name, cmd, applied_filters in strategies:
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

            # Start with capturing just stderr to check progress
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=outdir
            )
            
            # Track time and last progress
            start_time = time.time()
            last_progress = start_time
            last_stderr = ""
            
            while True:
                # Check if process has finished
                retcode = process.poll()
                if retcode is not None:
                    break
                    
                # Read stderr without blocking
                stderr = process.stderr.readline()
                if stderr:
                    last_stderr = stderr
                    # Common progress indicators in datasets tool output:
                    progress_indicators = ['%', '=', 'downloading', 'fetching', 'MB', 'GB', 'bytes']
                    
                    # Log the stderr for debugging
                    logger.debug("Progress output: %s", stderr.strip())
                    
                    # If we see any progress indicator, update the last_progress time
                    if any(indicator.lower() in stderr.lower() for indicator in progress_indicators):
                        last_progress = time.time()
                        logger.debug("Progress detected, updating last_progress time")
                
                # Check timeout conditions:
                # 1. Less than 30 minutes total time, continue
                # 2. If more than 30 mins but progress in last 5 mins, continue
                # 3. Otherwise, timeout
                current_time = time.time()
                total_time = current_time - start_time
                time_since_progress = current_time - last_progress
                
                if total_time > 1800 and time_since_progress > 300:  # 30 mins total and 5 mins no progress
                    process.kill()
                    raise subprocess.TimeoutExpired(cmd, 1800)
                
                time.sleep(0.1)  # Prevent CPU spin
            
            stdout, stderr = process.communicate()
            result = subprocess.CompletedProcess(
                args=cmd,
                returncode=retcode,
                stdout=stdout,
                stderr=stderr
            )
            
            # Check if the command was successful
            if result.returncode == 0 and os.path.exists(zip_path):
                file_size = os.path.getsize(zip_path)
                logger.info("✅ %s successful: %s (%.2f MB)", 
                           strategy_name, os.path.basename(zip_path), file_size / 1024 / 1024)
                
                # Log any important output from the datasets CLI
                if result.stdout:
                    logger.debug("datasets CLI output: %s", result.stdout.strip())
                
                # Check which filters from the original request weren't applied in this strategy
                requested_filters = []
                if complete_only: requested_filters.append("complete-only")
                if lineage: requested_filters.append(f"lineage={lineage}")
                if host: requested_filters.append(f"host={host}")
                if annotated: requested_filters.append("annotated")
                
                missing_filters = [f for f in requested_filters if f not in applied_filters]
                if missing_filters:
                    logger.warning("⚠️ Some requested filters were not applied in successful strategy:")
                    logger.warning("   Filters applied: %s", ", ".join(applied_filters) if applied_filters else "none")
                    logger.warning("   Filters missing: %s", ", ".join(missing_filters))
                    logger.warning("   These filters will need to be applied through post-processing")
                
                return zip_path
            else:
                # Strategy failed, prepare error message
                error_msg = f"{strategy_name} failed with return code {result.returncode}"
                if result.stderr:
                    error_msg += f": {result.stderr.strip()}"
                logger.warning("❌ %s", error_msg)
                last_error = error_msg
                
                # If this was an accession download that failed, provide specific guidance
                if use_accession:
                    error_msg = (
                        f"Failed to download SARS-CoV-2 sequence for accession '{accession}'. "
                        f"Please verify that this is a valid SARS-CoV-2 accession number. "
                        f"If you're not sure if this is a SARS-CoV-2 sequence, try without setting is_sars_cov2=True."
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
            
        except FileNotFoundError as e:
            # datasets CLI not found - this is a critical error, don't continue
            raise RuntimeError(
                "datasets CLI not found. Please install NCBI datasets CLI tools. "
                "Installation guide: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/"
            ) from e
            
        except Exception as e:
            error_msg = f"{strategy_name} unexpected error: {e}"
            logger.warning("❌ %s", error_msg)
            last_error = error_msg
            continue
    
    logger.warning("🚨 All cached download strategies failed. Last error: %s", last_error)
    
    # Provide helpful guidance
    guidance_messages = [
        "🔧 TROUBLESHOOTING SUGGESTIONS:",
        "1. Check your internet connection",
        "2. Verify datasets CLI is properly installed and updated",
        "3. Try running the command manually to see detailed error messages:",
        f"   datasets download virus genome taxon SARS-CoV-2 --filename test.zip",
        "4. NCBI servers may be temporarily unavailable - try again later",
        "5. Consider using the general API method by removing SARS-CoV-2 specific terms from your query"
    ]
    
    for msg in guidance_messages:
        logger.info(msg)
    
    # Raise error with the last failure details
    raise RuntimeError(
        f"All SARS-CoV-2 cached download strategies failed. "
        f"Last error: {last_error}. "
        f"Consider using the general API method instead by modifying your virus query to avoid SARS-CoV-2 detection. Use Taxon ID without is_sars_cov2 flag instead."
    )


def download_sequences_by_accessions(accessions, outdir=None, batch_size=200):
    """
    Download virus genome sequences for a specific list of accession numbers.
    
    This function downloads sequences for a pre-filtered list of accessions,
    using NCBI E-utilities API with batching to avoid URL length limitations.
    Large requests are automatically split into smaller batches.
    
    Args:
        accessions (list): List of accession numbers to download
        outdir (str): Output directory for downloaded files
        batch_size (int): Maximum number of accessions per batch (default: 200)
        
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
    fasta_path = os.path.join(outdir, f"virus_sequences_{timestamp}_{random_suffix}.fasta")
    logger.debug("Saving sequences to: %s", fasta_path)
    
    # Use NCBI E-utilities to fetch FASTA sequences with batching
    # base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    # If we have many accessions, split into batches to avoid URL length limits
    if len(accessions) > batch_size:
        logger.info("Large request detected (%d accessions). Using batch processing with batch size %d", 
                   len(accessions), batch_size)
        return _download_sequences_batched(accessions, NCBI_EUTILS_BASE, fasta_path, batch_size)
    
    # For smaller requests, use single request
    return _download_sequences_single_batch(accessions, NCBI_EUTILS_BASE, fasta_path)


def _download_sequences_single_batch(accessions, NCBI_EUTILS_BASE, fasta_path):
    """Download sequences in a single E-utilities request."""
    
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
        response = requests.get(NCBI_EUTILS_BASE, params=params, timeout=300)
        
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
        
        # Check for specific URL length error
        if "414" in str(e) or "Request-URI Too Long" in str(e):
            logger.info("URL too long error detected. Retrying with batch processing...")
            # Retry with smaller batches
            return _download_sequences_batched(accessions, NCBI_EUTILS_BASE, fasta_path, batch_size=100)
        raise RuntimeError(f"Failed to download virus sequences via E-utilities: {e}") from e
    except IOError as e:
        logger.error("Failed to save FASTA file: %s", e)
        raise RuntimeError(f"Failed to save downloaded sequences: {e}") from e


def _download_sequences_batched(accessions, NCBI_EUTILS_BASE, fasta_path, batch_size):
    """Download sequences using multiple batched E-utilities requests."""
    
    # Split accessions into batches
    batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
    logger.info("Downloading %d accessions in %d batches of size %d", 
               len(accessions), len(batches), batch_size)
    
    all_sequences = []
    total_downloaded = 0
    
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
            response = requests.get(NCBI_EUTILS_BASE, params=params, timeout=300)
            response.raise_for_status()
            
            # Verify we got FASTA data
            if not response.text.strip().startswith('>'):
                logger.warning("Batch %d returned invalid FASTA data: %s", 
                             batch_num, response.text[:100])
                continue
            
            # Count sequences in this batch
            batch_sequence_count = response.text.count('>')
            total_downloaded += batch_sequence_count
            
            # Store sequences from this batch
            all_sequences.append(response.text)
            
            logger.info("Batch %d: Downloaded %d sequences (%.2f MB)", 
                       batch_num, batch_sequence_count, 
                       len(response.text.encode('utf-8')) / 1024 / 1024)
            
            # Add small delay between requests to be respectful to NCBI servers
            if batch_num < len(batches):  # Don't delay after the last batch
                time.sleep(0.5)
                
        except requests.exceptions.RequestException as e:
            logger.error("Batch %d failed: %s", batch_num, e)
            
            # Check for URL length error even in batch mode
            if "414" in str(e) and batch_size > 50:
                logger.warning("URL still too long for batch size %d. Retrying batch %d with smaller size...", 
                             batch_size, batch_num)
                # Recursively retry this batch with smaller size by splitting it further
                temp_batch_path = f"temp_batch_{batch_num}.fasta"
                try:
                    _download_sequences_batched(
                        batch_accessions, NCBI_EUTILS_BASE, temp_batch_path, batch_size // 2
                    )
                    # Read the temporary file and add to all_sequences
                    with open(temp_batch_path, 'r', encoding='utf-8') as f:
                        batch_content = f.read()
                        all_sequences.append(batch_content)
                        # Count sequences in this recovered batch
                        recovered_count = batch_content.count('>')
                        total_downloaded += recovered_count
                        logger.info("Recovered batch %d with smaller size: %d sequences", 
                                   batch_num, recovered_count)
                    os.remove(temp_batch_path)  # Clean up temp file
                except Exception as file_error:
                    logger.warning("Failed to recover batch %d with smaller size: %s", batch_num, file_error)
                    continue
            else:
                logger.warning("Batch %d failed and will be skipped: %s", batch_num, e)
                continue
    
    if not all_sequences:
        raise RuntimeError("All batches failed. No sequences were downloaded.")
    
    # Combine all sequences and write to output file
    try:
        with open(fasta_path, 'w', encoding='utf-8') as f:
            for sequence_data in all_sequences:
                f.write(sequence_data)
                if not sequence_data.endswith('\n'):
                    f.write('\n')  # Ensure proper line endings between batches
        
        file_size = os.path.getsize(fasta_path)
        logger.info("Successfully saved %d sequences to: %s (%.2f MB)", 
                   total_downloaded, fasta_path, file_size / 1024 / 1024)
        return fasta_path
        
    except IOError as e:
        logger.error("Failed to save combined FASTA file: %s", e)
        raise RuntimeError(f"Failed to save downloaded sequences: {e}") from e


def unzip_file(zip_file_path, extract_to_path):
    """Unzips a ZIP file to a specified directory."""
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
                "location": report.get("location", {}).get("geographic_location", pd.NA),  # Geographic location details
                "region": report.get("location", {}).get("geographic_region", pd.NA),  # Broad region
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
            # logger.debug("Loaded metadata for accession: %s", accession)
            # logger.debug("Metadata details: %s", metadata)
            
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


def filter_sequences(
    fna_file,
    metadata_dict,
    max_ambiguous_chars=None,
    has_proteins=None,
    proteins_complete=False,
    # **kwargs,
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
        
        # Get sequence length
        # seq_length = len(record.seq)
        
        # # Apply sequence length filters
        # if min_seq_length is not None and seq_length < min_seq_length:
        #     filter_stats['seq_length'] += 1
        #     record_passes = False
        #     continue
            
        # if max_seq_length is not None and seq_length > max_seq_length:
        #     filter_stats['seq_length'] += 1
        #     record_passes = False
        #     continue
            
        # Count ambiguous characters (N's)
        if max_ambiguous_chars is not None:
            ambiguous_count = record.seq.upper().count('N')
            if ambiguous_count > max_ambiguous_chars:
                filter_stats['ambiguous_chars'] += 1
                record_passes = False
                continue
                
        # TODO: Add protein feature analysis if required (has_proteins/proteins_complete)
        
        # If sequence passed all filters, keep it and its metadata
        if record_passes:
            filtered_sequences.append(record)
            filtered_metadata.append(metadata_dict.get(record.id, {}))

    # Log filtering results
    logger.info("Sequence filter results:")
    logger.info("- Total sequences processed: %d", total_sequences)
    logger.info("- Failed length filter: %d", filter_stats['seq_length'])
    logger.info("- Failed ambiguous char filter: %d", filter_stats['ambiguous_chars'])
    logger.info("- Failed protein requirements: %d", filter_stats['proteins'])
    logger.info("- Sequences passing all filters: %d", len(filtered_sequences))
    
    return filtered_sequences, filtered_metadata, protein_headers
    
    # Log comprehensive filtering statistics
    # num_seqs = len(filtered_sequences)
    # logger.info("Filtering complete:")
    # logger.info("  Total sequences processed: %d", total_sequences)
    # logger.info("  Sequences with metadata: %d", sequences_with_metadata)
    # logger.info("  Sequences passing all filters: %d", num_seqs)
    
    # # Log detailed filter statistics if any sequences were filtered out
    # total_filtered = sum(filter_stats.values())
    # if total_filtered > 0:
    #     logger.info("Filter statistics (sequences excluded):")
    #     for filter_name, count in filter_stats.items():
    #         if count > 0:
    #             logger.info("  %s: %d sequences", filter_name, count)
    
    # # Validate results and return
    # if num_seqs > 0:
    #     # Perform consistency checks
    #     if num_seqs != len(filtered_metadata):
    #         logger.warning(
    #             "Number of sequences (%d) and number of metadata entries (%d) do not match.", 
    #             num_seqs, len(filtered_metadata)
    #         )
    #     if num_seqs != len(protein_headers):
    #         logger.warning(
    #             "Number of sequences (%d) and number of protein headers (%d) do not match.", 
    #             num_seqs, len(protein_headers)
    #         )
            
    #     logger.debug(
    #         "Final counts - sequences: %d, metadata: %d, protein headers: %d", 
    #         num_seqs, len(filtered_metadata), len(protein_headers)
    #     )
        
    #     return filtered_sequences, filtered_metadata, protein_headers
    # else:
    #     logger.warning("No sequences passed the provided filters.")
    #     return None, None, None


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
        # logger.debug("Raw metadata for record %d: %s", i+1, metadata)
        # location_info = metadata.get("location", {})
        # logger.debug("Raw location info: %s", location_info)
        
        # Filter out empty values and reverse the order for hierarchical display
        # location_values = [v for v in location_info.values() if v and v != ""]
        # location_values.reverse()
        
        # Create a colon-separated geographic string (e.g., "Africa:South Africa:Cape Town")
        # geo_info = ":".join(location_values) if location_values else pd.NA
        
        # Extract geographic region (broadest level, e.g., "Africa", "Europe")
        # try:
        #     geo_region = location_values[0] if location_values else pd.NA
        # except IndexError:
        #     geo_region = pd.NA
            
        # Extract specific geographic location (e.g., "South Africa", "Germany")
        # try:
        #     geo_loc = location_values[1] if len(location_values) > 1 else pd.NA
        # except IndexError:
        #     geo_loc = pd.NA

        # logger.debug("Geographic info for record %d: region=%s, location=%s, full=%s", 
        #             i+1, geo_region, geo_loc, geo_info)

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
            # "Geographic Region": geo_region,
            # "Geographic Location": geo_loc,
            "Geographic Region": metadata.get("region", pd.NA),
            "Geographic Location": metadata.get("location", pd.NA),
            # "Geo String": geo_info,
            
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


def fetch_genbank_metadata(accessions, batch_size=200, delay=0.5):
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
    
    if not accessions:
        raise ValueError("No accessions provided for GenBank metadata retrieval")
    
    logger.info("Fetching GenBank metadata for %d accessions using E-utilities", len(accessions))
    logger.debug("First 5 accessions: %s", accessions[:5])
    
    # Dictionary to store all metadata results
    all_metadata = {}
    
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
            batch_metadata = _fetch_genbank_batch(batch_accessions)
            all_metadata.update(batch_metadata)
            
            logger.info("Batch %d: Successfully retrieved metadata for %d accessions", 
                       batch_num, len(batch_metadata))
            
            # Add delay between requests to be respectful to NCBI servers
            if batch_num < len(batches) and delay > 0:
                logger.debug("Adding %.1f second delay before next batch", delay)
                time.sleep(delay)
                
        except Exception as e:
            logger.error("Batch %d failed: %s", batch_num, e)
            # Continue with remaining batches rather than failing completely
            continue
    
    logger.info("GenBank metadata retrieval complete: %d/%d accessions processed", 
                len(all_metadata), len(accessions))
    
    if not all_metadata:
        logger.warning("No GenBank metadata was successfully retrieved")
    
    return all_metadata


def _fetch_genbank_batch(accessions):
    """
    Fetch GenBank metadata for a single batch of accessions.
    
    Args:
        accessions (list): List of accession numbers for this batch
        
    Returns:
        dict: Dictionary mapping accession numbers to metadata dictionaries
        
    Raises:
        RuntimeError: If the E-utilities request fails or XML parsing fails
    """
    
    # Build E-utilities efetch URL for GenBank XML format
    # base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    accession_string = ",".join(accessions)
    
    params = {
        'db': 'nucleotide',           # Nucleotide database (includes virus sequences)
        'id': accession_string,       # Comma-separated accession numbers
        'rettype': 'gb',              # GenBank format
        'retmode': 'xml',             # XML output for structured parsing
        'tool': 'gget',               # Identify our application
        'email': 'gget@example.com'   # Required by NCBI guidelines
    }
    
    try:
        logger.debug("Making E-utilities request for %d accessions", len(accessions))
        logger.debug("Request URL: %s", NCBI_EUTILS_BASE)
        logger.debug("Request parameters: %s", {k: v[:50] + '...' if len(str(v)) > 50 else v for k, v in params.items()})
        
        # Make the HTTP request with extended timeout for large batches
        response = requests.get(NCBI_EUTILS_BASE, params=params, timeout=300)
        response.raise_for_status()
        
        # Verify we got XML data
        if not response.text.strip().startswith('<?xml') and not response.text.strip().startswith('<'):
            raise RuntimeError(f"Invalid XML response: {response.text[:200]}")
        
        logger.debug("Received XML response: %d characters", len(response.text))
        
        # Parse the GenBank XML and extract metadata
        metadata_dict = _parse_genbank_xml(response.text)
        
        logger.debug("Successfully parsed metadata for %d records", len(metadata_dict))
        return metadata_dict
        
    except requests.exceptions.RequestException as e:
        logger.error("E-utilities request failed: %s", e)
        raise RuntimeError(f"Failed to fetch GenBank metadata via E-utilities: {e}") from e
    except Exception as e:
        logger.error("GenBank XML parsing failed: %s", e)
        raise RuntimeError(f"Failed to parse GenBank XML response: {e}") from e


def _parse_genbank_xml(xml_content):
    """
    Parse GenBank XML response and extract relevant metadata fields.
    
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
    
    import xml.etree.ElementTree as ET
    
    try:
        # Parse the XML content
        root = ET.fromstring(xml_content)
        logger.debug("XML parsing successful, root element: %s", root.tag)
        
    except ET.ParseError as e:
        logger.error("XML parsing failed: %s", e)
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
                ref_data['authors'] = authors_elem.text if authors_elem is not None else ""
                
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
            metadata['genbank_data']['country'] = source_feature.get('country', '')
            metadata['genbank_data']['geographic_location'] = source_feature.get('geo_loc_name', '')
            metadata['genbank_data']['host'] = source_feature.get('host', '')
            metadata['genbank_data']['isolation_source'] = source_feature.get('isolation_source', '')
            metadata['genbank_data']['strain'] = source_feature.get('strain', '')
            metadata['genbank_data']['isolate'] = source_feature.get('isolate', '')
            metadata['genbank_data']['collected_by'] = source_feature.get('collected_by', '')
            metadata['genbank_data']['specimen_voucher'] = source_feature.get('specimen_voucher', '')
            
            # Store all features for potential future use
            metadata['genbank_data']['all_features'] = features_data
            
            # Extract comment field (often contains additional metadata)
            comment_elem = gbseq.find('GBSeq_comment')
            comment_text = comment_elem.text if comment_elem is not None else ""
            metadata['genbank_data']['comment'] = comment_text
            
            # Parse assembly name from comment if present (used in some studies)
            assembly_name = ""
            if comment_text:
                import re
                assembly_match = re.search(r'Assembly Name :: (\S+)', comment_text)
                if assembly_match:
                    assembly_name = assembly_match.group(1)
            metadata['genbank_data']['assembly_name'] = assembly_name
            
            # Store the metadata for this accession
            metadata_dict[accession] = metadata
            
            logger.debug("Extracted GenBank metadata for %s: organism=%s, collection_date=%s, country=%s", 
                        accession, 
                        metadata['genbank_data']['organism'],
                        metadata['genbank_data']['collection_date'],
                        metadata['genbank_data']['country'])
            
        except Exception as e:
            logger.warning("Failed to parse GenBank record %s: %s", 
                          accession if 'accession' in locals() else 'unknown', e)
            continue
    
    logger.info("Successfully parsed GenBank metadata for %d records", len(metadata_dict))
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
        "country",
        "geographic_location",
        "strain",
        "isolate",
        
        # Host and source information
        "host",
        "isolation_source",
        "collected_by",
        "specimen_voucher",
        
        # Database and version information
        "create_date",
        "update_date",
        "assembly_name",
        
        # Taxonomic information
        "taxonomy",
        
        # Publication information
        "first_author",
        "first_title",
        "first_journal",
        "first_pubmed_id",
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
            "country": genbank_data.get('country', ''),
            "geographic_location": genbank_data.get('geographic_location', ''),
            "strain": genbank_data.get('strain', ''),
            "isolate": genbank_data.get('isolate', ''),
            
            "host": genbank_data.get('host', ''),
            "isolation_source": genbank_data.get('isolation_source', ''),
            "collected_by": genbank_data.get('collected_by', ''),
            "specimen_voucher": genbank_data.get('specimen_voucher', ''),
            
            "create_date": genbank_data.get('create_date', ''),
            "update_date": genbank_data.get('update_date', ''),
            "assembly_name": genbank_data.get('assembly_name', ''),
            
            "taxonomy": genbank_data.get('taxonomy', ''),
            
            "first_author": first_ref.get('authors', ''),
            "first_title": first_ref.get('title', ''),
            "first_journal": first_ref.get('journal', ''),
            "first_pubmed_id": first_ref.get('pubmed_id', ''),
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
                df = pd.merge(df, virus_df, on='accession', how='left', suffixes=('_genbank', '_virus'))
                logger.info("Successfully merged GenBank and virus metadata")
            else:
                logger.warning("Cannot merge: virus metadata missing 'accession' column")
        except Exception as e:
            logger.warning("Failed to merge GenBank and virus metadata: %s", e)
    
    # Write DataFrame to CSV file
    try:
        df.to_csv(output_file, index=False, encoding='utf-8')
        logger.info("GenBank metadata saved to: %s", output_file)
        logger.info("CSV file contains %d rows and %d columns", len(df), len(df.columns))
        
        # Log some summary statistics
        non_empty_collection_dates = df['collection_date'].astype(str).str.strip().ne('').sum()
        non_empty_countries = df['country'].astype(str).str.strip().ne('').sum()
        non_empty_hosts = df['host'].astype(str).str.strip().ne('').sum()
        
        logger.info("GenBank metadata summary:")
        logger.info("  Records with collection dates: %d/%d (%.1f%%)", 
                   non_empty_collection_dates, len(df), 100 * non_empty_collection_dates / len(df))
        logger.info("  Records with country information: %d/%d (%.1f%%)", 
                   non_empty_countries, len(df), 100 * non_empty_countries / len(df))
        logger.info("  Records with host information: %d/%d (%.1f%%)", 
                   non_empty_hosts, len(df), 100 * non_empty_hosts / len(df))
        
    except Exception as e:
        logger.error("Failed to save GenBank metadata CSV: %s", e)
        raise RuntimeError(f"Failed to save GenBank metadata to {output_file}: {e}") from e
    

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
    # Sequence-dependent filters are deferred
    # max_ambiguous_chars=None,
    # has_proteins=None,
    # proteins_complete=False,
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
        # 'geographic_region': 0,
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
    
    # if max_ambiguous_chars is not None or has_proteins is not None:
    #     logger.info("Note: Sequence-dependent filters (max_ambiguous_chars, has_proteins) will be applied after download")
    
    return filtered_accessions, filtered_metadata_list


def ncbi_virus(
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
    # host_taxid=None,
    lab_passaged=None,
    # geographic_region=None,
    geographic_location=None,
    submitter_country=None,
    min_collection_date=None,
    max_collection_date=None,
    annotated=None,
    refseq_only=False,
    source_database=None,
    min_release_date=None,
    max_release_date=None,
    min_mature_peptide_count=None,
    max_mature_peptide_count=None,
    min_protein_count=None,
    max_protein_count=None,
    max_ambiguous_chars=None,
    is_sars_cov2=False,
    lineage=None,
    genbank_metadata=False,
    genbank_batch_size=200,
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
        lineage (str): Virus lineage filter (SARS-CoV-2 specific)
        genbank_metadata (bool): Whether to fetch detailed GenBank metadata (default: False)
        genbank_batch_size (int): Batch size for GenBank API requests (default: 200)
        is_sars_cov2 (bool): Flag to indicate if the accession is for SARS-CoV-2, enabling optimized download method (default: False)

    Returns:
        None: Files are saved to the output directory
        
    Note:
        When genbank_metadata=True, an additional CSV file with detailed GenBank
        metadata will be saved alongside the standard output files. This includes
        collection dates, geographic information, host details, publication
        references, and other fields extracted from GenBank records.
    """
    logger.info("Starting NCBI virus data retrieval process...")
    logger.info("Query parameters: virus='%s', is_accession=%s, outfolder='%s'", 
                virus, is_accession, outfolder)
    logger.debug("Applied filters: host=%s, seq_length=(%s-%s), gene_count=(%s-%s), completeness=%s, annotated=%s, refseq_only=%s, lab_passaged=%s, geo_location=%s, submitter_country=%s, collection_date=(%s-%s), source_db=%s, release_date=(%s-%s), protein_count=(%s-%s), peptide_count=(%s-%s), max_ambiguous=%s, has_proteins=%s, proteins_complete=%s, genbank_metadata=%s, genbank_batch_size=%s",
    host, min_seq_length, max_seq_length, min_gene_count, max_gene_count, nuc_completeness, annotated, refseq_only, lab_passaged, geographic_location, submitter_country, min_collection_date, max_collection_date,source_database, min_release_date, max_release_date, min_protein_count, max_protein_count, min_mature_peptide_count, max_mature_peptide_count, max_ambiguous_chars, has_proteins, proteins_complete, genbank_metadata, genbank_batch_size)

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

    if refseq_only is not None and not isinstance(refseq_only, bool):
        raise TypeError(
            "Argument 'refseq_only' must be a boolean (True or False)."
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
        if genbank_batch_size > 500:
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
    
    # SECTION 2.5: SARS-CoV-2 CACHED DATA PROCESSING
    # For SARS-CoV-2 queries, use cached data packages with hierarchical fallback
    if is_sars_cov2 or is_sars_cov2_query(virus, is_accession):
        logger.info("=" * 60)
        logger.info("DETECTED SARS-CoV-2 QUERY - USING CACHED DATA PACKAGES")
        logger.info("=" * 60)
        logger.info("SARS-CoV-2 queries will use NCBI's optimized cached data packages")
        logger.info("with hierarchical fallback from specific to general cached files.")
        
        # Use the download_sars_cov2_optimized function which handles fallback strategies internally
        params = {
            'host': host,
            'complete_only': (nuc_completeness == "complete"),
            'annotated': annotated,
            'outdir': outfolder,
            'lineage': lineage,
            'accession': virus,  # Pass the virus parameter as either accession or taxon
            'use_accession': is_accession  # Tell the function whether to use accession or taxon endpoint
        }
            
        zip_file = download_sars_cov2_optimized(**params)
            
        if zip_file and os.path.exists(zip_file):
            # Extract directory path should come from the zip file name to match what download_sars_cov2_optimized uses
            extract_dir = os.path.splitext(zip_file)[0]
            unzip_file(zip_file, extract_dir)
            
            # Process the cached download
            if os.path.exists(extract_dir):
                logger.info("🔬 PROCESSING CACHED DATA...")
                logger.info("Extracted cached data to: %s", extract_dir)
            
            # Process cached data pipeline
            try:
                # 1. Find FASTA and metadata files in extracted directory
                fasta_files = []
                metadata_files = []
                
                for root, dirs, files in os.walk(extract_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        if file.endswith(('.fasta', '.fa', '.fna')):
                            fasta_files.append(file_path)
                        elif file.endswith(('.json', '.jsonl', '.csv')):
                            metadata_files.append(file_path)
                
                # Process all available FASTA files
                all_cached_sequences = []
                for fasta_file in fasta_files:
                    try:
                        sequences = list(FastaIO.parse(fasta_file, "fasta"))
                        all_cached_sequences.extend(sequences)
                        logger.info("Loaded %d sequences from %s", len(sequences), fasta_file)
                    except Exception as e:
                        logger.warning("Failed to load FASTA file %s: %s", fasta_file, e)
                        continue
                
                if not all_cached_sequences:
                    logger.warning("No valid sequences found in cached data. Falling back to API workflow.")
                    raise RuntimeError("No valid sequences found in cached data")
                
                # 3. For cached data, create minimal metadata (accession-based)
                cached_metadata = []
                for seq in all_cached_sequences:
                    # Extract basic info from FASTA header
                    metadata = {
                        'accession': seq.id.split('.')[0],  # Remove version if present
                        'description': seq.description,
                        'length': len(seq.seq),
                        'source': 'cached_data'
                    }
                    cached_metadata.append(metadata)
                
                logger.info("Created basic metadata for %d cached sequences", len(cached_metadata))
                
                # 4. Apply basic filtering to cached data
                # Note: Some filters may not work fully with cached data due to limited metadata
                logger.info("Applying basic filters to cached data...")
                
                # Apply sequence length filters if specified
                if min_seq_length is not None or max_seq_length is not None:
                    filtered_seqs = []
                    filtered_meta = []
                    for seq, meta in zip(all_cached_sequences, cached_metadata):
                        seq_len = len(seq.seq)
                        if min_seq_length is not None and seq_len < min_seq_length:
                            continue
                        if max_seq_length is not None and seq_len > max_seq_length:
                            continue
                        filtered_seqs.append(seq)
                        filtered_meta.append(meta)
                    all_cached_sequences = filtered_seqs
                    cached_metadata = filtered_meta
                    logger.info("After length filtering: %d sequences remain", len(all_cached_sequences))
                
                # Apply accession filter if specified
                if accession and isinstance(accession, (list, tuple)):
                    accession_set = set(str(acc).upper() for acc in accession)
                    filtered_seqs = []
                    filtered_meta = []
                    for seq, meta in zip(all_cached_sequences, cached_metadata):
                        if meta['accession'].upper() in accession_set:
                            filtered_seqs.append(seq)
                            filtered_meta.append(meta)
                    all_cached_sequences = filtered_seqs
                    cached_metadata = filtered_meta
                    logger.info("After accession filtering: %d sequences remain", len(all_cached_sequences))
                
                # 5. Save filtered cached results
                if all_cached_sequences:
                    logger.info("🎉 CACHED DATA PROCESSING SUCCESSFUL!")
                    logger.info("Using %d sequences from cached data", len(all_cached_sequences))
                    
                    # Save FASTA sequences
                    output_fasta_path = os.path.join(outfolder, f"{virus}_sequences.fasta")
                    FastaIO.write(all_cached_sequences, output_fasta_path, "fasta")
                    logger.info("Saved cached sequences to: %s", output_fasta_path)
                    
                    # Save metadata CSV
                    if cached_metadata:
                        output_csv_path = os.path.join(outfolder, f"{virus}_metadata.csv")
                        import pandas as pd
                        df = pd.DataFrame(cached_metadata)
                        df.to_csv(output_csv_path, index=False)
                        logger.info("Saved cached metadata to: %s", output_csv_path)
                    
                    logger.info("✅ Cached data processing completed successfully!")
                    return output_fasta_path  # Return early with cached results
                else:
                    logger.warning("No sequences remain after filtering cached data")
                    raise RuntimeError("No sequences passed the filter criteria")
                    
            except Exception as cache_error:
                logger.warning("Cached data processing failed: %s", cache_error)
                logger.info("Proceeding with regular API workflow...")
                logger.info("=" * 60)
    
    # Create temporary directory for intermediate processing
    # This will be cleaned up at the end regardless of success or failure
    temp_dir = os.path.join(outfolder, f"tmp_{timestamp}_{random_suffix}")
    os.makedirs(temp_dir, exist_ok=True)
    logger.debug("Created temporary processing directory: %s", temp_dir)
    
    try:
        # SECTION 3: METADATA RETRIEVAL WHILE APPLYING SERVER-SIDE FILTERS
        logger.info("=" * 60)
        logger.info("STEP 2: Fetching virus metadata from NCBI API")
        logger.info("=" * 60)

        api_annotated_filter = annotated if annotated is True else None
        api_complete_filter = True if nuc_completeness=="complete" else False

        logger.debug("Applying server-side filters: host=%s, geo_location=%s, annotated=%s, complete_only=%s, min_release_date=%s, refseq_only=%s", host, geographic_location, annotated, api_complete_filter, min_release_date, refseq_only)
        
        api_reports = fetch_virus_metadata(
            virus,
            accession=is_accession,
            host=host,
            geographic_location=geographic_location,
            # geographic_region=geographic_region,
            annotated=api_annotated_filter,
            complete_only=api_complete_filter,
            min_release_date=min_release_date,
            refseq_only=refseq_only,
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
        logger.info("STEP 3: Applying metadata-only filters")
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

        all_none_except_nuc = all(
            v is None for k, v in filters.items() if k != "nuc_completeness"
        )

        if all_none_except_nuc and filters["nuc_completeness"]!="partial":
            logger.info("No metadata-only filters specified, skipping this step.")
            filtered_accessions = list(metadata_dict.keys())
            # logger.debug("All accessions from metadata hereeee: %d", len(filtered_accessions))
            # logger.debug(filtered_accessions)
            filtered_metadata = list(metadata_dict.values())
            # logger.debug("All metadata from metadata hereeee: %d", len(filtered_metadata))
            # logger.debug(filtered_metadata)
            logger.info("All %d sequences will proceed to sequence download", len(filtered_accessions))
        else:
            filtered_accessions, filtered_metadata = filter_metadata_only(metadata_dict, **filters)
            if not filtered_accessions:
                logger.warning("No sequences passed metadata-only filters. Skipping sequence download.")
                return
            # logger.debug("Filtered accessions hereeee!!!: %d", len(filtered_accessions))
            # logger.debug(filtered_accessions)
            # logger.debug("Filtered metadata hereeee!!!: %d", len(filtered_metadata))
            # logger.debug(filtered_metadata)


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
        logger.info("STEP 4: Downloading sequences for filtered accessions")
        logger.info("=" * 60)

        fna_file = download_sequences_by_accessions(filtered_accessions, outdir=temp_dir)
        if not os.path.exists(fna_file):
            raise RuntimeError(f"Download failed: FASTA file not found at {fna_file}")
        logger.info("Downloaded FASTA file: %s (%.2f MB)", fna_file, os.path.getsize(fna_file) / 1024 / 1024)

        # SECTION 6: SEQUENCE-DEPENDENT FILTERING AND SAVING
        logger.info("=" * 60)
        logger.info("STEP 5: Applying sequence-dependent filters and saving results")
        logger.info("=" * 60)

        filters_seq={
            "max_ambiguous_chars": max_ambiguous_chars,
            "has_proteins": has_proteins,
            "proteins_complete": proteins_complete,
        }

        all_none_except_nuc = all(
            v is None for k, v in filters_seq.items() if k != "nuc_completeness"
        )

        if all_none_except_nuc:
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

        # SECTION 7: GENBANK METADATA RETRIEVAL (OPTIONAL)
        logger.info("=" * 60)
        logger.info("STEP 6: Fetching detailed GenBank metadata")
        logger.info("=" * 60)
        logger.info("GenBank metadata retrieval requested - fetching detailed information...")
            
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
                    
                    # Fetch GenBank metadata
                    genbank_data = fetch_genbank_metadata(
                        accessions=list(set(final_accessions)),  # Remove duplicates
                        batch_size=genbank_batch_size,
                        delay=0.5  # Be respectful to NCBI servers
                    )
                    
                    if genbank_data:
                        # Save GenBank metadata to CSV
                        genbank_csv_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata.csv")
                        save_genbank_metadata_to_csv(
                            genbank_metadata=genbank_data,
                            output_file=genbank_csv_path,
                            virus_metadata=filtered_metadata_final
                        )
                        
                        logger.info("✓ GenBank metadata CSV saved: %s (%.2f MB)", 
                                    genbank_csv_path, os.path.getsize(genbank_csv_path) / 1024 / 1024)
                    else:
                        logger.warning("No GenBank metadata was retrieved")
                else:
                    logger.warning("No accession numbers found for GenBank metadata lookup")
                    
            except Exception as genbank_error:
                logger.error("GenBank metadata retrieval failed: %s", genbank_error)
                logger.warning("Continuing without GenBank metadata - standard output files are still available")
                # Don't raise the error - continue with the rest of the process
            
            logger.info("GenBank metadata processing completed")
        else:
            logger.info("GenBank metadata retrieval not requested, skipping this step")

        if filtered_sequences:
            # SECTION 8: FINAL SUMMARY
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

            # Check if GenBank metadata CSV was created
            genbank_csv_path = os.path.join(outfolder, f"{virus_clean}_genbank_metadata.csv")
            if genbank_metadata and os.path.exists(genbank_csv_path):
                logger.info("  🧬 GenBank metadata:  %s", os.path.basename(genbank_csv_path))

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
                virus_clean = virus.replace(' ', '_').replace('/', '_')
                cmd_parts = [f"gget.ncbi_virus('{virus}'"]
                
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
                virus_clean = virus.replace(' ', '_').replace('/', '_')
                cmd_parts = [f"gget.ncbi_virus('{virus}'", "host='human'"]
                
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
                import traceback
                logger.debug("Full traceback:\n%s", traceback.format_exc())
        
        raise
        
    finally:
        # SECTION 9: CLEANUP
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


def fetch_virus_with_genbank_metadata(
    virus,
    accession=False,
    host=None,
    geographic_location=None,
    annotated=None,
    complete_only=None,
    min_release_date=None,
    refseq_only=False,
    outdir=None,
    save_genbank_csv=True,
    genbank_batch_size=200,
    **kwargs
):
    """
    Fetch virus data and detailed GenBank metadata in a single integrated workflow.
    
    This function combines the standard virus metadata retrieval with comprehensive
    GenBank metadata extraction. It provides a complete picture of each virus sequence
    including collection dates, geographic information, host details, publication
    references, and other fields from GenBank records.
    
    Args:
        virus (str): Virus taxon name/ID or accession number
        accession (bool): Whether virus parameter is an accession number
        host (str): Host organism name filter
        geographic_location (str): Geographic location filter
        annotated (bool): Filter for annotated genomes only
        complete_only (bool): Filter for complete genomes only
        min_release_date (str): Minimum release date filter (YYYY-MM-DD format)
        refseq_only (bool): Limit to RefSeq genomes only
        outdir (str): Output directory for files
        save_genbank_csv (bool): Whether to save GenBank metadata to CSV
        genbank_batch_size (int): Batch size for GenBank API requests
        **kwargs: Additional filtering parameters passed to the main virus function
        
    Returns:
        tuple: (filtered_sequences, filtered_metadata, genbank_metadata, output_files)
            - filtered_sequences: List of sequence records that passed filters
            - filtered_metadata: List of virus metadata dictionaries
            - genbank_metadata: Dictionary of detailed GenBank metadata
            - output_files: Dictionary with paths to generated output files
            
    Raises:
        RuntimeError: If virus data retrieval or GenBank fetching fails
        
    Example:
        >>> sequences, metadata, genbank_data, files = fetch_virus_with_genbank_metadata(
        ...     virus="SARS-CoV-2",
        ...     complete_only=True,
        ...     host="human",
        ...     outdir="./output"
        ... )
        >>> print(f"Retrieved {len(sequences)} sequences with GenBank metadata")
        >>> print(f"GenBank CSV saved to: {files['genbank_csv']}")
    """
    
    logger.info("=== INTEGRATED VIRUS + GENBANK METADATA RETRIEVAL ===")
    logger.info("Virus query: %s (accession: %s)", virus, accession)
    logger.info("Filters: host=%s, location=%s, complete=%s, annotated=%s, refseq=%s", 
                host, geographic_location, complete_only, annotated, refseq_only)
    
    # Determine output directory
    if not outdir:
        outdir = os.getcwd()
        logger.info("Using current directory for output: %s", outdir)
    else:
        os.makedirs(outdir, exist_ok=True)
        logger.info("Using specified output directory: %s", outdir)
    
    output_files = {}
    
    try:
        # STEP 1: Fetch standard virus data using the main function
        logger.info("🦠 STEP 1: Fetching virus sequences and metadata...")
        
        # Create a temporary output directory for ncbi_virus function
        virus_temp_dir = os.path.join(outdir, f"temp_{virus.replace(' ', '_')}")
        os.makedirs(virus_temp_dir, exist_ok=True)
        
        # Call the main virus retrieval function with all parameters
        # Note: ncbi_virus saves files but doesn't return data, so we need to read the files
        ncbi_virus(
            virus=virus,
            accession=accession,
            outfolder=virus_temp_dir,  # Note: parameter name is 'outfolder' not 'outdir'
            host=host,
            geographic_location=geographic_location,
            annotated=annotated,
            nuc_completeness='complete' if complete_only else None,
            min_release_date=min_release_date,
            **kwargs
        )
        
        # Read the generated files from ncbi_virus
        virus_name_clean = virus.replace(' ', '_').replace('/', '_').replace('-', '_')
        expected_fasta = None
        expected_csv = None
        expected_jsonl = None
        
        # Look for output files with various naming patterns
        for filename in os.listdir(virus_temp_dir):
            if filename.endswith('.fasta'):
                expected_fasta = os.path.join(virus_temp_dir, filename)
            elif filename.endswith('_metadata.csv'):
                expected_csv = os.path.join(virus_temp_dir, filename)
            elif filename.endswith('_metadata.jsonl'):
                expected_jsonl = os.path.join(virus_temp_dir, filename)
        
        if not expected_fasta or not expected_csv:
            available_files = os.listdir(virus_temp_dir) if os.path.exists(virus_temp_dir) else []
            logger.error("Required output files not found in %s", virus_temp_dir)
            logger.error("Available files: %s", available_files)
            raise RuntimeError("ncbi_virus function did not generate expected output files")
        
        # Read sequences from FASTA file
        logger.info("Reading sequences from: %s", expected_fasta)
        filtered_sequences = list(FastaIO.parse(expected_fasta, "fasta"))
        
        # Read metadata from CSV file (convert back to list of dictionaries)
        logger.info("Reading metadata from: %s", expected_csv)
        import pandas as pd
        metadata_df = pd.read_csv(expected_csv)
        filtered_metadata = metadata_df.to_dict('records')
        
        if not filtered_sequences or not filtered_metadata:
            logger.warning("No virus sequences retrieved - skipping GenBank metadata fetch")
            return None, None, {}, output_files
        
        logger.info("✓ Retrieved %d virus sequences with metadata", len(filtered_sequences))
        
        # STEP 2: Extract accession numbers for GenBank lookup
        logger.info("🔍 STEP 2: Extracting accession numbers for GenBank lookup...")
        
        # Get accession numbers from the filtered sequences
        accession_numbers = []
        for seq_record in filtered_sequences:
            # Extract accession from sequence ID (usually the first part before spaces)
            acc = seq_record.id.split()[0] if hasattr(seq_record, 'id') else str(seq_record)
            if acc:
                accession_numbers.append(acc)
        
        if not accession_numbers:
            logger.warning("No accession numbers found in sequence records")
            return filtered_sequences, filtered_metadata, {}, output_files
        
        logger.info("✓ Found %d unique accession numbers", len(set(accession_numbers)))
        logger.debug("Sample accessions: %s", accession_numbers[:5])
        
        # STEP 3: Fetch detailed GenBank metadata
        logger.info("📚 STEP 3: Fetching detailed GenBank metadata...")
        
        genbank_metadata = fetch_genbank_metadata(
            accessions=list(set(accession_numbers)),  # Remove duplicates
            batch_size=genbank_batch_size,
            delay=0.5  # Be respectful to NCBI servers
        )
        
        logger.info("✓ Retrieved GenBank metadata for %d accessions", len(genbank_metadata))
        
        # STEP 4: Save GenBank metadata to CSV if requested
        if save_genbank_csv and genbank_metadata:
            logger.info("💾 STEP 4: Saving GenBank metadata to CSV...")
            
            genbank_csv_path = os.path.join(outdir, f"{virus}_genbank_metadata.csv")
            save_genbank_metadata_to_csv(
                genbank_metadata=genbank_metadata,
                output_file=genbank_csv_path,
                virus_metadata=filtered_metadata
            )
            
            output_files['genbank_csv'] = genbank_csv_path
            logger.info("✓ GenBank metadata saved to: %s", genbank_csv_path)
        
        # STEP 5: Create summary report
        logger.info("📊 STEP 5: Generating summary report...")
        
        _log_integration_summary(
            virus=virus,
            num_sequences=len(filtered_sequences),
            num_virus_metadata=len(filtered_metadata),
            num_genbank_metadata=len(genbank_metadata),
            output_files=output_files
        )
        
        logger.info("🎉 INTEGRATION COMPLETE: All data successfully retrieved and processed")
        
        return filtered_sequences, filtered_metadata, genbank_metadata, output_files
        
    except Exception as e:
        logger.error("❌ INTEGRATION FAILED: %s", e)
        logger.error("Error occurred during integrated virus + GenBank retrieval")
        
        if logger.getEffectiveLevel() <= logging.DEBUG:
            import traceback
            logger.debug("Full traceback:\n%s", traceback.format_exc())
        
        # Provide helpful troubleshooting guidance
        logger.error("TROUBLESHOOTING SUGGESTIONS:")
        logger.error("1. Check internet connection and NCBI server availability")
        logger.error("2. Verify virus query parameters are correct")
        logger.error("3. Try reducing batch size for GenBank requests")
        logger.error("4. Check output directory permissions")
        
        raise RuntimeError(f"Integrated virus + GenBank metadata retrieval failed: {e}") from e


def _log_integration_summary(virus, num_sequences, num_virus_metadata, num_genbank_metadata, output_files):
    """Log a comprehensive summary of the integration results."""
    
    logger.info("=" * 80)
    logger.info("INTEGRATION SUMMARY REPORT")
    logger.info("=" * 80)
    logger.info("Query: %s", virus)
    logger.info("Sequences retrieved: %d", num_sequences)
    logger.info("Virus metadata records: %d", num_virus_metadata)
    logger.info("GenBank metadata records: %d", num_genbank_metadata)
    
    if output_files:
        logger.info("Output files generated:")
        for file_type, file_path in output_files.items():
            logger.info("  %s: %s", file_type, file_path)
    
    # Calculate data completeness metrics
    if num_genbank_metadata > 0 and num_virus_metadata > 0:
        coverage = min(num_genbank_metadata / num_virus_metadata, 1.0) * 100
        logger.info("GenBank metadata coverage: %.1f%%", coverage)
        
        if coverage < 90:
            logger.info("ℹ️  Note: Some sequences may not have complete GenBank metadata")
            logger.info("   This is normal for newer sequences or specific databases")
    
    logger.info("=" * 80)


# Example usage and integration helper
# def example_genbank_integration():
#     """
#     Example demonstrating how to use the integrated GenBank metadata functionality.
    
#     This function shows the recommended workflow for retrieving both virus sequence
#     data and detailed GenBank metadata in a single operation.
#     """
    
#     # Example 1: SARS-CoV-2 with GenBank metadata
#     print("Example 1: SARS-CoV-2 sequences with GenBank metadata")
#     sequences, virus_meta, genbank_meta, files = fetch_virus_with_genbank_metadata(
#         virus="SARS-CoV-2",
#         complete_only=True,
#         host="human",
#         outdir="./sars_cov2_output",
#         save_genbank_csv=True
#     )
    
#     print(f"Retrieved {len(sequences)} sequences")
#     print(f"GenBank CSV: {files.get('genbank_csv', 'Not saved')}")
    
#     # Example 2: Specific accession with GenBank details
#     print("\nExample 2: Specific accession with GenBank details")
#     sequences, virus_meta, genbank_meta, files = fetch_virus_with_genbank_metadata(
#         virus="NC_045512.2",
#         accession=True,
#         outdir="./reference_output",
#         save_genbank_csv=True
#     )
    
#     # Show sample GenBank metadata
#     if genbank_meta:
#         sample_acc = list(genbank_meta.keys())[0]
#         sample_data = genbank_meta[sample_acc]['genbank_data']
#         print(f"Sample GenBank data for {sample_acc}:")
#         print(f"  Collection date: {sample_data.get('collection_date', 'Not available')}")
#         print(f"  Country: {sample_data.get('country', 'Not available')}")
#         print(f"  Host: {sample_data.get('host', 'Not available')}")
#         print(f"  Reference count: {len(sample_data.get('references', []))}")


if __name__ == "__main__":
    # Main module entry point
    pass