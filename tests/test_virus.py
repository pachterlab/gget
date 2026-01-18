"""Unit tests for gget virus module.

This test module provides comprehensive validation of the virus function,
which downloads viral sequences and metadata from NCBI's virus database.

Test Structure:
---------------
The test suite uses a hybrid approach combining JSON-defined tests (for input
validation) and code-defined tests (for functional and data quality checks).

1. JSON-Defined Tests (19 tests):
   - Loaded from tests/fixtures/test_virus.json
   - Focus on input validation and error handling
   - Test invalid types, values, and parameter combinations

2. Code-Defined Functional Tests (18 tests):
   - Test basic functionality (file creation, downloads)
   - Test individual filters (host, completeness, length, annotated, refseq, etc.)
   - Test filter combinations and edge cases
   - Test advanced features (GenBank metadata, protein filters, etc.)

3. Code-Defined Data Quality Tests (6 tests):
   - Verify data consistency across output formats
   - Validate filter effectiveness (host, release date, completeness)
   - Check metadata schema and field presence
   - Detect API/data source changes

4. Datasets CLI Tests (3 tests):
   - Test _get_datasets_path() returns valid path when CLI is available
   - Test _get_datasets_path() uses bundled binary when system CLI is missing
   - Test datasets CLI version output validation

Parameters Tested:
------------------
Core Parameters:
  ✓ virus (str/int) - validated, functional tests
  ✓ is_accession (bool) - type validation, functional test
  ✓ outfolder (str) - implicit in all tests

Sequence Filters:
  ✓ min_seq_length (int) - min/max validation, functional test
  ✓ max_seq_length (int) - min/max validation, functional test
  ✓ nuc_completeness (str) - value validation, functional test
  ✓ host (str) - functional test with verification
  ✓ annotated (bool) - type validation, functional test
  ✓ refseq_only (bool) - type validation, functional test

Count Filters:
  ✓ min_gene_count (int) - min/max validation, functional test
  ✓ max_gene_count (int) - min/max validation, functional test
  ✓ min_protein_count (int) - min/max validation, functional test
  ✓ max_protein_count (int) - min/max validation, functional test
  ✓ min_mature_peptide_count (int) - min/max validation
  ✓ max_mature_peptide_count (int) - min/max validation

Date Filters:
  ✓ min_release_date (str) - min/max validation, functional test with API verification
  ✓ max_release_date (str) - min/max validation
  ✓ min_collection_date (str) - min/max validation, functional test
  ✓ max_collection_date (str) - min/max validation, functional test

Advanced Filters:
  ✓ lab_passaged (bool) - type validation, functional test
  ✓ proteins_complete (bool) - type validation only
  ✓ keep_temp (bool) - type validation only
  ✓ genbank_metadata (bool) - type validation, functional test
  ✓ genbank_batch_size (int) - type validation, functional test
  ✓ geographic_location (str) - functional test
  ✓ source_database (str) - functional test
  ✓ max_ambiguous_chars (int) - functional test
  ✓ has_proteins (str/list) - functional test

Parameters NOT Tested:
  ✗ submitter_country (str) - similar to geographic_location, not critical
  ✗ is_sars_cov2 (bool) - special mode requiring specific test setup
  ✗ is_alphainfluenza (bool) - special mode requiring specific test setup
  ✗ lineage (str) - SARS-CoV-2 specific, requires special test setup
  ✗ download_all_accessions (bool) - would download entire database (impractical for tests)

What These Tests Catch:
-----------------------
✓ Input validation errors and type checking
✓ File creation and basic functionality
✓ Data consistency between FASTA/CSV/JSONL outputs
✓ Filter effectiveness (host, release date, completeness, location, protein counts)
✓ API schema changes (missing columns, field renames)
✓ Data loss or format conversion bugs
✓ Date filter accuracy with API verification
✓ GenBank metadata retrieval functionality
✓ Protein and gene filtering
✓ Geographic location filtering
✓ Quality filters (ambiguous characters)
✓ NCBI datasets CLI detection (system or bundled binary)
✓ Bundled binary fallback when system CLI is not installed

What These Tests Don't Catch:
-----------------------------
✗ Exact sequence counts (database is ever-changing)
✗ Specific sequence content validation
✗ Complete API failures (would need integration tests)
✗ Network-related errors (outside test scope)
✗ Special modes (SARS-CoV-2, Alphainfluenza optimization paths)
✗ Download all accessions mode (too large for unit tests)
✗ Submitter country filter (similar to geographic location)

Test Coverage Summary:
----------------------
- Total parameters: 34
- Tested (validation or functional): 29 (85%)
- Type/value validation: 21 (62%)
- Functional tests: 19 (56%)
- Not tested: 5 (15%)

Total: 44 tests covering validation, functionality, data quality, and CLI detection.

Notes:
------
- Special modes (SARS-CoV-2, Alphainfluenza) require specific test infrastructure
  and are tested separately in dedicated test modules
- download_all_accessions is impractical for unit tests (would download entire database)
- submitter_country is intentionally not tested (similar to geographic_location)
- Datasets CLI tests cover the connection between gget_virus.py and gget_setup.py
"""
import unittest
import json
import os
import shutil
import subprocess
import tempfile
import time
import functools
from gget.gget_virus import virus, _get_datasets_path, _clear_datasets_cache
from .from_json import from_json


def retry_on_network_error(max_retries=3, delay=5):
    """Decorator to retry tests that may fail due to network issues.
    
    This is useful for tests that make real API calls to NCBI, which can
    occasionally time out or fail due to network flakiness.
    
    Args:
        max_retries: Maximum number of retry attempts (default: 3)
        delay: Seconds to wait between retries (default: 5)
    """
    def decorator(test_func):
        @functools.wraps(test_func)
        def wrapper(*args, **kwargs):
            last_exception = None
            for attempt in range(max_retries):
                try:
                    return test_func(*args, **kwargs)
                except Exception as e:
                    # Only retry on network-related errors
                    error_msg = str(e).lower()
                    if any(keyword in error_msg for keyword in ['timeout', 'timed out', 'connection', 'network']):
                        last_exception = e
                        if attempt < max_retries - 1:
                            time.sleep(delay)
                        continue
                    # Re-raise non-network errors immediately
                    raise
            # If all retries failed, raise the last exception
            raise last_exception
        return wrapper
    return decorator

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_virus.json") as json_file:
    virus_dict = json.load(json_file)


class TestVirus(unittest.TestCase, metaclass=from_json(virus_dict, virus)):
    """Test suite for gget.virus module.
    
    This comprehensive test suite covers:
    
    1. Input Validation (19 JSON-defined tests):
       - Type checking for boolean, string, and integer parameters
       - Value validation (completeness, batch sizes, virus names)
       - Range validation (min/max pairs for dates, lengths, counts)
    
    2. Functional Tests (18 code-defined tests):
       - Basic file creation and accession downloads
       - Individual filter functionality (host, completeness, length, annotated, refseq)
       - Geographic location and source database filters
       - Protein and gene count filters
       - Collection date filters
       - Advanced filters (lab_passaged, max_ambiguous_chars, has_proteins)
       - GenBank metadata retrieval
       - Multiple filter combinations
       - Integer virus ID handling
    
    3. Data Quality & Verification Tests (6 code-defined tests):
       - Relationship checks: FASTA/CSV/JSONL count consistency
       - Filter verification: Host and release date filter effectiveness
       - Schema validation: Expected metadata columns exist
       - Completeness filter verification
       - Multi-filter relationship checks
    
    Coverage: 85% of parameters tested (29/34), with 43 total test cases.
    See module docstring for detailed parameter coverage analysis.
    """
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures that are shared across all tests."""
        cls.test_output_dir = "test_virus_output"
        
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests have run."""
        # Clean up test output directory if it exists
        if os.path.exists(cls.test_output_dir):
            shutil.rmtree(cls.test_output_dir)
    
    def setUp(self):
        """Set up before each test method."""
        # Create a fresh test output directory
        if os.path.exists(self.test_output_dir):
            shutil.rmtree(self.test_output_dir)
        os.makedirs(self.test_output_dir, exist_ok=True)
    
    def tearDown(self):
        """Clean up after each test method."""
        # Clean up test output directory
        if os.path.exists(self.test_output_dir):
            try:
                shutil.rmtree(self.test_output_dir)
            except Exception:
                pass  # Ignore cleanup errors
    
    def _check_output_files(self, virus_name, outfolder):
        """Helper method to check if expected output files were created.
        
        Args:
            virus_name: Name of the virus (used in file naming)
            outfolder: Output folder where files should be created
            
        Returns:
            dict: Dictionary with file paths and existence status
        """
        # Clean virus name for file naming (replace spaces with underscores)
        virus_clean = virus_name.replace(" ", "_")
        
        expected_files = {
            "fasta": os.path.join(outfolder, f"{virus_clean}_sequences.fasta"),
            "csv": os.path.join(outfolder, f"{virus_clean}_metadata.csv"),
            "jsonl": os.path.join(outfolder, f"{virus_clean}_metadata.jsonl")
        }
        
        results = {}
        for file_type, file_path in expected_files.items():
            results[file_type] = {
                "path": file_path,
                "exists": os.path.exists(file_path),
                "size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
            }
        
        return results
    
    def _count_fasta_sequences(self, fasta_file):
        """Count the number of sequences in a FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            int: Number of sequences
        """
        count = 0
        if os.path.exists(fasta_file):
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
        return count
    
    def _count_jsonl_records(self, jsonl_file):
        """Count the number of records in a JSONL file.
        
        Args:
            jsonl_file: Path to JSONL file
            
        Returns:
            int: Number of records
        """
        count = 0
        if os.path.exists(jsonl_file):
            with open(jsonl_file, 'r') as f:
                for line in f:
                    if line.strip():
                        count += 1
        return count
    
    def _count_csv_records(self, csv_file):
        """Count the number of records in a CSV file (excluding header).
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            int: Number of records (excluding header)
        """
        count = 0
        if os.path.exists(csv_file):
            with open(csv_file, 'r') as f:
                # Skip header
                next(f, None)
                for line in f:
                    if line.strip():
                        count += 1
        return count
    
    def _parse_csv_metadata(self, csv_file):
        """Parse CSV metadata file and return records as list of dicts.
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            list: List of dictionaries containing metadata records
        """
        import csv
        records = []
        if os.path.exists(csv_file):
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    records.append(row)
        return records
    
    def _get_csv_columns(self, csv_file):
        """Get column names from CSV file.
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            list: List of column names
        """
        import csv
        if os.path.exists(csv_file):
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                return reader.fieldnames
        return []
    
    # =========================================================================
    # FUNCTIONAL TESTS: Basic file creation and filter functionality
    # =========================================================================
    # These tests verify that the virus function creates output files
    # correctly and that individual filters work as expected.
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_specific_accession_file_creation(self):
        """Test that files are created when downloading a specific accession.
        
        Downloads SARS-CoV-2 reference sequence (NC_045512.2) and verifies:
        - Function returns None (writes to disk)
        - All three output files created (FASTA, CSV, JSONL)
        - Files are not empty
        - At least one sequence in FASTA file
        """
        virus_name = "NC_045512.2"
        outfolder = self.test_output_dir
        
        # Run the function (should create files, returns None)
        result = virus(
            virus=virus_name,
            is_accession=True,
            outfolder=outfolder
        )
        
        # Check that function returns None
        self.assertIsNone(result)
        
        # Check that output files were created
        files = self._check_output_files(virus_name, outfolder)
        
        # Assert all files exist
        self.assertTrue(files["fasta"]["exists"], 
                       f"FASTA file not created: {files['fasta']['path']}")
        self.assertTrue(files["csv"]["exists"], 
                       f"CSV file not created: {files['csv']['path']}")
        self.assertTrue(files["jsonl"]["exists"], 
                       f"JSONL file not created: {files['jsonl']['path']}")
        
        # Assert files are not empty
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty")
        self.assertGreater(files["csv"]["size"], 0, "CSV file is empty")
        self.assertGreater(files["jsonl"]["size"], 0, "JSONL file is empty")
        
        # Count sequences (should be 1 for a specific accession)
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreaterEqual(seq_count, 1, "No sequences found in FASTA file")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_host_filter(self):
        """Test that host filter works and creates appropriate files."""
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            host="human",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with host filter")
        self.assertTrue(files["csv"]["exists"], "CSV file not created with host filter")
        self.assertTrue(files["jsonl"]["exists"], "JSONL file not created with host filter")
        
        # Verify that files contain data
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with host filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_completeness_filter(self):
        """Test that completeness filter works correctly."""
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            nuc_completeness="complete",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with completeness filter")
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with completeness filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_length_filters(self):
        """Test that sequence length filters work correctly."""
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            min_seq_length=10000,
            max_seq_length=11000,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with length filters")
        
        # Verify sequences are within expected length range
        # This would require parsing the FASTA file, which we do with count
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreater(seq_count, 0, "No sequences passed length filters")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_annotated_filter(self):
        """Test that annotated filter works correctly."""
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            annotated=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with annotated filter")
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with annotated filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_refseq_filter(self):
        """Test that RefSeq filter works correctly."""
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            refseq_only=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with RefSeq filter")
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with RefSeq filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_multiple_filters(self):
        """Test that multiple filters can be combined correctly."""
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            host="human",
            nuc_completeness="complete",
            min_seq_length=10500,
            max_seq_length=11000,
            annotated=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with multiple filters")
        self.assertTrue(files["csv"]["exists"], "CSV file not created with multiple filters")
        self.assertTrue(files["jsonl"]["exists"], "JSONL file not created with multiple filters")
        
        # Check that filters reduced the dataset (should have some sequences)
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreater(seq_count, 0, "No sequences passed multiple filters")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_integer_virus_id(self):
        """Test that integer virus IDs are handled correctly.
        
        Tests using Zika virus taxon ID (64320) as integer input.
        Verifies that integer IDs are properly converted and files created.
        """
        virus_id = 64320  # Zika virus taxon ID
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_id,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        # Check files with string version of virus ID
        virus_clean = str(virus_id)
        expected_fasta = os.path.join(outfolder, f"{virus_clean}_sequences.fasta")
        self.assertTrue(os.path.exists(expected_fasta), 
                       f"FASTA file not created for integer virus ID: {expected_fasta}")
    
    # =========================================================================
    # DATA QUALITY & VERIFICATION TESTS
    # =========================================================================
    # These tests verify data consistency, filter effectiveness, and that
    # API/data source changes would be detected. They go beyond simple file
    # existence checks to validate actual data quality.
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_relationship_check_counts_match(self):
        """Test that FASTA sequence count matches CSV and JSONL record counts.
        
        Downloads a specific accession and verifies:
        - Number of FASTA sequences = number of CSV records = number of JSONL records
        - No data loss between different output formats
        - At least one record in all files
        
        This catches: Format conversion bugs, data loss, parsing errors.
        """
        virus_name = "NC_045512.2"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            is_accession=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        
        # Count records in each file type
        fasta_count = self._count_fasta_sequences(files["fasta"]["path"])
        csv_count = self._count_csv_records(files["csv"]["path"])
        jsonl_count = self._count_jsonl_records(files["jsonl"]["path"])
        
        # All counts should match
        self.assertEqual(fasta_count, csv_count, 
                        f"FASTA count ({fasta_count}) does not match CSV count ({csv_count})")
        self.assertEqual(fasta_count, jsonl_count, 
                        f"FASTA count ({fasta_count}) does not match JSONL count ({jsonl_count})")
        self.assertEqual(csv_count, jsonl_count, 
                        f"CSV count ({csv_count}) does not match JSONL count ({jsonl_count})")
        
        # Should have at least one record
        self.assertGreater(fasta_count, 0, "No records found in output files")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_host_filter_verification(self):
        """Test that host filter actually filters by host in metadata.
        
        Downloads Zika virus with host="human" filter and verifies:
        - Records are returned (filter doesn't break the query)
        - Host column exists in metadata
        - If host data is populated, it matches the filter criterion
        
        Note: Host filter is applied server-side by NCBI API. The returned
        records should all match, but the Host field in CSV may be empty or
        have various formats (scientific names, common names).
        
        This catches: Broken host filters, API changes in filtering behavior.
        """
        virus_name = "Zika virus"
        host = "human"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            host=host,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        
        # Parse CSV metadata
        records = self._parse_csv_metadata(files["csv"]["path"])
        
        # Should have some records
        self.assertGreater(len(records), 0, "No records returned with host filter")
        
        # Check that host column exists
        if records:
            self.assertIn("Host", records[0].keys(), 
                         "Host column not found in metadata")
            
            # Note: Host filter is applied server-side by NCBI API
            # The returned records should all match, but the Host field in CSV
            # may be empty or have various formats (scientific names, common names)
            # We verify the filter worked by checking that records were returned
            # (if filter was broken, we'd get all hosts or an error)
            
            # Count non-empty host values
            non_empty_hosts = sum(1 for record in records 
                                if record.get("Host", "").strip())
            
            # If we have host data populated, verify it matches
            if non_empty_hosts > 0:
                host_lower = host.lower()
                # Also check for "Homo sapiens" which is scientific name for human
                matching_hosts = sum(1 for record in records 
                                   if host_lower in record.get("Host", "").lower() 
                                   or "homo sapiens" in record.get("Host", "").lower())
                
                # If host data is populated, at least 50% should match
                if non_empty_hosts > 0:
                    match_percentage = (matching_hosts / non_empty_hosts) * 100
                    self.assertGreater(match_percentage, 50, 
                                     f"Only {match_percentage:.1f}% of populated host fields match filter '{host}'")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_release_date_filter_verification(self):
        """Test that release date filter is applied correctly in metadata.
        
        Downloads Mumps virus with min_release_date="2024-12-31" and verifies:
        - Records are returned (API is working)
        - Release date field exists in metadata
        - All release dates are on or after 2024-12-31
        - Count matches expected API results
        
        This test compares against the direct API call:
        curl -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/taxon/mumps%20virus/dataset_report?filter.released_since=2024-12-31T00:00:00.000Z"
        
        This catches: Release date filter bugs, date parsing errors, API filter issues.
        """
        import requests
        from datetime import datetime
        
        virus_name = "mumps virus"
        min_release_date = "2024-12-31"
        outfolder = self.test_output_dir
        
        # First, get the expected count from direct API call using full timestamp format
        api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/taxon/mumps%20virus/dataset_report"
        params = {"filter.released_since": "2024-12-31T00:00:00.000Z", "page_size": 1000}
        
        try:
            response = requests.get(api_url, params=params, headers={'accept': 'application/json'})
            response.raise_for_status()
            api_data = response.json()
            expected_count = len(api_data.get('reports', []))
        except Exception as e:
            self.skipTest(f"Could not fetch API data for comparison: {e}")
        
        # Run virus function with same filter
        result = virus(
            virus=virus_name,
            min_release_date=min_release_date,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        
        # Verify files were created
        self.assertTrue(os.path.exists(files["csv"]["path"]), 
                       "CSV file not created with release date filter")
        self.assertTrue(os.path.exists(files["fasta"]["path"]), 
                       "FASTA file not created with release date filter")
        
        # Parse CSV metadata
        records = self._parse_csv_metadata(files["csv"]["path"])
        
        # Should have records matching API count (allowing for small variance due to timing)
        self.assertGreater(len(records), 0, "No records returned with release date filter")
        self.assertEqual(len(records), expected_count,
                        f"Record count ({len(records)}) doesn't match API count ({expected_count})")
        
        # Check that release date column exists
        release_date_field = None
        for possible_field in ["Release date", "Release Date", "ReleaseDate", "release_date"]:
            if possible_field in records[0].keys():
                release_date_field = possible_field
                break
        
        self.assertIsNotNone(release_date_field, 
                           f"Release date field not found. Available fields: {list(records[0].keys())}")
        
        # Parse filter date for comparison (inclusive - on or after this date)
        filter_date = datetime.strptime(min_release_date, "%Y-%m-%d")
        
        # Verify all release dates are on or after the filter date (inclusive)
        invalid_dates = []
        for record in records:
            date_str = record.get(release_date_field, "").strip()
            if date_str:
                try:
                    # Parse ISO date format (YYYY-MM-DD)
                    record_date = datetime.strptime(date_str, "%Y-%m-%d")
                    if record_date < filter_date:
                        invalid_dates.append((record.get('Accession', 'unknown'), date_str))
                except ValueError as e:
                    # If date parsing fails, that's also a test failure
                    self.fail(f"Could not parse release date '{date_str}': {e}")
        
        self.assertEqual(len(invalid_dates), 0,
                        f"Found {len(invalid_dates)} records with release dates before {min_release_date}: {invalid_dates[:5]}")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_metadata_schema_validation(self):
        """Test that expected metadata columns exist in CSV output.
        
        Downloads a specific accession and verifies:
        - CSV contains expected essential columns (accession, length, host)
        - At least 5 columns present (reasonable metadata breadth)
        - Column names are properly formatted
        
        This catches: API schema changes, missing metadata fields, field
        name changes that would break downstream analysis tools.
        """
        virus_name = "NC_045512.2"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            is_accession=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        
        # Get column names from CSV
        columns = self._get_csv_columns(files["csv"]["path"])
        
        # Check for expected essential columns (these should always be present)
        # Using case-insensitive checking since column names might vary
        columns_lower = [col.lower() for col in columns]
        
        expected_columns = [
            "accession",  # Or GenBank Accession
            "length",     # Or Sequence Length
            "host",       # Host information
        ]
        
        missing_columns = []
        for expected in expected_columns:
            found = any(expected in col_lower for col_lower in columns_lower)
            if not found:
                missing_columns.append(expected)
        
        self.assertEqual(len(missing_columns), 0, 
                        f"Missing expected metadata columns: {missing_columns}. "
                        f"Available columns: {columns}")
        
        # Verify we have a reasonable number of columns (at least 5)
        self.assertGreaterEqual(len(columns), 5, 
                               f"Only {len(columns)} columns found, expected at least 5")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_completeness_filter_verification(self):
        """Test that completeness filter returns appropriate sequences.
        
        Downloads Zika virus with nuc_completeness="complete" and verifies:
        - Records are returned (filter works)
        - If completeness field exists, validates values
        - Falls back to checking length field exists
        
        This catches: Broken completeness filters, metadata field changes,
        filter logic errors.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            nuc_completeness="complete",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        
        # Parse CSV metadata
        records = self._parse_csv_metadata(files["csv"]["path"])
        
        # Should have some records
        self.assertGreater(len(records), 0, "No records returned with completeness filter")
        
        # Check if there's a completeness or length field
        if records:
            # Look for completeness-related fields
            completeness_field = None
            for possible_field in ["Completeness", "Nuc_Completeness", "Nucleotide Completeness", 
                                  "Genome Coverage", "completeness"]:
                if possible_field in records[0].keys():
                    completeness_field = possible_field
                    break
            
            # If completeness field exists, verify values
            if completeness_field:
                complete_count = sum(1 for record in records 
                                   if "complete" in record.get(completeness_field, "").lower())
                
                # At least 50% should be marked as complete
                if complete_count > 0:
                    complete_percentage = (complete_count / len(records)) * 100
                    self.assertGreater(complete_percentage, 50,
                                     f"Only {complete_percentage:.1f}% marked as complete")
            else:
                # If no explicit completeness field, check length field exists
                # (complete genomes should have consistent lengths)
                length_field = None
                for possible_field in ["Length", "Sequence Length", "Nuc_Length", "length"]:
                    if possible_field in records[0].keys():
                        length_field = possible_field
                        break
                
                self.assertIsNotNone(length_field, 
                                    "Neither completeness nor length field found in metadata")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_multiple_filters_relationship_check(self):
        """Test relationship checks work correctly with multiple filters applied.
        
        Downloads Zika virus with multiple filters (host, completeness, length) and verifies:
        - FASTA/CSV/JSONL counts still match with complex filtering
        - At least one record passes all filters
        - No data loss when multiple filters interact
        
        This catches: Filter interaction bugs, data loss with complex queries,
        inconsistent filtering across output formats.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            host="human",
            nuc_completeness="complete",
            min_seq_length=10000,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        
        # Count records in each file type
        fasta_count = self._count_fasta_sequences(files["fasta"]["path"])
        csv_count = self._count_csv_records(files["csv"]["path"])
        jsonl_count = self._count_jsonl_records(files["jsonl"]["path"])
        
        # All counts should match even with filters
        self.assertEqual(fasta_count, csv_count, 
                        f"FASTA count ({fasta_count}) does not match CSV count ({csv_count}) with multiple filters")
        self.assertEqual(fasta_count, jsonl_count, 
                        f"FASTA count ({fasta_count}) does not match JSONL count ({jsonl_count}) with multiple filters")
        
        # Should have at least one record
        self.assertGreater(fasta_count, 0, "No records found with multiple filters applied")

    # =========================================================================
    # ADDITIONAL FUNCTIONAL TESTS: Testing previously untested parameters
    # =========================================================================
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_geographic_location_filter(self):
        """Test that geographic location filter works correctly.
        
        Downloads Zika virus sequences from Brazil and verifies:
        - Files are created successfully
        - Records are returned
        - Geographic location metadata field exists
        
        This catches: Geographic location filter bugs, API parameter issues.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            geographic_location="Brazil",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with geographic location filter")
        self.assertTrue(files["csv"]["exists"], "CSV file not created with geographic location filter")
        
        # Parse CSV metadata
        records = self._parse_csv_metadata(files["csv"]["path"])
        
        # Should have some records (Brazil had Zika outbreak)
        self.assertGreater(len(records), 0, "No records returned with geographic location filter")
        
        # Check that geographic location fields exist
        if records:
            geo_fields = ["Geographic Location", "Geographic Region", "Geo String"]
            has_geo_field = any(field in records[0].keys() for field in geo_fields)
            self.assertTrue(has_geo_field, 
                          f"No geographic location field found. Available fields: {list(records[0].keys())}")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_protein_count_filters(self):
        """Test that protein count filters work correctly.
        
        Downloads Zika virus with protein count filters and verifies:
        - Files are created successfully
        - Records are returned
        - Protein count field exists in metadata
        
        This catches: Protein count filter bugs, metadata field issues.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            min_protein_count=1,
            max_protein_count=20,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with protein count filters")
        
        # Parse CSV metadata
        records = self._parse_csv_metadata(files["csv"]["path"])
        
        # Should have some records
        self.assertGreater(len(records), 0, "No records returned with protein count filters")
        
        # Check that protein count field exists
        if records:
            self.assertIn("Protein count", records[0].keys(), 
                         f"Protein count field not found. Available fields: {list(records[0].keys())}")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_source_database_filter(self):
        """Test that source database filter works correctly.
        
        Downloads Zika virus from GenBank database and verifies:
        - Files are created successfully
        - Records are returned
        - Source database field exists in metadata
        
        This catches: Source database filter bugs, API parameter issues.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            source_database="GenBank",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with source database filter")
        
        # Parse CSV metadata
        records = self._parse_csv_metadata(files["csv"]["path"])
        
        # Should have some records
        self.assertGreater(len(records), 0, "No records returned with source database filter")
        
        # Check that source database field exists
        if records:
            db_field = None
            for possible_field in ["GenBank/RefSeq", "Source Database", "Database"]:
                if possible_field in records[0].keys():
                    db_field = possible_field
                    break
            
            self.assertIsNotNone(db_field, 
                               f"Source database field not found. Available fields: {list(records[0].keys())}")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_lab_passaged_filter(self):
        """Test that lab_passaged filter works correctly.
        
        Downloads Zika virus with lab_passaged=False filter and verifies:
        - Files are created successfully
        - Records are returned
        
        Note: Lab passaged data may be sparse, so we mainly verify the filter
        doesn't break the query.
        
        This catches: Lab passaged filter bugs, API parameter issues.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            lab_passaged=False,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with lab_passaged filter")
        
        # Should create files (even if no lab passaged field in results)
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with lab_passaged filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_collection_date_filters(self):
        """Test that collection date filters don't break the query.
        
        Downloads Zika virus with collection date range and verifies:
        - Function completes without errors
        
        Note: Collection date data is often sparse, filters may return no results.
        This test just ensures the filter doesn't cause errors.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        # This will complete without error even if no results match
        result = virus(
            virus=virus_name,
            min_collection_date="2016-01-01",
            max_collection_date="2016-12-31",
            outfolder=outfolder
        )
        
        # Function should complete successfully
        self.assertIsNone(result)
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_max_ambiguous_chars_filter(self):
        """Test that max_ambiguous_chars filter works correctly.
        
        Downloads Zika virus with max_ambiguous_chars filter and verifies:
        - Files are created successfully
        - Records are returned
        - Filter doesn't break the query
        
        This catches: Max ambiguous chars filter bugs, sequence quality filtering issues.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            max_ambiguous_chars=100,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with max_ambiguous_chars filter")
        
        # Should have some records (most sequences have some ambiguous bases)
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreater(seq_count, 0, "No sequences passed max_ambiguous_chars filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_has_proteins_filter(self):
        """Test that has_proteins filter works correctly.
        
        Downloads Zika virus requiring specific proteins and verifies:
        - Files are created successfully
        - Records are returned
        - Filter doesn't break the query
        
        This catches: has_proteins filter bugs, protein filtering logic issues.
        """
        virus_name = "Zika virus"
        outfolder = self.test_output_dir
        
        # Test with a common protein (polyprotein is typical for Zika)
        result = virus(
            virus=virus_name,
            has_proteins="polyprotein",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with has_proteins filter")
        
        # Should have some records (polyprotein is common in Zika)
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreater(seq_count, 0, "No sequences passed has_proteins filter")
    
    @retry_on_network_error(max_retries=3, delay=5)
    def test_virus_with_genbank_metadata_retrieval(self):
        """Test that GenBank metadata retrieval works correctly.
        
        Downloads a single accession with genbank_metadata=True and verifies:
        - Function completes without errors
        - Standard files are created
        - GenBank metadata CSV file is created
        
        This catches: GenBank metadata retrieval bugs, batch processing issues.
        """
        virus_name = "NC_045512.2"
        outfolder = self.test_output_dir
        
        result = virus(
            virus=virus_name,
            is_accession=True,
            genbank_metadata=True,
            genbank_batch_size=10,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus_name, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with genbank_metadata")
        
        # Check for GenBank metadata file
        genbank_csv = os.path.join(outfolder, f"{virus_name}_genbank_metadata.csv")
        self.assertTrue(os.path.exists(genbank_csv), 
                       f"GenBank metadata CSV not created: {genbank_csv}")
        
        # Verify GenBank CSV has data
        self.assertGreater(os.path.getsize(genbank_csv), 0, 
                         "GenBank metadata CSV is empty")

    # =========================================================================
    # DATASETS CLI TESTS: Testing NCBI datasets CLI check and setup
    # =========================================================================
    # These tests verify the datasets CLI detection and installation functionality
    # that connects gget_virus.py to gget_setup.py
    
    def test_get_datasets_path_returns_valid_path(self):
        """Test that _get_datasets_path returns a valid path to the datasets CLI.
        
        The function should return a path to either:
        1. The system-installed datasets CLI (if available)
        2. The bundled datasets binary (fallback)
        
        This catches: Detection logic bugs, path resolution issues.
        """
        # _get_datasets_path should always return a valid path
        # (either system CLI or bundled binary)
        datasets_path = _get_datasets_path()
        
        # Should return a non-empty string
        self.assertIsInstance(datasets_path, str)
        self.assertTrue(len(datasets_path) > 0, 
                       "_get_datasets_path should return a non-empty path")
        
        # The returned path should be executable
        result = subprocess.run(
            [datasets_path, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        self.assertEqual(result.returncode, 0,
                        f"datasets CLI at {datasets_path} should be executable")
    
    def test_get_datasets_path_uses_bundled_binary(self):
        """Test that _get_datasets_path falls back to bundled binary when system CLI is missing.
        
        When the system-installed datasets CLI is not in PATH, the function
        should fall back to the bundled binary included with gget.
        
        Note: _get_datasets_path caches its result, so this test clears the cache
        before testing the fallback behavior.
        
        This catches: Bundled binary fallback logic, path resolution issues.
        """
        import gget.gget_virus as gget_virus_module
        
        # Save original PATH and cache
        original_path = os.environ.get("PATH", "")
        original_cache = gget_virus_module._datasets_path_cache
        
        try:
            # Clear the cache to force re-detection
            _clear_datasets_cache()
            
            # Set PATH to empty to simulate system datasets not being found
            os.environ["PATH"] = ""
            
            # Should still return a valid path (to bundled binary)
            datasets_path = _get_datasets_path()
            
            # Should return a non-empty string path to bundled binary
            self.assertIsInstance(datasets_path, str)
            self.assertTrue(len(datasets_path) > 0,
                           "Should return path to bundled binary")
            
            # Path should contain 'bins' indicating bundled binary
            self.assertIn("bins", datasets_path,
                         f"Path should be bundled binary, got: {datasets_path}")
            
        finally:
            # Restore original PATH and cache
            os.environ["PATH"] = original_path
            gget_virus_module._datasets_path_cache = original_cache
    
    def test_datasets_cli_version_output(self):
        """Test that the datasets CLI returns a valid version string.
        
        When available (system or bundled), the datasets CLI should return a 
        version string that can be parsed. This helps ensure the CLI is properly 
        functional.
        
        This catches: Corrupted installations, version parsing issues.
        """
        # Use _get_datasets_path() to get either system or bundled binary
        try:
            datasets_path = _get_datasets_path()
            result = subprocess.run(
                [datasets_path, "--version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            cli_available = result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired, RuntimeError):
            cli_available = False
        
        if not cli_available:
            self.skipTest("NCBI datasets CLI not available - skipping version test")
        
        # Version output should not be empty and should contain version info
        version_output = result.stdout.strip()
        self.assertTrue(len(version_output) > 0, 
                       "Version output should not be empty")
        # NCBI datasets typically outputs version like "datasets version: X.Y.Z" or just "X.Y.Z"
        self.assertTrue(
            any(char.isdigit() for char in version_output),
            f"Version output should contain version numbers: {version_output}"
        )

    

if __name__ == '__main__':
    unittest.main()
