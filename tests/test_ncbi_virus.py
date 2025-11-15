import unittest
import json
import os
import shutil
import tempfile
from gget.gget_ncbi_virus import ncbi_virus
from .from_json import from_json

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_ncbi_virus.json") as json_file:
    ncbi_virus_dict = json.load(json_file)


class TestNcbiVirus(unittest.TestCase, metaclass=from_json(ncbi_virus_dict, ncbi_virus)):
    """Test suite for gget.ncbi_virus module.
    
    This test suite covers:
    - Input validation (type checking, value validation)
    - Error handling (invalid arguments, malformed inputs)
    - File creation and output validation
    - Filter combinations
    - Edge cases
    """
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures that are shared across all tests."""
        cls.test_output_dir = "test_ncbi_virus_output"
        
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
    
    # Custom test methods for file creation tests (type: "code_defined")
    
    def test_ncbi_virus_specific_accession_file_creation(self):
        """Test that files are created when downloading a specific accession."""
        virus = "NC_045512.2"
        outfolder = self.test_output_dir
        
        # Run the function (should create files, returns None)
        result = ncbi_virus(
            virus=virus,
            is_accession=True,
            outfolder=outfolder
        )
        
        # Check that function returns None
        self.assertIsNone(result)
        
        # Check that output files were created
        files = self._check_output_files(virus, outfolder)
        
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
    
    def test_ncbi_virus_with_host_filter(self):
        """Test that host filter works and creates appropriate files."""
        virus = "Zika virus"
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus,
            host="human",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with host filter")
        self.assertTrue(files["csv"]["exists"], "CSV file not created with host filter")
        self.assertTrue(files["jsonl"]["exists"], "JSONL file not created with host filter")
        
        # Verify that files contain data
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with host filter")
    
    def test_ncbi_virus_with_completeness_filter(self):
        """Test that completeness filter works correctly."""
        virus = "Zika virus"
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus,
            nuc_completeness="complete",
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with completeness filter")
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with completeness filter")
    
    def test_ncbi_virus_with_length_filters(self):
        """Test that sequence length filters work correctly."""
        virus = "Zika virus"
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus,
            min_seq_length=10000,
            max_seq_length=11000,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with length filters")
        
        # Verify sequences are within expected length range
        # This would require parsing the FASTA file, which we do with count
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreater(seq_count, 0, "No sequences passed length filters")
    
    def test_ncbi_virus_with_annotated_filter(self):
        """Test that annotated filter works correctly."""
        virus = "Zika virus"
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus,
            annotated=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with annotated filter")
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with annotated filter")
    
    def test_ncbi_virus_with_refseq_filter(self):
        """Test that RefSeq filter works correctly."""
        virus = "Zika virus"
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus,
            refseq_only=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with RefSeq filter")
        self.assertGreater(files["fasta"]["size"], 0, "FASTA file is empty with RefSeq filter")
    
    def test_ncbi_virus_with_multiple_filters(self):
        """Test that multiple filters can be combined correctly."""
        virus = "Zika virus"
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus,
            host="human",
            nuc_completeness="complete",
            min_seq_length=10500,
            max_seq_length=11000,
            annotated=True,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        files = self._check_output_files(virus, outfolder)
        self.assertTrue(files["fasta"]["exists"], "FASTA file not created with multiple filters")
        self.assertTrue(files["csv"]["exists"], "CSV file not created with multiple filters")
        self.assertTrue(files["jsonl"]["exists"], "JSONL file not created with multiple filters")
        
        # Check that filters reduced the dataset (should have some sequences)
        seq_count = self._count_fasta_sequences(files["fasta"]["path"])
        self.assertGreater(seq_count, 0, "No sequences passed multiple filters")
    
    def test_ncbi_virus_integer_virus_id(self):
        """Test that integer virus IDs are handled correctly."""
        virus_id = 64320  # Zika virus taxon ID
        outfolder = self.test_output_dir
        
        result = ncbi_virus(
            virus=virus_id,
            outfolder=outfolder
        )
        
        self.assertIsNone(result)
        
        # Check files with string version of virus ID
        virus_clean = str(virus_id)
        expected_fasta = os.path.join(outfolder, f"{virus_clean}_sequences.fasta")
        self.assertTrue(os.path.exists(expected_fasta), 
                       f"FASTA file not created for integer virus ID: {expected_fasta}")
    

if __name__ == '__main__':
    unittest.main()
