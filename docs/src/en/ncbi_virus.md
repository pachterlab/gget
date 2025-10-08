# NCBI Virus Retrieval Workflow

## Overview

The `gget.ncbi_virus()` function implements an optimized 6-step workflow for retrieving virus sequences and associated metadata from NCBI. The system is designed to minimize download overhead by filtering metadata first, then downloading only the sequences that pass initial filters, with optional detailed GenBank metadata retrieval.

## Architecture

```
┌─────────────────────────────┐
│           Users             │
│                             │
│  • Virus Query (Taxon/Acc)  │
│  • Filter Criteria          │
│    (Host, Dates, Length...) │
│  • Output Flags             │
│    (`--genbank_metadata`)   │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│   API & Pre-Filtering       │
│                             │
│  • Calls NCBI Datasets API  │
│  • Applies server-side      │
│    filters (host, refseq)   │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│ Local Metadata Filtering &  │
│     Sequence Acquisition    │
│                             │
│  • Applies remaining local  │
│    filters (date ranges,    │
│    gene counts, etc.)       │
│  • Generates final list of  │
│    accession numbers        │
│  • Downloads FASTA sequences│
│    via E-utilities API      │
└──────────────┬──────────────┘
               │
   ┌───────────┴──────────────────────────────────────────┐
   │                                                      │
   ▼                                                      ▼
┌─────────────────────────────┐      ┌───────────────────────────────────┐
│   Final Processing          │      │   GenBank Sideload (Optional)     │
│                             │      │                                   │
│  • Applies sequence-level   │      │ • Uses final accession list to    │
│    filters (e.g., max N's)  │      │   fetch detailed GenBank records  │
│  • Formats standard metadata│      │   via E-utilities API             │
└──────────────┬──────────────┘      └──────────────────┬────────────────┘
               │                                        │
               └──────────────────┬─────────────────────┘
                                  │
                                  ▼
                    ┌───────────────────────────────┐
                    │            Results            │
                    │                               │
                    │  • _sequences.fasta           │
                    │  • _metadata.csv & .jsonl     │
                    │  • _genbank_metadata.csv      │
                    └───────────────────────────────┘
```

## Workflow Steps

### Step 1: Input Validation & Setup (Silent)
- **Function**: `ncbi_virus()` main function
- **Purpose**: Validate all user parameters and configure logging
- **Key Operations**:
  - Validate virus taxon/accession format
  - Check filter parameter ranges and formats
  - Set up output directory structure
  - Configure logging based on verbosity level
  - Check for SARS-CoV-2 optimization opportunities

### Step 2: Metadata Retrieval
- **Function**: `fetch_virus_metadata()`
- **Purpose**: Retrieve metadata from NCBI Datasets API with server-side filtering
- **Key Operations**:
  - Apply server-side filters (host, geographic location, release date, completeness)
  - Handle API pagination with connection pooling
  - Implement exponential backoff with jitter for retries
  - Parse JSON responses with streaming for large datasets
  - Store metadata in structured format with validation

### Step 3: Metadata-Only Filtering
- **Function**: `filter_metadata_only()`
- **Purpose**: Apply local filters that don't require sequence data
- **Key Operations**:
  - Filter by date ranges with smart date parsing
  - Filter by genome completeness and quality indicators
  - Apply numeric range filters (gene/protein counts)
  - Handle missing or malformed metadata gracefully
  - Generate optimized accession list for targeted download

### Step 4: Targeted Sequence Download
- **Function**: `download_sequences_by_accessions()`
- **Purpose**: Download FASTA sequences only for filtered accessions
- **Key Operations**:
  - Use E-utilities API with batch optimization
  - Implement configurable batch sizes (default: 200)
  - Stream large responses to manage memory
  - Handle download retries with exponential backoff
  - Return path to downloaded FASTA file

### Step 5: GenBank Metadata Retrieval (Optional)
- **Function**: `fetch_genbank_metadata()`
- **Purpose**: Fetch detailed GenBank records for sequences
- **Key Operations**:
  - Retrieve comprehensive GenBank records
  - Extract 23+ metadata fields per record
  - Process in configurable batch sizes
  - Implement rate limiting and retries
  - Parse and validate GenBank XML

### Step 6: Sequence-Dependent Filtering & Output
- **Function**: `filter_sequences()`
- **Purpose**: Apply final filters requiring sequence analysis and save results
- **Key Operations**:
  - Parse FASTA sequences and calculate sequence metrics
  - Filter by sequence length, ambiguous character count
  - Filter by protein completeness indicators
  - Save filtered sequences and metadata to output files

## Function Dependencies

```
ncbi_virus()
├── is_sars_cov2_query()
│   └── SARS-CoV-2 detection logic
├── download_sars_cov2_optimized()  [For SARS-CoV-2 queries]
│   ├── NCBI datasets CLI calls
│   └── Cached package downloads
├── fetch_virus_metadata()
│   ├── NCBI Datasets API client
│   ├── Connection pooling
│   ├── Pagination handling
│   └── Retry logic with backoff
├── filter_metadata_only()
│   ├── Date parsing utilities
│   ├── Numeric validation
│   └── Missing data handling
├── download_sequences_by_accessions()
│   ├── E-utilities API client
│   ├── Batch processing
│   └── Stream handling
├── fetch_genbank_metadata() [Optional]
│   ├── E-utilities API client
│   ├── XML parsing utilities
│   ├── Rate limiting
│   └── Batch processing
├── filter_sequences()
│   ├── FastaIO parser
│   ├── Sequence validation
│   └── Output generation
└── unzip_file()  [For cached downloads]
    └── ZIP extraction utilities
```

## Optimization Features

### 1. **Server-Side Filtering**
- Applies filters at the NCBI API level to reduce data transfer
- Supported filters: host, geographic location, release date, genome completeness
- Automatic validation of filter compatibility and values

### 2. **Multi-Stage Filtering**
- **Stage 1**: Metadata-only filters (fast, no sequence download)
- **Stage 2**: Sequence-dependent filters (pre-filtered set)
- **Stage 3**: GenBank metadata integration and filtering
- **Stage 4**: Final validation and quality checks

### 3. **Optimized Downloads**
- Configurable batch sizes for different data types
- Connection pooling for improved performance
- Stream handling for large downloads
- Rate limiting and retry mechanisms

### 4. **SARS-CoV-2 Optimization**
- Special handling for SARS-CoV-2 queries using NCBI's cached data packages
- Automatic detection and fallback to standard API
- Optimized GenBank metadata retrieval for SARS-CoV-2

### 5. **Efficient Data Structures**
- Accession-based dictionaries for O(1) lookups
- Streaming parsers for JSON and XML
- Memory-efficient FASTA handling
- Optimized metadata merging

## Output Files

### 1. **FASTA Sequences** (`{virus}_sequences.fasta`)
- Contains nucleotide sequences for filtered results
- Standard FASTA format with detailed headers
- Original orientation from NCBI preserved
- Optional protein/segment annotations in headers

### 2. **CSV Metadata** (`{virus}_metadata.csv`)
- Tabular format for spreadsheet analysis
- Standardized column structure
- Geographic and taxonomic information
- Collection and submission details
- Quality metrics and annotations

### 3. **GenBank Metadata** (`{virus}_genbank_metadata.csv`) [Optional]
- 23+ detailed metadata columns
- Publication references
- Feature annotations
- Cross-references to other databases
- Strain and isolate details

### 4. **JSONL Metadata** (`{virus}_metadata.jsonl`)
- JSON Lines format for full data structure
- Streaming-friendly format
- Complete nested data preserved (protein and annotated data)
- Combined virus and GenBank data when available

## Performance Characteristics

### Scalability
- **Small datasets** (< 1,000 sequences): Near-instantaneous processing
- **Medium datasets** (1,000 - 10,000 sequences): Minutes to complete
- **Large datasets** (> 10,000 sequences): Optimized pagination and filtering

### Memory Usage
- Streaming processing minimizes memory footprint
- Metadata cached in memory for filtering operations
- Large FASTA files processed in chunks

### Network Efficiency
- Minimal API calls due to server-side filtering
- Targeted downloads reduce bandwidth usage
- Automatic retry with exponential backoff

## Error Handling

### API Failures
- Smart retry strategy with exponential backoff and jitter
- Server-side error detection with specific guidance:
  - Timeout handling for large datasets
  - Geographic filter optimization suggestions
  - Batch size adjustments for GenBank metadata
- Connection pooling and session management
- Detailed error logging with troubleshooting steps

### Data Validation
- Comprehensive input parameter validation:
  - Type checking for all parameters
  - Range validation for numeric values
  - Date format and range validation
  - Boolean parameter normalization
- Sequence integrity verification:
  - FASTA format validation
  - Ambiguous character detection
  - Protein/gene completeness checks
- Metadata consistency validation:
  - Required field presence checks
  - Data type validation
  - Cross-reference validation
  - GenBank record validation

### Recovery Mechanisms
- Automatic temporary file cleanup
- Partial result preservation:
  - Intermediate metadata saving
  - Progressive filtering state saving
  - GenBank metadata caching
- Hierarchical fallback strategies:
  - SARS-CoV-2 optimized packages
  - Cached data fallback
  - API-based retrieval fallback
- Detailed error reporting:
  - Root cause analysis
  - Alternative command suggestions
  - Filter relaxation recommendations
  - Performance optimization tips

## Usage Examples

### Command Line Examples

```bash
# Get help and see all available parameters
$ gget ncbi_virus --help

$ gget ncbi_virus "Nipah virus"

# Download Zika virus sequences with basic filtering (API + metadata filtering)
$ gget ncbi_virus "Zika virus" --host human --min_seq_length 10000 --max_seq_length 11000

# Download with metadata and sequence filtering
$ gget ncbi_virus "Ebolavirus" --max_seq_length 20000 --genbank_metadata -o ./ebola_data

# Download SARS-CoV-2 with cached optimization
$ gget ncbi_virus "SARS-CoV-2" --host dog --nuc_completeness complete

# Download Influenza A with post-download sequence filtering (warning: big data size)
$ gget ncbi_virus "Influenza A virus" --host human --max_ambiguous_chars 50 --has_proteins spike

# Using accession ID to get data
$ gget ncbi_virus -a "MK947457" --host deer --min_collection_date "2020-01-01"
```

### Python Examples

```python
import gget

# Basic download with GenBank metadata
result = gget.ncbi_virus(
    "Zika virus",
    host="human",
    genbank_metadata=True,
    outfolder="zika_data"
)

# Access different data types
sequences = result['sequences']           # List of SeqRecord objects
virus_metadata = result['metadata']       # List of virus metadata dicts
genbank_data = result['genbank_metadata'] # Dict of GenBank records

# Print GenBank metadata summary
for acc, data in genbank_data.items():
    print(f"Sequence: {acc}")
    print(f"  Length: {data['sequence_length']} bp")
    print(f"  Host: {data.get('host', 'Unknown')}")
    print(f"  Location: {data.get('country', 'Unknown')}")
    print(f"  Collection date: {data.get('collection_date', 'Unknown')}")

# Advanced filtering with GenBank data
result = gget.ncbi_virus(
    "SARS-CoV-2", 
    host="human",
    min_seq_length=29000,
    max_seq_length=30000,
    min_collection_date="2020-03-01",
    max_collection_date="2020-03-31",
    geographic_region="North America",
    genbank_metadata=True,
    genbank_batch_size=200,
    outfolder="covid_march2020"
)

# Process and analyze results
import pandas as pd

# Read virus metadata
virus_df = pd.read_csv("covid_march2020/SARS-CoV-2_metadata.csv")
print(f"Total sequences: {len(virus_df)}")
print(f"Unique hosts: {virus_df['Host'].nunique()}")
print(f"Date range: {virus_df['Collection Date'].min()} to {virus_df['Collection Date'].max()}")

# Read GenBank metadata for detailed analysis
genbank_df = pd.read_csv("covid_march2020/SARS-CoV-2_genbank_metadata.csv")
print(f"Sequences with GenBank data: {len(genbank_df)}")
print("
Publication summary:")
print(genbank_df['reference_count'].describe())

# Custom filtering workflow
from Bio import SeqIO

# Read sequences
sequences = list(SeqIO.parse("covid_march2020/SARS-CoV-2_sequences.fasta", "fasta"))

# Custom sequence analysis
for record in sequences:
    gc_content = (str(record.seq).count('G') + str(record.seq).count('C')) / len(record.seq)
    print(f"{record.id}: GC content = {gc_content:.2%}")

# Merge metadata sources
merged_df = pd.merge(
    virus_df,
    genbank_df,
    on='accession',
    how='left',
    suffixes=('_virus', '_genbank')
)

# Save merged analysis
merged_df.to_csv("covid_march2020/combined_analysis.csv", index=False)
```

### Analysis Strategy Examples

The examples above demonstrate different analysis approaches:

1. **Basic GenBank Integration**: Fetch sequences with GenBank metadata for comprehensive analysis
2. **Advanced Filtering**: Combine virus metadata and GenBank data with custom filters
3. **Custom Analysis**: Process sequences and metadata using BioPython and Pandas
4. **Data Integration**: Merge virus and GenBank metadata for detailed analysis

### Programmatic Access
```python
# Access filtered metadata and sequences
metadata_file = "covid_data/SARS-CoV-2_metadata.jsonl"
sequences_file = "covid_data/SARS-CoV-2_sequences.fasta"

# Process results with custom analysis
import json
with open(metadata_file) as f:
    for line in f:
        record = json.loads(line)
        # Custom analysis here
```

