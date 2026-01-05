#!/bin/bash

# Test commands for gget ncbi_virus with different scenarios
# Logs are saved to corresponding .log files

uv pip install .

echo "Starting gget ncbi_virus test runs..."
echo "======================================"
echo ""

# Test 1: COVID-19 variants from California, March 2020
echo "Test 1: COVID-19 from California (March 2020)..."
gget ncbi_virus "SARS-CoV-2" --geographic_location "USA: California" --min_collection_date "2020-03-01" --max_collection_date "2020-03-31" -o covid_ca_march2020 > covid_ca_march2020.log 2>&1
echo "✓ Test 1 complete. Log: covid_ca_march2020.log"
echo ""

# Test 2: Alphainfluenza downloads
echo "Test 2: Alphainfluenza download..."
gget ncbi_virus "Alphainfluenza" -g -o alphainfluenza_download > alphainfluenza_download.log 2>&1
echo "✓ Test 2 complete. Log: alphainfluenza_download.log"
echo ""

# Download all human viruses released in 2024 with complete genomes
echo "Download all human viruses released in 2024 with complete genomes..."
gget ncbi_virus --download_all_accessions --host human --nuc_completeness complete --min_release_date 2024-01-01 --min_seq_length 300 -o all_viruses_2024 > all_viruses_2024.log 2>&1
echo "✓ Smoke Test complete. Log: all_viruses_2024.log"
echo ""

# Test 3: Download ALL viruses (WARNING: EXTREMELY large - entire NCBI Virus database)
# Uncomment the following lines to run this test (use with caution!)
echo "Test 3: Download ALL viruses (WARNING: This is HUGE!)..."
gget ncbi_virus --download_all_accessions --refseq_only --nuc_completeness complete -o all_viruses > all_viruses_download.log 2>&1
echo "✓ Test 3 complete. Log: all_viruses_download.log"
echo ""

# Test 4: COVID-19 with filters and GenBank metadata from USA 2020
echo "Test 4: COVID-19 USA 2020 with complete genomes and GenBank metadata..."
gget ncbi_virus "SARS-CoV-2" --geographic_location "USA" --min_collection_date "2020-01-01" --max_collection_date "2020-12-31" -g --nuc_completeness complete -o covid_usa_2020_complete > covid_usa_2020_complete.log 2>&1
echo "✓ Test 4 complete. Log: covid_usa_2020_complete.log"
echo ""

# Test 5: Alphainfluenza with specific filters (USA only, with GenBank metadata)
echo "Test 5: Alphainfluenza from USA with GenBank metadata..."
gget ncbi_virus "Alphainfluenza" --geographic_location "USA" --min_seq_length "2000" --max_ambiguous_chars "10" -g -o alphainfluenza_usa_2000 > alphainfluenza_usa_2000.log 2>&1
echo "✓ Test 5 complete. Log: alphainfluenza_usa_2000.log"
echo ""

# Test 6: Small COVID-19 test from Wuhan origin period
echo "Test 6: COVID-19 from Wuhan (Dec 2019 - Jan 2020)..."
gget ncbi_virus "SARS-CoV-2" --geographic_location "China: Wuhan" --min_collection_date "2019-12-01" --max_collection_date "2020-01-31" --has_proteins "spike" -g -o covid_wuhan_origin > covid_wuhan_origin.log 2>&1
echo "✓ Test 6 complete. Log: covid_wuhan_origin.log"
echo ""

echo "======================================"
echo "All tests complete!"
echo ""
echo "To view logs in real-time during a run, use:"
echo "  tail -f <logfile>.log"
echo ""
echo "To check download results:"
echo "  ls -lh */"
echo "  wc -l */*.csv"
echo "  grep -c '>' */*.fasta"
