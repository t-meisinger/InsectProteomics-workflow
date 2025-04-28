#######################################################################################################################################################
# Proteomic insights into Novel Food Insects: Homology-based proteome characterization and allergenicity considerations for EU-regulated insect species 
# Tobias Meisinger, Hannes Planatscher, Albert Braeuning, Eva-Maria Ladenburger, Dieter Stoll, Cristiano Garino, Hermann Broll, Oliver Poetz
# 2025

# Script for preparation of the arthropod reference database: Removal of clusters that consist only of UniParc entries. 
# See section "Databases" of the manuscript.
# Instructions: 1. Install all required packages, as defined in the import section
#               2. Run from CMD (Windows) or bash (Linux): 
#                   python fasta_remove_UniParc.py -f <example.fasta>
#######################################################################################################################################################

# Imports
import argparse
from Bio import SeqIO
import re

def process_fasta(input_file, output_UniParc_file, output_cleaned_file):
    # Regular expression to match 'RepID=UPI' followed by any characters
    repid_pattern = re.compile(r'RepID=UPI\w*')

    # Open the input and output files
    with open(input_file, 'r') as infile, \
         open(output_UniParc_file, 'w') as outfile_UniParc, \
         open(output_cleaned_file, 'w') as outfile_cleaned:
        
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(infile, 'fasta'):
            header = record.description
            
            # Check if the header contains 'RepID=UPI' and 'n=1'
            if 'n=1' in header and repid_pattern.search(header):
                # Write the record to the UniParc output file
                SeqIO.write(record, outfile_UniParc, 'fasta')
            else:
                # Write the record to the cleaned output file
                SeqIO.write(record, outfile_cleaned, 'fasta')

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Process a multi-FASTA file to separate entries based on specific criteria.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Input multi-FASTA file')
    args = parser.parse_args()
    
    # Define output file names
    input_file = args.file
    output_UniParc_file = 'UniParc_entries.fasta'
    output_cleaned_file = 'cleaned_entries.fasta'
    
    # Call the function to process the FASTA files
    process_fasta(input_file, output_UniParc_file, output_cleaned_file)

if __name__ == '__main__':
    main()
