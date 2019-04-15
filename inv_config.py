"""
Configuration file for inv_repeat_analysis

Author: Jake Bourgeois
Email: jacob.bourgeois@tufts.edu
Affiliation: Camilli Laboratory
License: MIT

"""
import os

# NCBI Info - needed for scraping information without too much consequence from NCBI
email = 'your_email@goes.here'
api_key = 'ncbi_key'

# inversion_detection.py variables
data_file = 'ncbi_data.txt'  # Output table from ncbi_scrape.py

# Various file locations for data from inversion_detection.py. Assumes CWD. May change if desired to external drive.
acc_save_path = os.path.join(os.getcwd(), 'Entrez')
fasta_save_path = os.path.join(os.getcwd(), 'FASTA')
run_save_path = os.path.join(os.getcwd(), 'SRA')
sam_save_path = os.path.join(os.getcwd(), 'SAM')
sor_save_path = os.path.join(os.getcwd(), 'SOR Data')
sclip_save_path = os.path.join(os.getcwd(), 'sCLIP Data')
script_path = os.path.join(os.getcwd(), 'detect_inversion_clusters.py')
acc_list_path = os.path.join(os.getcwd(), 'accession_list.txt')
cluster_save_path = os.path.join(os.getcwd(), 'cluster_data')

# Required amount of non-zero point SOR reads to include in final analysis.
sor_read_threshold = 50


# Analyze_clusters info

# One - switches to be considered for database management. These are switches you already know exist that you may
# compare novel IR sequences to. The CDIF switches are listed here. Modify as needed.
switches = [
    ('TTNTAATTCTAAAGGNTACTT', 'cdi1'),
    ('GTTGTAAAANNGTT', 'cdi2'),
    ('CATTTCTTGTAAAATGGATAGTTT', 'cdi3'),
    ('AAGTTNCTATNTTACAAAAAA', 'cdi4'),
    ('GTTGTAAAANNGTTA', 'cdi5'),
    ('TTAAAGTNNCCNNTTGTTGGANAANGGA', 'cdi6'),
    ('TTAACTTTTGA', 'cdi7')
    ]

# Two - you can compare novel proteins to a reference database. For our example, we have R20291 as a model Cdif
# Place a protein fasta file to make a blastp database from
ref_protein_database = 'r20291.txt'

# Number of nt to look around inversion sites for more IR sequences if none detected the first time around
ir_nt_buffer = 50


# Other variables

use_threads = 8  # amount of threads to use with any multithreaded command, such as bowtie or SAMTools.

# Loaction of the ascp client ssh key. This is the default installation location. Change if necessary.
ascp_ssh_key = os.path.expanduser('~/Applications/Aspera CLI/etc/asperaweb_id_dsa.openssh')
