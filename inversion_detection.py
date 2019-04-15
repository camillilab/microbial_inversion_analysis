"""
Playground for downloading GenBank files, wrapping in bowtie2, sCLIP/SOR prep
"""

from Bio import Entrez
import csv
import re
import requests
import subprocess
import os
from Bio import SeqIO
import shutil
from urllib.error import HTTPError
import time
import inv_config as config

# TODO Remove .gbk files at the end that don't have SOR/sCLIP to save disc space?

# my ncbi info - please don't share haha. You can get your own easily by logging in to NCBI and requesting an API Key
Entrez.email = config.email
Entrez.api_key = config.api_key


def ascp(accession_num, save_path=os.getcwd()):
    url = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=read_run&fields=fastq_ftp'.format(
        accession_num)

    # Send the request to ENA
    print("Requesting ENA for FASTQ FTP link for accession run {0}...".format(accession_num))

    r = requests.get(url)

    # from the text of the request, grab the fastq link(s)
    fastq_finder = re.compile('ftp.*?fastq.gz')

    print("FASTQ FTP links found:")
    fastq_links = fastq_finder.findall(r.text)

    if len(fastq_links) < 2:
        print("Insufficient links found! Please check accession number or if the accession has submitted FASTQ files.")
        return False

    for link in fastq_links:
        print(link)

    # Alright, now for each link, build an ascp command

    # Modify as needed, but should be default
    ascp_openssh_file = config.ascp_ssh_key

    print("Retrieving files by ascp...")
    for link in fastq_links:
        # build the ASCP file path
        ascp_path = 'era-fasp@fasp.sra.ebi.ac.uk:/' + link[18:]

        # build the ASCP command
        cmd = 'ascp -QT -l300M -P33001 -i "{0}" {1} {2}'.format(ascp_openssh_file, ascp_path, save_path)

        # subprocess
        try:
            subprocess.run(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            print("Error:", err)
    return True


# uses bowtie2-build to make a reference index
def bowtie2_build(ref, ind):
    subprocess.run(['bowtie2-build', ref, ind], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return


# align using bowtie2
# note bowtie2 uses stderr for output, oddly enough
# Use local mode to capture clipping events
def bowtie2(ind, fq1, fq2, sam, use_threads=4):

    # Very rarely, there is a poorly formatted FASTQ file that catches. Return a fail.
    try:
        subprocess.run(['bowtie2', '-x', ind, '-1', fq1, '-2', fq2, '-S', sam, '-p', str(use_threads), '--local'],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        return True
    except subprocess.CalledProcessError as e:
        print("Bowtie2 error! {0}".format(e))
        return False


def get_file_len(a):
    with open(a, 'r') as b:
        lines = len(b.readlines())
        return lines


# get SOR files using awk
def dump_sor(reads, outfile):
    headers = 'colors,POS,CIGAR,TLEN\n'

    # CMD has been modified to exclude zero-based TLEN
    cmd = "awk 'BEGIN { OFS = \",\" } $2 ~ /113|177|65|129/ $9 !~ /0/ {print $2, $4, $6, $9}'"
    with open(outfile, 'w') as o:
        # add headers
        o.write(headers)
        o.flush()
        subprocess.run(cmd, stdin=reads, shell=True, stdout=o)
    return


# get sCLIP reads using awk
def dump_sclip(reads, outfile):
    headers = 'colors,POS,CIGAR,TLEN\n'
    cmd = "awk 'BEGIN {OFS=\",\"} ($2 ~ /147|83/ && $6 ~ /^..?S/) || ($2 ~ /99|163/ && $6 ~ /S$/) {next;} $6 ~ /^..?S/ {print $2, $4, $6, $9 }'"
    with open(outfile, 'w') as o:
        o.write(headers)
        o.flush()
        # maybe using rb will allow us to properly read samtools bam
        subprocess.run(cmd, stdin=reads, shell=True, stdout=o)
    return


# bam, sort, and index using samtools
def bamify(sam_file, bam_file, use_threads=8):

    # first, convert sam to bam
    print("Convering to BAM...")
    subprocess.run(['samtools', 'view', '-u', '-b', sam_file, '-o', 'tmp.bam', '-@', str(use_threads)])

    # then, sort the bam file
    print("Sorting...")
    subprocess.run(['samtools', 'sort', 'tmp.bam', '-o', bam_file, '-@', str(use_threads)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # now, index the bam file
    print("Indexing...")
    subprocess.run(['samtools', 'index', bam_file])

    os.remove('tmp.bam')

    return


# Extract reads by accession
# Keep in memory to pass to awk commands
def extract_reads(acc, bam_file, sor_file, sclip_file, use_threads=8):

    print("Extracting reads from {0}...".format(acc))
    # Unfortunately, assumes .1 for accession version...could alternatively re-do the grab so it maintains version
    with open('tmp.sam', 'w') as o:
        subprocess.run(['samtools', 'view', bam_file, acc+'.1', '-@', str(use_threads)], stdout=o, encoding='utf-8')

    with open('tmp.sam', 'r') as i:
        #print("Extracting SOR reads...")
        dump_sor(i, sor_file)

    with open('tmp.sam', 'r') as i:
        #print("Extracting sCLIP reads...")
        dump_sclip(i, sclip_file)

    os.remove('tmp.sam')

    return


### MAIN ###
# Firstly, load up the data table
data_file = config.data_file
acc_save_path = config.acc_save_path
fasta_save_path = config.fasta_save_path
run_save_path = config.run_save_path
sam_save_path = config.sam_save_path
sor_save_path = config.sor_save_path
sclip_save_path = config.sclip_save_path
script_path = config.script_path
acc_list_path = config.acc_list_path
use_threads = config.use_threads
sor_read_threshold = config.sor_read_threshold
sclip_read_threshold = 100  # Usually not a problem.
max_error = 10  # Timeout error threshold


with open(data_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')

    for row in reader:

        my_paths = [acc_save_path, fasta_save_path, run_save_path, sam_save_path, sor_save_path, sclip_save_path]
        for path in my_paths:
            if not os.path.exists(path):
                os.mkdir(path)

        acc_with_genes = []  # List of accessions we want to run SOR mapping for
        good_acc = []  # List of accession numbers we ultimately want to pass to detect_inversion_clusters

        # Now, we want a couple things: The SRA Accession, the Biosample Accession, and the list of NUCCORE Accessions
        sra_accession = row['SRA']
        biosample_accession = row['Biosample']
        nuccore_accessions = row['RefSeq Accessions'].split(',')

        print("Now processing: {0}".format(biosample_accession))

        # First, we need to bowtie2 the reads and reference together. So let's first grab the reads:
        code = ascp(sra_accession, save_path=run_save_path)
        if not code:
            print("ASCP error! Skipping...")
        else:

            # Then, let's grab the accessions using Entrez Efetch - we want separate gbk files and a combined FASTA
            # Sometimes there are literally too many accessions, and we get an HTTP 414 error. Break it up by 50?

            batch_size = 50
            for i in range(0, len(nuccore_accessions), batch_size):
                end = min(len(nuccore_accessions), i + batch_size)
                print('Retrieving gbk records {0} to {1}...'.format(i, end))

                current_acc = nuccore_accessions[i:end]
                nuccore_acc_query = ','.join(current_acc)

                num_attempts = 1

                while num_attempts < max_error:
                    try:
                        handle = Entrez.efetch(db='nuccore', id=nuccore_acc_query, rettype='gbwithparts', retmode='text')
                        num_attempts = max_error + 1
                    except HTTPError as err:
                        if 500 <= err.code <= 599:
                            print("Received error from server: {0}".format(err.code))
                            print("Attempt {0} of {1}".format(num_attempts, max_error))
                            num_attempts += 1
                            time.sleep(15)
                        else:
                            raise
                    for record in SeqIO.parse(handle, format='gb'):
                        name = os.path.join(acc_save_path, record.name + '.gb')

                        # Verify CDS info. If none, exclude from further analysis
                        elements = record.features
                        num_cds = 0
                        for element in elements:
                            if element.type == "CDS":
                                num_cds += 1

                        if num_cds == 0:
                            print("No gene data detected for {0}. Removing from analysis...".format(record.name))
                        else:
                            acc_with_genes.append(record.name)

                        with open(name, 'w') as out_handle:
                            SeqIO.write(record, out_handle, "gb")
                    handle.close()
                print("{0} accessions with gene data detected.".format(len(acc_with_genes)))

                print("Retrieving fasta records {0} to {1}...".format(i, end))
                fasta_output = os.path.join(fasta_save_path, biosample_accession + '.fasta')
                fasta_records = []
                with Entrez.efetch(db='nuccore', id=nuccore_acc_query, rettype='fasta', retmode='text') as handle:
                    for record in SeqIO.parse(handle, format='fasta'):
                        fasta_records.append(record)

            # Now write all the fasta to a single combined reference record.
            with open(fasta_output, 'w') as out_handle:
                SeqIO.write(fasta_records, out_handle, "fasta")

            # Now let's use bowtie2 to align the reads
            ref_path = os.path.join(fasta_save_path, biosample_accession + '.fasta')
            f1 = os.path.join(run_save_path, sra_accession + '_1.fastq.gz')
            f2 = os.path.join(run_save_path, sra_accession + '_2.fastq.gz')
            sam_output = os.path.join(sam_save_path, biosample_accession + '.sam')

            print("Aligning {0} to read set {1} using bowtie2...".format(ref_path, sra_accession))
            bowtie2_build(ref_path, ind='INDEX')
            code = bowtie2(ind='INDEX', fq1=f1, fq2=f2, sam=sam_output, use_threads=use_threads)

            if not code:
                print("Bowtie2 encountered an error! Skipping {0}...".format(biosample_accession))
            else:
                bam_output = os.path.join(sam_save_path, biosample_accession + '.bam')
                # print("Indexing and sorting using SAMtools...")
                bamify(sam_output, bam_output, use_threads=use_threads)


                # Now, for each accession, let's extract the reads
                for acc in acc_with_genes:

                    sor_file = os.path.join(sor_save_path, acc + '_sor.csv')
                    sclip_file = os.path.join(sclip_save_path, acc + '_sclip.csv')

                    # Extract reads for the accession
                    extract_reads(acc, bam_output, sor_file, sclip_file, use_threads=use_threads)

                    # Check to make sure the SOR and sCLIP files aren't empty. For ones that aren't add to accessions_list.txt
                    sor_lines = get_file_len(sor_file)
                    sclip_lines = get_file_len(sclip_file)
                    if (sor_lines >= sor_read_threshold) and (sclip_lines >= sclip_read_threshold):
                        good_acc.append(acc)
                    else:
                        print("{0} has insufficient SOR/sCLIP reads! Excluding from analysis. (SOR={1}, sCLIP={2})".format(acc, sor_lines, sclip_lines))

                # Create a list of accessions with actual SOR/sCLIP data to feed to detect_inversion_clusters
                with open(acc_list_path, 'w') as acc_list:
                    for acc in good_acc:
                        acc_list.write(acc+'\n')

                # Now, with everything in place, run the detect inversions script, placing output in a Biosample folder
                # Also check to make sure SOR and sCLIP have data in them first!

                print("Executing detection script...")
                subprocess.run(['python3', '{0}'.format(script_path), biosample_accession])

                # Now remove all the data we no longer need to save hard disk space.
                print("Cleaning SAM/BAM files for {0}...".format(biosample_accession))
                shutil.rmtree(sam_save_path, ignore_errors=True)
                print("Cleaning Run files for {0}...".format(sra_accession))
                shutil.rmtree(run_save_path, ignore_errors=True)
                print("Cleaning sCLIP files for {0}...".format(biosample_accession))
                shutil.rmtree(sclip_save_path, ignore_errors=True)
                print("Cleaning SOR files for {0}...".format(biosample_accession))
                shutil.rmtree(sor_save_path, ignore_errors=True)
                print("Cleaning Entrez FASTA files for {0}...".format(biosample_accession))
                shutil.rmtree(fasta_save_path, ignore_errors=True)
                print("Done!")


