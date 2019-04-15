"""
Script for biosample-level analysis of SOR/sCLIP inversions

Outputs the following:
- General statistics (total number samples, total number inversions, distinct inverted repeats, etc.
- Rough histogram of protein functions
  - For hypothetical proteins, use identical proteins group Entrez database to group with a putative pfam
- Rough histogram of detected inverted repeats, or "Nothing Found"
  - Consider running an extra IR detection script of nothing found inversions - complex
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import os
import csv
from collections import defaultdict
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline as bioblastp
from Bio.Blast.Applications import NcbiblastnCommandline as bioblastn
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio import SearchIO
import re
from urllib.error import HTTPError
import time
import inv_config

# NCBI info for IPG/efetch
Entrez.email = inv_config.email
Entrez.api_key = inv_config.api_key

# ------ANALYSES-------


# Needs modification to assess an overall alignment - there are SNPs in many
# Maybe make a local BLAST db?
def create_switch_inversion_blastdb():

    my_switches = inv_config.switches
    my_seqs = []

    for switch, name in my_switches:
        switch_seqrecord = SeqRecord(seq=Seq(switch, generic_nucleotide), id=name)
        my_seqs.append(switch_seqrecord)

    with open('switches.fasta', 'w') as o:
        SeqIO.write(my_seqs, o, format='fasta')
    make_local_blast_db('switches.fasta', dbname='SWITCH', dbtype='nucl')
    return


def append_known_inversions(ir):

    # BLAST against the db
    ir_seq = SeqRecord(seq=Seq(ir, generic_nucleotide).upper(), id='IRSEQ')
    with open('ir.fasta', 'w') as o:
        SeqIO.write(ir_seq, o, 'fasta')
    blast('ir.fasta', dbname='SWITCH', out_xml='ir.xml', prog='blastn')
    blast_result = SearchIO.read('ir.xml', 'blast-xml')
    best_hit = 'None detected'
    if blast_result.hits:
        bh = blast_result.hits[0]
        hsp = bh.hsps[0]
        best_hit_name = bh.blast_id
        evalue = hsp.evalue
        hit_strand = hsp.hit_strand
        hit_range = hsp.hit_range

        mod = ''
        if hit_strand == -1:
            mod = '(reverse strand)'

        best_hit = '{0}: {1} eval:{3} {2}'.format(best_hit_name, hit_range, mod, evalue)

        # print(blast_result.hits[0])
        # print(blast_result.hsps[0])

    return best_hit


# Create a local blast database of provided FASTA.
# Saves to /usr/local/ncbi/blast/db/
def make_local_blast_db(fasta, dbname, dbtype='nucl'):
    # dbname_path = os.path.join('/usr/local/ncbi/blast/db', dbname)
    subprocess.run(['makeblastdb', '-in', fasta, '-parse_seqids', '-title', dbname, '-dbtype', dbtype, '-out', dbname],
                   stdout=subprocess.DEVNULL)
    return


# Blasts a sequence over a local database.
def blast(seq, dbname, out_xml, prog):

    # Run the blastp command and save the xml
    if prog == 'blastp':
        blast_prog = bioblastp(query=seq, db=dbname, outfmt=5, out=out_xml, num_threads=4)
        _, _ = blast_prog()
    elif prog == 'blastn':
        blast_prog = bioblastn(query=seq, db=dbname, outfmt=5, out=out_xml, num_threads=4, task='blastn-short')
        _, _ = blast_prog()
    else:
        print("Program {0} not recognized.".format(prog))
    return


# Analysis: make frequency histogram of inverted repeats detected.
def inverted_repeat_analysis(biosamples, analysis_file):

    # Create the switch database
    create_switch_inversion_blastdb()

    # Make a defaultdict of type int
    inverted_repeat_frequency_dict = defaultdict(int)

    # Make a defaultdict of type list for uniquq product descriptors
    inverted_repeat_upstream_dict = defaultdict(list)
    inverted_repeat_downstream_dict = defaultdict(list)
    inverted_repeat_within_dict = defaultdict(list)
    inverted_repeat_enc_dict = defaultdict(list)

    # Cycle through the inversions
    for i in biosamples:
        for k in i.accessions():
            for j in k.inversions:
                inverted_repeat = j.ir_seq
                nearby_genes = j.nearby_genes

                inverted_repeat_frequency_dict[inverted_repeat] += 1

                for p in nearby_genes:
                    if p == 'within':
                        if nearby_genes['within'].qualifiers['product'][0] not in inverted_repeat_within_dict[inverted_repeat]:
                            inverted_repeat_within_dict[inverted_repeat].append(nearby_genes['within'].qualifiers['product'][0])
                    if p == 'upstream':
                        if nearby_genes['upstream'].qualifiers['product'][0] not in inverted_repeat_upstream_dict[inverted_repeat]:
                            inverted_repeat_upstream_dict[inverted_repeat].append(nearby_genes['upstream'].qualifiers['product'][0])
                    if p == 'downstream':
                        if nearby_genes['downstream'].qualifiers['product'][0] not in inverted_repeat_downstream_dict[inverted_repeat]:
                            inverted_repeat_downstream_dict[inverted_repeat].append(nearby_genes['downstream'].qualifiers['product'][0])
                    if p == 'encompass':
                        if nearby_genes['encompass'].qualifiers['product'][0] not in inverted_repeat_enc_dict[inverted_repeat]:
                            inverted_repeat_enc_dict[inverted_repeat].append(nearby_genes['encompass'].qualifiers['product'][0])

    with open(analysis_file, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(("Inverted Repeat", "Known Cdi", "Number detected", "Upstream Descriptors",
                         "Downstream Descriptors", "Within Descriptors", "Encompass Descriptors"))
        for ir in inverted_repeat_frequency_dict:
            known_ir = append_known_inversions(ir)
            # known_ir = 'skip'
            writer.writerow((ir, known_ir, inverted_repeat_frequency_dict[ir], inverted_repeat_upstream_dict[ir],
                            inverted_repeat_downstream_dict[ir], inverted_repeat_within_dict[ir],
                             inverted_repeat_enc_dict[ir]))

    return


# Analysis: make descriptor:IR repeat(s) if an inverted repeat is associated
def descriptor_analysis(biosamples, descriptor, descriptor_file):

    # Make a blastdb for reference
    master_ref_prot_file = inv_config.ref_protein_database
    db_name = 'REF_PROT'
    make_local_blast_db(fasta=master_ref_prot_file, dbname=db_name, dbtype='prot')

    # First search for all instances of inversions with the descriptor
    biosamples_with_descriptor = list()
    for b in biosamples:
        for a in b.accessions():
            for i in a.inversions:
                for g in i.nearby_genes:
                    cds = i.nearby_genes[g]
                    if descriptor in cds.qualifiers['product'][0]:

                        try:
                            # BLASTp the cds to find the R20291 equivalent
                            cds_protein_name = cds.qualifiers['locus_tag'][0]
                            cds_protein_translation = cds.qualifiers['translation'][0]
                            cds_protein = SeqRecord(seq=Seq(cds_protein_translation, generic_protein), id=cds_protein_name,
                                                    description=cds.qualifiers['product'][0])
                            with open('cds.fasta', 'w') as out_handle:
                                SeqIO.write(cds_protein, out_handle, 'fasta')
                            blast('cds.fasta', dbname=db_name, out_xml='cds.xml', prog='blastp')

                            blast_result = SearchIO.read('cds.xml', 'blast-xml')
                            # Grab the first hit and grab the R20291 locus tag
                            r20291_hit = blast_result[0]
                            r20291_locus_tag = re.search('(?<=locus_tag=).*?(?=\])', r20291_hit.description).group(0)

                        except KeyError:
                            # No translation likely available. Skip.
                            r20291_locus_tag = 'No translation of protein.'

                        biosamples_with_descriptor.append((b.accession,
                                                           b.name,
                                                           b.strain,
                                                           a.accession,
                                                           i.pos_start,
                                                           i.pos_end,
                                                           i.ir_seq,
                                                           append_known_inversions(i.ir_seq),
                                                           g,
                                                           cds.qualifiers['locus_tag'][0],
                                                           cds.qualifiers['product'][0],
                                                           r20291_locus_tag))

    with open(descriptor_file, 'w') as o:
        writer = csv.writer(o)
        writer.writerow(('Biosample Accession',
                         'Biosample Name',
                         'Biosample Strain',
                         'RefSeq Accession',
                         'Inversion Start',
                         'Inversion End',
                         'Inverted Repeat',
                         'Annotated Repeat',
                         'Relative position to Inversion',
                         'Locus Tag',
                         'Description',
                         'Reference Locus Tag BLASTp Hit'))
        for b in biosamples_with_descriptor:
            writer.writerow(b)

    return

# Analysis:

# ----CLASSES AND PROCESSING----


# Match scoring for IR finder
def match_score(s, match, mismatch):

    score = 0
    for c in s:
        if c.isupper():
            score += match
        else:
            score += mismatch

    return score


# Function parses the analysis folder output into pythonic classes
def process_analysis_folder(analysis_folder_path, entrez_data_path, download=False):

    print("Obtaining Cluster data from {0}...".format(analysis_folder_path))
    my_samples = []    # List of biosamples to return

    # Get a list of the biosample folders
    with os.scandir(analysis_folder_path) as biosamples:
        for biosample in biosamples:
            if not biosample.name.startswith('.'):
                with os.scandir(biosample) as sample_files:
                    my_accessions = []
                    for accession in sample_files:
                        if not accession.name.startswith('.'):
                            with os.scandir(accession) as accession_files:
                                    for accession_file in accession_files:
                                        if accession_file.name[-14:] == 'inversions.csv':
                                            my_inversions = []
                                            with open(accession_file, 'r') as inv_file:
                                                inv_data = inv_file.readlines()
                                                try:
                                                    invs = inv_data[6:]
                                                    for i in invs:
                                                        parts = cluster_parse(i)
                                                        my_inversions.append(Inversion(*parts,
                                                                                       parent=accession.name))
                                                except IndexError:
                                                    pass

                                            # Now collect the data into Accession class
                                            a = Accession(accession=accession.name, parent=biosample.name,
                                                          children=my_inversions)
                                            my_accessions.append(a)

                    # Now collect the data into Biosample class
                    b = Biosample(accession=biosample.name, children=my_accessions)
                    my_samples.append(b)

    # Now let's fill in the blanks using Entrez

    # Fill in biosample data
    print("Querying for biosample metadata...")
    fill_in_biosample_data(my_samples)

    # If we have the Entrez folder still, we can use that. Otherwise, use Entrez again, then fill in the data
    if download:
        print("Querying for nuccore data, saving to {0}...".format(entrez_data_path))
        download_nuccore_entrez_data(my_samples, entrez_data_path)
    print("Obtaining nuccore data from {0}...".format(entrez_data_path))
    fill_in_nuccore_data_by_files(my_samples, entrez_data_path)

    return my_samples


# parses the line in the inversion analysis file - see process_analysis_folder
def cluster_parse(line):
    data = line.split(',')
    return [int(data[1]), int(data[3]), data[2], data[4], data[6]]


# Function uses Entrez NCBI to get specifics for each biosample
def fill_in_biosample_data(biosamples):

    # First, we can use Biosample Esummary to update the biosample

    # Build a list of all biosample accessions to make the call to esearch to get all the uids
    biosample_accessions = []
    for biosample in biosamples:
        biosample_accessions.append(biosample.accession)

    # Now let's make the Esearch call - Let's post all the relevant biosample uids to the history server
    search_accessions = [k+'[accn]' for k in biosample_accessions]
    search_query = ' OR '.join(search_accessions)

    with Entrez.esearch(db='biosample', term=search_query, usehistory='y') as esearch_handle:
        esearch_data = Entrez.read(esearch_handle)
    esearch_webenv = esearch_data['WebEnv']
    esearch_querykey = esearch_data['QueryKey']

    # Now let's use Esummary to get a dictionary for name and strain
    biosample_data_dict = dict()
    with Entrez.esummary(db='biosample', query_key=esearch_querykey, webenv=esearch_webenv) as esummary_handle:
        esummary_data = Entrez.read(esummary_handle)

        for record in esummary_data['DocumentSummarySet']['DocumentSummary']:
            b_id = record.attributes['uid']
            acc = record['Accession']
            org = record['Organism']
            strain = record['Infraspecies']

            biosample_data_dict[acc] = (b_id, org, strain)

    # Finally, let's append the data to the biosamples
    for biosample in biosamples:
        biosample.esummary_update(*biosample_data_dict[biosample.accession])

    return


# Fills in nuccore data by Entrez requests
# This takes entirely too long to loop...need to either multiprocess it or use sets. Or both.
def fill_in_nuccore_data_by_entrez(biosamples):

    # Alright, now let's move onto the nuccore accessions
    # Make a list of all accessions we've seen
    all_accessions = []
    accession_container = []
    biosample_nuccore_accessions = [j.accessions() for j in biosamples]
    for b in biosample_nuccore_accessions:
        all_accessions += [k.accession for k in b]
        accession_container += b

    search_query = ','.join(all_accessions)
    with Entrez.esearch(db='nuccore', term=search_query, usehistory='y') as esearch_handle:
        esearch_data = Entrez.read(esearch_handle)
    esearch_webenv = esearch_data['WebEnv']
    esearch_querykey = esearch_data['QueryKey']

    # Let's retrieve records 20 at a time
    batch_size = 20

    for x in range(0, len(all_accessions), batch_size):
        y = min(len(all_accessions), x + batch_size)

        print("Annotating gbk info for {0} to {1} of {2}...".format(x, y, len(all_accessions)))
        # Send the query

        max_attempts = 10
        num_attempts = 1
        while num_attempts <= max_attempts:
            print("Attempt {0} of {1}".format(num_attempts, max_attempts))
            try:
                efetch_handle = Entrez.efetch(db='nuccore', retstart=x, retmax=batch_size, rettype='gbwithparts', retmode='text',
                           webenv=esearch_webenv, query_key=esearch_querykey)
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server: {0}".format(err))
                    print("Attempt {0} of {1}".format(num_attempts, max_attempts))
                    num_attempts += 1
                    time.sleep(15)
                else:
                    raise

        efetch_data = SeqIO.parse(efetch_handle, format='gb')
        for record in efetch_data:

            name = record.name
            seq = record.seq
            desc = record.description

            cds = [k for k in record.features if k.type == 'CDS']
            payload = [desc, seq, cds]

            # Now, find out accession. If it's ours, append it as necessary
            for a in accession_container:
                if a.accession == name:
                    a.add_genbank(*payload)

                    # Remove from container to make next loop shorter
                    accession_container.remove(a)

    return


# If we have the Entrez gb files still, we can use those instead of making Entrez requests
# This takes entirely too long to loop...need to multiprocess.
def fill_in_nuccore_data_by_files(biosamples, entrez_data_path):

    all_accessions = []
    accession_container = []
    biosample_nuccore_accessions = [j.accessions() for j in biosamples]
    for b in biosample_nuccore_accessions:
        all_accessions += [k.accession for k in b]
        accession_container += b

    entrez_files = [k for k in os.listdir(entrez_data_path) if k[:-3] in all_accessions]

    # Accession container is smaller
    for entrez_file in entrez_files:
        org_name = entrez_file[:-3]
        org_path = os.path.join(entrez_data_path, entrez_file)

        # A while loop will be faster: can jump out when we find the match

        i = 0
        a = accession_container[i]
        while org_name != a.accession:
            i += 1
            a = accession_container[i]
        # print("Found {0}".format(a.accession))
        with open(org_path, 'r') as efile:
            for record in SeqIO.parse(efile, format='gb'):

                seq = record.seq
                desc = record.description
                cds = [k for k in record.features if k.type == 'CDS']
                payload = [desc, seq, cds]
                a.add_genbank(*payload)
                a.data_path = org_path

        for inversion in a.inversions:
            inversion.set_nuccore_accession(a.accession)

        accession_container.remove(a)

    return


# Download entrez records - that way, you can use fill_in_nuccore_by_files for subsequent analyses.
def download_nuccore_entrez_data(biosamples, entrez_data_path):
    # Alright, now let's move onto the nuccore accessions
    # Make a list of all accessions we've seen
    all_accessions = []
    accession_container = []
    biosample_nuccore_accessions = [j.accessions() for j in biosamples]
    for b in biosample_nuccore_accessions:
        all_accessions += [k.accession for k in b]
        accession_container += b

    search_query = ','.join(all_accessions)
    with Entrez.esearch(db='nuccore', term=search_query, usehistory='y') as esearch_handle:
        esearch_data = Entrez.read(esearch_handle)
    esearch_webenv = esearch_data['WebEnv']
    esearch_querykey = esearch_data['QueryKey']

    # Let's retrieve records 20 at a time
    batch_size = 20
    max_attempts = 10

    for x in range(0, len(all_accessions), batch_size):
        y = min(len(all_accessions), x + batch_size)

        print("Downloading gbk for {0} to {1} of {2} records...".format(x, y, len(all_accessions)))
        num_attempts = 1
        while num_attempts <= max_attempts:
            num_attempts += 1
            try:
                efetch_handle = Entrez.efetch(db='nuccore', retstart=x, retmax=batch_size, rettype='gbwithparts', retmode='text',
                           webenv=esearch_webenv, query_key=esearch_querykey)
                num_attempts = max_attempts + 1
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print('Network error {0}'.format(err))
                    time.sleep(15)
                else:
                    raise
        efetch_data = SeqIO.parse(efetch_handle, format='gb')
        for record in efetch_data:
            name = os.path.join(entrez_data_path, record.name + '.gb')
            with open(name, 'w') as out_handle:
                SeqIO.write(record, out_handle, "gb")
    return


# Fills in hypothetical gene data using IPG
def annotate_by_ipg(protein_accessions, inv_map_dict):
    # First, build an esearch call to get accessions to uids
    search_accessions = [k + '[accn]' for k in protein_accessions]
    search_query = ' OR '.join(search_accessions)

    with Entrez.esearch(db='protein', term=search_query, usehistory='y') as esearch_handle:
        esearch_data = Entrez.read(esearch_handle)
    esearch_webenv = esearch_data['WebEnv']
    esearch_querykey = esearch_data['QueryKey']

    # Now let's use Esummary to get a dictionary for name and strain
    with Entrez.efetch(db='protein', retmode='xml', rettype='ipg', query_key=esearch_querykey,
                       webenv=esearch_webenv) as efetch_handle:
        efetch_data = Entrez.read(efetch_handle)

        for record in efetch_data:
            ipg_list = record['ProteinList'][0]
            query_protein_acc = record.attributes['product_acc']  # use this to connect back to our inversion instance
            ipg_functions = defaultdict(int)  # Record number of occurrences we see a description
            for ipg in ipg_list:
                ipg_fxn = ipg.attributes['name']
                if 'hypothetical' not in ipg_fxn and 'Uncharacterized' not in ipg_fxn and 'Uncharacterised' not in ipg_fxn and 'uncharacterized'not in ipg_fxn and 'uncharacterised' not in ipg_fxn:
                    ipg_functions[ipg_fxn] += 1

            # Find the best IPG description
            if ipg_functions:
                best_ipg = max(ipg_functions, key=lambda key: ipg_functions[key])

                # add to inversion instance description, overwrite
                to_modify = inv_map_dict[query_protein_acc]
                for inv, loc in to_modify:
                    inv.nearby_genes[loc].qualifiers['product'] = ['IPG: '+best_ipg]


# class container for Biosample
class Biosample:

    def __init__(self, name=None, strain=None, accession=None, children=None):

        # Core parameters
        self.accession = accession
        self.name = name
        self.strain = strain

        # Linking parameters
        self.nuccore_accessions = children
        self.uid = None

    # Method to add nuccore accessions to Biosample
    def add_nuccore_accessions(self, acc):
        for n in acc:
            self.nuccore_accessions.append(n)
        return

    # Method to return accessions:
    def accessions(self):
        return self.nuccore_accessions

    # Method using Esummary to update the instance
    # Pass the parsed handle
    def esummary_update(self, uid, org, strain):
        self.uid = uid
        self.name = org
        self.strain = strain
        return


# class container for Nuccore Accessions
class Accession:

    def __init__(self, name=None, accession=None, parent=None, children=None):

        # Core parameters
        self.name = name
        self.accession = accession
        self.strain = None

        # Parameters from loading a gbkwithparts file
        self.data_path = None   # Path to gbk file
        self.sequence = None    # Organism sequence
        self.cds = []           # List of coding sequences

        # Linking parameters
        self.biosample_accession = parent   # Biosample accession
        self.inversions = children          # Associated inversions

        return

    # Method to add inversions
    def add_inversions(self, invs):
        for i in invs:
            self.inversions.append(i)
        return

    # Method to append genbank file info
    def add_genbank(self, strain, sequence, cds):
        self.strain = strain
        self.sequence = sequence
        self.cds = cds
        return

    # Finds the closest upstream and downstream loci of the inversion. If applicable, also detect if inside a gene or
    # encompassing a gene.
    def detect_nearest_genes(self):

        # First see if we have gene data. Else, raise.
        if not self.inversions:
            return
        if not self.cds:
            raise AttributeError('No gene data detected on Accession instance {0}!'.format(self.accession))

        # For each inversion:
        for i in self.inversions:

            # Inversion attributes
            inv_start = i.pos_start
            inv_end = i.pos_end

            # First, let's get a list of CDS that's on the positive strand
            pos_strand_cds = [k for k in self.cds if k.strand == 1]
            neg_strand_cds = [k for k in self.cds if k.strand == -1]

            # Define our vars
            closest_genes = dict()

            # Iterate through the positive gene tracts and identify the closest upstream and downstream genes
            d_up_min = 9999999
            d_down_min = 9999999

            # A for loop is not ideal, but sometimes we reach the end of upstream/downstream without fulfilling the position condition in a while loop
            # TODO list comprehension conditions may fix this
            for gene in pos_strand_cds:

                cds_start = gene.location.start
                cds_end = gene.location.end
                locus_tag = gene.qualifiers['locus_tag'][0]
                product = gene.qualifiers['product'][0]
                # transl = pos_strand_cds[j].qualifiers['translation'][0]
                # protein_id = pos_strand_cds[j].qualifiers['protein_id'][0]

                # Case one - we are upstream of the inversion. Here, the inversion start is greater than the end of the CDS
                if inv_start >= cds_end:
                    d_up = inv_start - cds_end
                    if d_up < d_up_min:
                        closest_genes['upstream'] = gene
                        d_up_min = d_up

                # Case two - the inversion encompasses the cds
                if inv_start <= cds_start and inv_end >= cds_end:
                    closest_genes['encompass'] = gene

                # Case three - the inversion is within the cds
                if inv_start >= cds_start and inv_end <= cds_end:
                    closest_genes['within'] = gene

                # Case four - we are downstream of the inversion
                if inv_end <= cds_start:
                    d_down = cds_start - inv_end
                    if d_down < d_down_min:
                        closest_genes['downstream'] = gene
                        d_down_min = d_down

            # Repeat for genes in the negative orientation. Here, we consider the boundaries differently.
            for gene in neg_strand_cds:

                cds_start = gene.location.start
                cds_end = gene.location.end
                locus_tag = gene.qualifiers['locus_tag'][0]
                product = gene.qualifiers['product'][0]
                # transl = pos_strand_cds[j].qualifiers['translation'][0]
                # protein_id = pos_strand_cds[j].qualifiers['protein_id'][0]

                # Case one - we are upstream of the inversion
                if inv_end <= cds_start:
                    d_up = cds_start - inv_end
                    if d_up < d_up_min:
                        closest_genes['upstream'] = gene
                        d_up_min = d_up

                # Case two - the inversion encompasses the cds
                if inv_start <= cds_start and inv_end >= cds_end:
                    closest_genes['encompass'] = gene

                # Case three - the inversion is within the cds
                if (inv_start >= cds_start and inv_end <= cds_end)\
                        or (inv_start <= cds_start <= inv_end <= cds_end) \
                        or (cds_start <= inv_start <= cds_end <= inv_end):
                    closest_genes['within'] = gene

                # Case four - we are downstream of the inversion
                if inv_start >= cds_end:
                    d_down = inv_start - cds_end
                    if d_down < d_down_min:
                        closest_genes['downstream'] = gene
                        d_down_min = d_down

            # Now apply the data to the inversion object
            i.nearby_genes = closest_genes

        return

    # Processing: if the inverted repeat is "Nothing Found", try harder. Consider EMBOSS IR finder.
    # EMBOSS IR finder will not work. Unable to set minimum inversion distance.
    # Ultimately, implementation is optimal in the original inversion detection. But here's a shot!
    def find_inverted_repeats(self, buffer=50, threshold=8, match=3, mismatch=-4, max_error=3):

        if not self.inversions:
            return None

        q = 0
        # for i in [a for a in self.inversions if a.ir_seq == 'Nothing found']:  # Alternative version where only those without previous annotations are described.
        for i in self.inversions:

            q += 1  # So we can report how many inversions were appended.

            ir_dict = dict()  # Container for all IRs found

            # Obtain the inversion positions
            inv_start = i.pos_start
            begin = max(inv_start - buffer, 0)

            inv_end = i.pos_end
            finish = min(inv_end + buffer, len(self.sequence))

            # Start looking for repeats upstream of buffer nt
            for left_pos in range(begin, inv_start + buffer):

                # Build a seed sequence of threshold nt long and it's RC
                left_seq = self.sequence[left_pos:left_pos + threshold]

                # Now iterate through at the end
                for right_pos in range(inv_end - buffer, finish):

                    # obtain the downstream seq and it's reverse complement
                    right_seq = self.sequence[right_pos:right_pos + threshold]
                    right_seq_rc = right_seq.reverse_complement()

                    # Do these strings match?
                    if left_seq == right_seq_rc:

                        # If so, try to continue building the sequence. Continue until we hit max errors
                        e, k = 0, 1
                        final_left_seq = left_seq
                        mismatch_check_seq = left_seq
                        while e <= max_error:

                            left_seq = self.sequence[left_pos: left_pos + threshold + k]
                            # For right_seq, look upstream for more effectors.
                            right_seq = self.sequence[right_pos - k: right_pos + threshold]
                            right_seq_rc = right_seq.reverse_complement()

                            # Have -1 to index so it grabs the right letter!!
                            if left_seq[-1] != right_seq_rc[-1]:
                                e += 1
                                final_left_seq += self.sequence[left_pos + threshold + k - 1].lower()
                            else:
                                final_left_seq += self.sequence[left_pos + threshold + k - 1]
                                if e == 0:
                                    mismatch_check_seq += self.sequence[left_pos + threshold + k - 1]
                            k += 1

                        # Remove any mismatches off the end and convert IRseq to string

                        final_ir_seq, mismatch_check_seq = str(final_left_seq), str(mismatch_check_seq)
                        while final_ir_seq[-1].islower():
                            final_ir_seq = final_ir_seq[:-1]

                        no_mismatch_score = match_score(mismatch_check_seq, match, mismatch)
                        all_score = match_score(final_ir_seq, match, mismatch)

                        # Now, add the better score to the ir_dict
                        if no_mismatch_score >= all_score:
                            ir_dict[mismatch_check_seq] = no_mismatch_score
                        else:
                            ir_dict[final_ir_seq] = all_score

            # Now let's return our best-scoring inversion and apply it to the inversion class
            # Also...should we revise the new IR positions? Run before gene looking...
            if not ir_dict:
                i.ir_seq = 'Nothing Found'
            else:
                # There's a possibility multiple IRs have the same score. In this case, I'd choose the pair
                # with the lowest absolute deviation from the original pair?
                best_ir = max(ir_dict, key=lambda key: ir_dict[key])

                i.ir_seq = best_ir
        return


class Inversion:

    def __init__(self, inv_start, inv_end, reads_start, reads_end, ir_seq, parent=None):

        # Core parameters
        self.pos_start = inv_start
        self.pos_end = inv_end
        self.reads_start = reads_start
        self.reads_end = reads_end
        self.ir_seq = ir_seq

        # Linking paramaters
        self.nuccore_accession = parent

        # Analysis parameters
        self.nearby_genes = None

    # Method to associate with a nuccore accession
    def set_nuccore_accession(self, acc):
        self.nuccore_accession = acc
        return

    # Method returns nearby gene descriptors if they are hypothetical.
    def hypotheticals(self):

        if not self.nearby_genes:
            return None

        hypothetical_genes = list()

        for category in self.nearby_genes:
            gene = self.nearby_genes[category]

            if 'hypothetical' in gene.qualifiers['product'][0]:
                try:
                    acc = gene.qualifiers['protein_id'][0]
                    hypothetical_genes.append((acc, category))
                except KeyError:
                    pass

        return hypothetical_genes

    # Processing: if the predicted protein is "Hypothetical", try harder. Consider HMMER tools or IPG database.
    def hmmer_time(self):

        # May be easiest to do an IPG Entrez search and then extract the next one that ain't hypothetical
        if not self.nearby_genes:
            return "No genes! Continuing..."

        for category in self.nearby_genes:
            gene = self.nearby_genes[category]

            # Continue with characterization if hypothetical
            if 'hypothetical' in gene.qualifiers['product'][0]:

                # Make an entrez search to ipg
                try:
                    acc = gene.qualifiers['protein_id'][0]
                    with Entrez.esearch(db='protein', term=acc) as esearch_handle:
                        esearch_data = Entrez.read(esearch_handle)
                        ipg_uid = esearch_data['IdList'][0]
                    with Entrez.efetch(db='protein', id=ipg_uid, retmode='xml', rettype='ipg') as efetch_handle:
                        data = Entrez.read(efetch_handle)

                        ipg_list = data[0]['ProteinList'][0]
                        fxns = defaultdict(int)
                        for ipg in ipg_list:
                            put_fxn = ipg.attributes['name']
                            if 'hypothetical' not in put_fxn and 'Uncharacterized' not in put_fxn and 'Uncharacterised' not in put_fxn:
                                fxns[put_fxn] += 1

                        best_ipg = max(fxns, key=lambda key: fxns[key])
                        # add to inversion instance description, overwrite

                        self.nearby_genes[category].qualifiers['product'] = [best_ipg]

                except:
                    pass

        return

#--------------------------------------------


### MAIN ###
def main():

    my_file_folder = inv_config.cluster_save_path
    my_entrez_folder = inv_config.acc_save_path

    my_samples = process_analysis_folder(my_file_folder, my_entrez_folder, download=False)

    print("{0} samples loaded.".format(len(my_samples)))

    # Now for further processing...
    ir_nt_buffer = inv_config.ir_nt_buffer
    print("Accession processing - aligning genes to inversion sites and detecting novel IR seqs using a {0}nt buffer..."
          .format(ir_nt_buffer), end='')
    hypothetical_inversions = defaultdict(list)
    acc_query = list()
    for b in my_samples:
        for a in b.accessions():
            # Processing 1 - align genes to inversion sites.
            # print("Aligning accession genes to inversion sites...", end='')
            a.detect_nearest_genes()
            # print("Done")

            # Processing 2 - find novel inversion sites
            # print("Detecting novel inverted repeats using a {0}nt buffer...".format(ir_nt_buffer), end='')
            a.find_inverted_repeats(buffer=ir_nt_buffer)
            # print("Done")

            # Processing 3 - query IPG to annotate hypothetical proteins - part 1
            # First build a list of inversions with hypothetical annotations
            for i in a.inversions:
                hyp_genes = i.hypotheticals()
                if hyp_genes:
                    for acc, loc in hyp_genes:
                        if acc not in acc_query:
                            acc_query.append(acc)
                        hypothetical_inversions[acc].append((i, loc))
    print("Done")

    # Now, let's fill in the gaps of the hypotheticals
    if acc_query:
        print("Querying IPG to flesh out hypothetical annotations of {0} genes...".format(len(acc_query)), end='')
        annotate_by_ipg(acc_query, hypothetical_inversions)
        print("Done!")

    # Now, let's start the analysis.

    # Analysis 1: Frequency histogram of inverted repeats.
    print("Analysis 1: Making inverted repeat histogram...", end='')
    ir_histogram_file = 'ir_histogram.csv'
    inverted_repeat_analysis(my_samples, ir_histogram_file)
    print("Complete!")

    # Requires  a description file.
    if os.path.exists('desc.txt'):
        descs = []
        with open('desc.txt', 'r') as desc_in:
            for line in desc_in:
                desc = line[:-1]
                descs.append(desc)

        for desc in descs:
            print("Analysis: Description analysis: {0}...".format(desc), end='')
            desc_file = desc.replace('/', '_') + '.csv'
            descriptor_analysis(my_samples, desc, desc_file)
            print("Done!")

    return


if __name__ == '__main__':

    start = time.time()
    main()
    end = time.time()
    print("Process completed in {0} minutes".format((end-start)/60))



