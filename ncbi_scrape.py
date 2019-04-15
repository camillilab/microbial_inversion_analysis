"""
NCBI scraper. Obtains SRA, Biosample, Assembly, RefSeq info for an organism query to Biosample.

Author: Jacob Bourgeois
Email: jacob.bourgeois@tufts.edu
Affiliation: Camilli Laboratory
License: MIT

Make sure you enter your NCBI credentials into inv_config.py
"""


from Bio import Entrez
import csv
from urllib.error import HTTPError
import time
import re
import requests
import inv_config
import argparse
import os


# grab uses RE to find search terms in ExpXml element in SRA string
def grab(s, q):
    regex = re.compile('(?<={0}=").*?(?=")'.format(q))
    m = regex.search(s)
    if m:
        return m.group()
    else:
        return None


# similar to grab - just looks for an element
def check(s, q):
    regex = re.compile('{0}'.format(q))
    m = regex.search(s)
    if m:
        return True
    else:
        return None


def entrez_query(db, type, retstart, query_key, webenv, max_attempts=3, batch_size=100):

        attempt = 0
        handle = 'FAIL'
        # Try to get the query
        while attempt < max_attempts:
            attempt += 1
            try:
                if type == 'elink':
                    handle = Entrez.elink(dbfrom=db, retstart=retstart, retmax=batch_size, query_key=query_key, webenv=webenv, cmd='acheck')
                    attempt = max_attempts + 1
            except HTTPError as e:
                if 500 <= e.code <= 599:
                    print("Received error from server: {0}".format(e))
                    print("Attempt {0} of {1}".format(attempt, max_attempts))
                    time.sleep(15)
                else:
                    raise
        return handle


# my ncbi info - please don't share haha
Entrez.email = inv_config.email
Entrez.api_key = inv_config.api_key

# handle = Entrez.einfo(db='biosample')
# handle = Entrez.efetch(db='biosample', id='2644363')
# handle = Entrez.elink(dbfrom='biosample', db='nuccore', id=biosample_id)
# handle = Entrez.elink(dbfrom='biosample', id='1766596', cmd='lcheck')
# result = Entrez.read(handle)
# print(result)

batch_size = 1000


# This handle query searches the SRA database for anything for your organism and posts it on the NCBI history server
# Make sure to grab the WebEnv and query key for subsequent esummary grabs

parser = argparse.ArgumentParser()

# REQUIRED ARGUMENTS
parser.add_argument('BioSample_Query', help='Biosample Organism Query')
args = parser.parse_args()
search_term = args.BioSample_Query

# Unfortunately, elink chokes using WebEnv. So I'm retrieving all UIDs at once.
with Entrez.esearch(db='biosample', term=search_term, retmax=100000) as handle:
    result = Entrez.read(handle)
    all_biosample_uids = result['IdList']
    num_results = int(result['Count'])

    print("{0} UIDs retrieved for Biosample query: {1}".format(num_results, search_term))


print("Now searching for SRA and Assembly links in Biosample UIDs...")

biosample_uids_with_sra_assembly = list()  # list of biosample uids with sra AND assembly info

for start in range(0, num_results, batch_size):

    uids = ','.join(all_biosample_uids[start:min(start + batch_size, num_results)])
    print("Requesting {0} of {1}...".format(start, min(start + batch_size, num_results)))

    elink_handle = Entrez.elink(dbfrom='biosample', id=uids, cmd='acheck')
    elink_data = Entrez.read(elink_handle)

    the_goods = elink_data[0]['IdCheckList']['IdLinkSet']  # Container for the iterable data
    for link_data in the_goods:

        biosample_uid = link_data['Id']

        links = link_data['LinkInfo']
        is_sra = False
        is_assembly = False

        for link in links:
            if link['HtmlTag'] == 'SRA':
                is_sra = True
            if link['HtmlTag'] == 'Assembly':
                is_assembly = True

        if is_sra and is_assembly:
            # print("Biosample UID {0} - SRA: {1}, Assembly: {2}".format(biosample_uid, is_sra, is_assembly))
            biosample_uids_with_sra_assembly.append(biosample_uid)

print("{0} biosamples with both SRA and assembly data detected.".format(len(biosample_uids_with_sra_assembly)))

# It does not seem possible to do linked parameters with multiple unique ids in Entrez Biopython. Construct my own URL

print("Now using eLink to obtain SRA UIDs and link together Biosample UIDs with SRA UIDs...")
biosample_uid_sra_dict = dict()
for start in range(0, len(biosample_uids_with_sra_assembly), 100):

    end = min(start + 100, len(biosample_uids_with_sra_assembly))
    print("Constructing eLink URL for SRAs {0} to {1}...".format(start, end))
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?'
    dbfrom = 'biosample'
    db = 'sra'
    current_uids = biosample_uids_with_sra_assembly[start:end]

    id_string = '&id='
    # id_string = ''
    for uid in current_uids:
        id_string = id_string + uid + '&id='
    id_string = id_string[:-4]  # lop off the final &id=

    url = '{0}dbfrom={1}&db={2}{3}&email={4}&api_key={5}'.format(base, dbfrom, db, id_string, Entrez.email, Entrez.api_key)

    params = ({'dbfrom':dbfrom, 'db':db, 'id':id_string, 'email':Entrez.email, 'api_key':Entrez.api_key})

    # post the request
    print("Posting request...")
    r = requests.post(url, params=params)

    print("Writing XML response to file...")
    with open('data.xml', 'w') as f:
        for line in r.text:
            f.write(line)

    # try reading the file in Entrez

    print("Extracting SRA UIDs...")
    with open('data.xml', 'r') as f:
        data = Entrez.read(f)
        for record in data:
            buid = record['IdList'][0]

            try:
                sra_uid = record['LinkSetDb'][0]['Link'][0]['Id']
                biosample_uid_sra_dict[buid] = sra_uid
            except IndexError:
                pass

            # NOTE: It may be possible for a single BioSample to be associated with more than one SRA UID.
            # For now, let's take the first one, for simplicity.

biosample_uid_sra_dict.pop('0', None)  # Removes the zero element obtained from XML Parser
all_sras = list(biosample_uid_sra_dict.values())  # get all SRA uids
print("Done. Retrieved {0} SRA UIDs".format(len(all_sras)))

print("Now retrieving SRA data using esummary to find datasets with paired end reads...")
paired_sra_uids = list()
sra_uid_acc_dict = dict()
sra_paired_dict = dict()
for start in range(0, len(all_sras), batch_size):

    end = min(len(all_sras), start + batch_size)
    print("Requesting eSummary {0} of {1}...".format(start, end))
    sra_uids = ','.join(all_sras[start:end])

    esummary_handle = Entrez.esummary(db='sra', id=sra_uids)
    esummary_data = Entrez.read(esummary_handle)

    for record in esummary_data:

        sra_id = record['Id']

        expxml = record['ExpXml']
        run_data = record['Runs']
        paired_info = check(expxml, 'PAIRED')
        exp_acc = grab(expxml, 'Experiment acc')
        run_acc = grab(run_data, 'Run acc')
        is_paired = 0
        if paired_info:
            is_paired = 1
            paired_sra_uids.append(sra_id)
        sra_uid_acc_dict[sra_id] = run_acc
        sra_paired_dict[sra_id] = is_paired

print("Of the {0} SRA UIDs, {1} contain paired data reads".format(len(all_sras), len(paired_sra_uids)))

# Create a subdict containing biosample uids with paired info
biosample_uids_paired = {k: v for k, v in biosample_uid_sra_dict.items() if v in paired_sra_uids}

print("Building Biosample:Assembly Dict...")

biosample_uid_assembly_dict = dict()
for start in range(0, len(biosample_uids_with_sra_assembly), 100):

    end = min(start + 100, len(biosample_uids_with_sra_assembly))
    print("Constructing eLink URL for Assembly UIDs {0} to {1}...".format(start, end))
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?'
    dbfrom = 'biosample'
    db = 'assembly'
    current_uids = biosample_uids_with_sra_assembly[start:end]

    id_string = '&id='
    # id_string = ''
    for uid in current_uids:
        id_string = id_string + uid + '&id='
    id_string = id_string[:-4]  # lop off the final &id=

    url = '{0}dbfrom={1}&db={2}{3}&email={4}&api_key={5}'.format(base, dbfrom, db, id_string, Entrez.email, Entrez.api_key)

    params = ({'dbfrom':dbfrom, 'db':db, 'id':id_string, 'email':Entrez.email, 'api_key':Entrez.api_key})

    # post the request
    print("Posting request...")
    r = requests.post(url, params=params)

    print("Writing XML response to file...")
    with open('data.xml', 'w') as f:
        for line in r.text:
            f.write(line)

    # try reading the file in Entrez

    print("Extracting Assembly UIDs...")
    with open('data.xml', 'r') as f:
        data = Entrez.read(f)
        for record in data:
            buid = record['IdList'][0]

            try:
                assembly_uid = record['LinkSetDb'][0]['Link'][0]['Id']
                biosample_uid_assembly_dict[buid] = assembly_uid
            except IndexError:
                pass

            # NOTE: It may be possible for a single BioSample to be associated with more than one SRA UID.
            # For now, let's take the first one, for simplicity.

biosample_uid_assembly_dict.pop('0', None)  # Removes the zero element obtained from XML Parser
all_assembly = list(biosample_uid_assembly_dict.values())  # get all assembly uids
print("Done. Retrieved {0} assembly UIDs".format(len(all_assembly)))


print("Building Assembly UID:Accession Dict...")
assembly_uid_acc_dict = dict()
assembly_status_dict = dict()  # dict contains the assembly status
for start in range(0, len(all_assembly), 100):

    end = min(len(all_assembly), start + 100)
    print("Requesting eSummary {0} of {1}...".format(start, end))
    assembly_uids = ','.join(all_assembly[start:end])

    esummary_handle = Entrez.esummary(db='assembly', id=assembly_uids)
    esummary_data = Entrez.read(esummary_handle)

    i = 0  # obnoxiously, metadata does not contain original UID. Thanks NCBI
    for record in esummary_data['DocumentSummarySet']['DocumentSummary']:

        current_uid = all_assembly[start:end][i]
        acc = record['AssemblyAccession']
        status = record['AssemblyStatus']

        i += 1

        assembly_uid_acc_dict[current_uid] = acc
        assembly_status_dict[current_uid] = status

print("Now retrieving nuccore UIDs associated with Assembly UIDs...")
assembly_nuccore_dict = dict()  # dict contains assembly UIDs keyed to a list of nuccore uids
all_nuccore_uids = list()


for start in range(0, len(all_assembly), 100):

    end = min(start + 100, len(all_assembly))
    print("Constructing eLink URL for Assembly UIDs {0} to {1}...".format(start, end))
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?'
    dbfrom = 'assembly'
    db = 'nuccore'
    current_uids = all_assembly[start:end]

    id_string = '&id='
    # id_string = ''
    for uid in current_uids:
        id_string = id_string + uid + '&id='
    id_string = id_string[:-4]  # lop off the final &id=

    url = '{0}dbfrom={1}&db={2}{3}&email={4}&api_key={5}'.format(base, dbfrom, db, id_string, Entrez.email, Entrez.api_key)

    params = ({'dbfrom':dbfrom, 'db':db, 'id':id_string, 'email':Entrez.email, 'api_key':Entrez.api_key})

    # post the request
    print("Posting request...")
    r = requests.post(url, params=params)

    print("Writing XML response to file...")
    with open('data.xml', 'w') as f:
        for line in r.text:
            f.write(line)

    # try reading the file in Entrez

    print("Extracting Nuccore UIDs...")
    with open('data.xml', 'r') as f:
        data = Entrez.read(f)
        for record in data:
            assembly_uid = record['IdList'][0]
            nuccore_uids = list()
            try:
                links = record['LinkSetDb']
                for link in links:
                    if link['LinkName'] == 'assembly_nuccore_refseq':  # The RefSeq link should contain all the juicy bits with annotations!
                    # if link['LinkName'] == 'assembly_nuccore_insdc':
                        for nuids in link['Link']:
                            nuid = nuids['Id']
                            all_nuccore_uids.append(nuid)
                            nuccore_uids.append(nuid)
                            assembly_nuccore_dict[assembly_uid] = nuccore_uids
            except IndexError:
                print("Something went wrong! {0}".format(assembly_uid))

            # Sometimes, there is no RefSeq, thus no annotation. Discard at end
            if len(nuccore_uids) == 0 and assembly_uid != '0':
                print("Warning! No RefSeq UIDs discovered for {0}...".format(assembly_uid))
                assembly_nuccore_dict.pop(assembly_uid, None)


assembly_nuccore_dict.pop('0', None)  # Removes the zero element obtained from XML Parser
os.remove('data.xml')
print("Done. Retrieved {0} nuccore uids.".format(len(all_nuccore_uids)))

print("Building BioSample UID:Accession Dict...")
biosample_uid_acc_dict = dict()  # Accession info
biosample_uid_desc_dict = dict()  # Project information
biosample_uid_org_dict = dict()  # Organism info
for start in range(0, len(biosample_uids_with_sra_assembly), 100):

    end = min(len(biosample_uids_with_sra_assembly), start + 100)
    print("Requesting eSummary {0} of {1}...".format(start, end))
    biosample_uids = ','.join(biosample_uids_with_sra_assembly[start:end])

    esummary_handle = Entrez.esummary(db='biosample', id=biosample_uids)
    esummary_data = Entrez.read(esummary_handle)

    for record in esummary_data['DocumentSummarySet']['DocumentSummary']:

        b_id = record.attributes['uid']
        acc = record['Accession']
        org = record['Organism']
        strain = record['Infraspecies']

        biosample_uid_acc_dict[b_id] = acc
        biosample_uid_org_dict[b_id] = org
        biosample_uid_desc_dict[b_id] = strain

# Now get the nuccore uid::acc dict
print("Building Nuccore UID:Accession Dict...")
nuccore_uid_acc_dict = dict()
for start in range(0, len(all_nuccore_uids), 1000):

    end = min(len(all_nuccore_uids), start + 1000)
    print("Requesting eSummary {0} of {1}...".format(start, end))
    nuccore_uids = ','.join(all_nuccore_uids[start:end])

    #DEBUG!
    # nuccore_uids = '531475993'

    esummary_handle = Entrez.esummary(db='nuccore', id=nuccore_uids)
    esummary_data = Entrez.read(esummary_handle)

    for record in esummary_data:

        nu_id = record['Id']
        acc = record['AccessionVersion']

        nuccore_uid_acc_dict[nu_id] = acc

    time.sleep(3)  # Needed because sometimes a bad response is given

print("Done.")

print("Now outputting final table...")
# Alright, so we want biosample accession, biosample description, biosample organism, assembly acccession, sra accession, and nuccore accessions

fieldnames = ['Biosample', 'Organism', 'Strain', 'SRA', 'Assembly', 'Assembly Status', 'RefSeq Accessions']
with open('ncbi_data.txt', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(fieldnames)

    for biosample_uid in biosample_uids_paired:

        try:

            accession = biosample_uid_acc_dict[biosample_uid]
            desc = biosample_uid_desc_dict[biosample_uid]
            org = biosample_uid_org_dict[biosample_uid]
            sra_uid = biosample_uid_sra_dict[biosample_uid]
            sra_acc = sra_uid_acc_dict[sra_uid]
            assembly_uid = biosample_uid_assembly_dict[biosample_uid]
            assembly_acc = assembly_uid_acc_dict[assembly_uid]
            status = assembly_status_dict[assembly_uid]
            nuccore_uids = assembly_nuccore_dict[assembly_uid]
            nuccore_accs = []
            for uid in nuccore_uids:
                nuccore_accs.append(nuccore_uid_acc_dict[uid])
            nuccore_accs = ','.join(nuccore_accs)

            if nuccore_uids != '':

                writer.writerow((accession, org, desc, sra_acc, assembly_acc, status, nuccore_accs))

        except KeyError:
            pass

