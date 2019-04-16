# microbial_inversion_analysis
Analysis of Inversion Switches across Microbial Genomes

This repository is a collection of scripts designed to work with analyze_clusters to extract BioSamples related to an organism of interest that have corresponding assembly and paired SRA data.

Then, a series of scripts automates pulling of SRA FASTQ files, alignment to the RefSeq assembly, and subsequent exrtraction of single-orientation reads (SOR) and 5'-clipped reads (sCLIP) to use in the detection of inverted repeats from sequencing data as described here:

Sekulovic O, Mathias Garrett E, Bourgeois J, Tamayo R, Shen A, et al. (2018) Genome-wide detection of conservative site-specific recombination in bacteria. PLOS Genetics 14(4): e1007332. https://doi.org/10.1371/journal.pgen.1007332

---------------------------------

# Requirements:

This suite requires:
- Python 3.7+
- ASCP file transfer protocol command line tools (https://developer.asperasoft.com/desktop-advance/command-line-client)
- BioPython
- BLAST command line tools, particularly blastn and blastp
- Requests module
- Seaborn
- ReportLab
- Numpy
- MatPlotLib

---------------------------------

# HOW TO:

0. Enter your NCBI email and api-key into <b>inv_config.py</b>. This is a courtesy to NCBI when making several queries to their servers, increases the amount of info per second, and prevents you from being ip-locked. The values are for email and api-key respectively.

1. Obtain the Biosample, Assembly, SRA, and Refseq information from NCBI using <b>ncbi_scrape.py</b>. Ensure the query is in double-quotes.

`python ncbi_scrape.py "BIOSAMPLE_QUERY"`

So, for example, for out work with C.dif:

`python ncbi_scrape.py "clostridioides difficile[organism]"`

Sometimes, you may receive a URL error, particularly during the RefSeq interrogation. In this case, simply re-running the program usually works.

2. Run the main script, <b>inversion_detection.py</b>, which automates SRA extraction of FASTQ, mapping of reads using bowtie2, and SAM file manipulation to get SOR/sCLIP reads. Then, feeds into the inversion detection script (the requisite files are included in this distribution).

`python inversion_detection.py`

This while take a WHILE. Consider increasing the thread count, splitting the computational load, running in the background, or running in several batches. For example,

`python inversion_detection.py | tee log.txt &`

3. Run the analysis script, <b>analyze_inversions.py</b>, which outputs a number of analysis files from the Cluster Data.

You may also run analyses over certain descriptors, such as 'GntR family transcriptional regulator'. Add single-line descriptions in a text file in the repository directory called <b>desc.txt</b>. The script will then pull out the biosamples that correspond to inversions involving those descriptors.

The output histogram, <b>ir_histogram.csv</b>, gives a breakdown of the detected inverted repeats. The analysis script is still in it's nascent stages, but you may use the classes available to write your own analyses.

We ask that you please cite out paper, CITATION, if you use our scripts in our publication.

Please direct any questions to jacob.bourgeois@tufts.edu or Ognjen.Sekulovic@tufts.edu .

Thank you!

----------------------

# Final Notes

One major limitation of the current implementation of this strategy is the relative paucity of BioSample records that have both NCBI Assembly and SRA data. For example, there are nearly 15,000 SRA records for "clostridioides difficile", but only ~530 (~3%) has Assembly data. We could do without the RefSeq data in Assembly, but we lose any nearby gene loci information that can be used to draw regulatory hypotheses.

We are currently working on adding a second pipeline that uses <i>de novo</i> assembly and command-line prokaryotic annotation to overcome this limitation. However, this will vastly increase runtime and using default parameters for assemblers and annotators may be inappropriate. Still, this could help include the other 97% of sequencing data available on the SRA database in the detection and validation of invertible switches. If you have any suggestions, please let me know!
