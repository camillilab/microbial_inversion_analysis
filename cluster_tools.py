#! usr/bin/python

"""

Toolkit for detection and analysis of genomic inversions using same-orientation reads

Author: Jacob Bourgeois
Email: jacob.bourgeois@tufts.edu
Organization: Tufts Sackler School of Biomedical Sciences; Camilli Lab

License: 3-clause BSD

"""

# IMPORTS

try:
    import config
    from write_tools import *
    import seaborn as sns
    from collections import defaultdict
    import matplotlib.pyplot as plt
    import numpy as np
    import csv
    import sys
    import os
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from Bio import Entrez
    from Bio import SeqIO
    from Bio.Seq import Seq, MutableSeq
    from Bio import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Graphics import GenomeDiagram
except ImportError as e:
    sys.exit("Import warning! Make sure dependencies are installed. {}".format(e))


# class BUG represents an organism and is the master organism class in our cluster detection algorithm
class Bug:

    def __init__(self, acc_num):

        # Characteristics of the organism:
        self.accession_num = acc_num    # Bug accession number needed to pull correct reference genome from NCBI
        self.genes = list()             # list of CDS detected in entrez file
        self.name = 'NONE'              # Name of organism
        self.sequence = Seq('')         # sequence of organism, needed for inverted repeat analysis

        # CLUSTERING VALUES

        # SOR data
        self.SOR_pos_freq_dict = defaultdict(int)       # contains position-frequency data of the SOR file
        self.SOR_pos_array = np.array([])               # array of unique positions in SOR
        self.SOR_ignored_positions = []                 # list of positions to ignore when loading SOR data
        self.SOR_read_sum = 0                           # sum of all read counts in SOR data
        self.SOR_pos_min = -1                           # minimum position in SOR data
        self.SOR_pos_max = 0                            # maximum position in SOR data
        self.SOR_read_cutoff = 0                        # density value of SOR histogram cutoff to be a cluster
        self.SOR_histogram = []                         # histogram of SOR data generated during thresholding

        # Initial screen clustering parameters
        self.SOR_bin_size = config.nbin_size            # size of bin in nt to represent the initial SOR histogram
        self.SOR_final_bin_size = 0                     # what bin size we got in the end

        # sCLIP data for clustering validation
        self.sCLIP_pos_freq_dict = defaultdict(int)     # contains pos-freq data of the sCLIP file
        self.sCLIP_pos_array = np.array([])             # array of unique positions in sCLIP
        self.sCLIP_ignored_positions = []               # list of positions to ignore when loading sCLIP data
        self.sCLIP_read_sum = 0                         # sum of all read counts in sCLIP data
        self.sCLIP_pos_min = -1                         # minimum position in sCLIP data
        self.sCLIP_pos_max = 0                          # maximum position in sCLIP data

        # Data from clustering analysis
        self.clusters = []                              # list of all clusters filtered from SOR histogram screen
        self.signals = []                               # list of signals identified from cluster analysis
        self.inversions = []                            # list of signals that have at least two identifiable positions
        self.spikes = []                                # list of signals that only have one significant position
        self.loci = []                                  # loci nearby inversions. should line up with inversion index
        self.is_intergenic = []                         # flags if intergenic. Index lines with inversion
        self.rcscores = []                              # list of reverse complement percentages

    # method load_sor loads SOR file from default directory and creates SOR_pos_freq_dict and SOR_pos_array
    def load_SOR(self, sor_file):

        print("Loading {0} as SOR data...".format(sor_file), end='')

        with open(sor_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:

                    # ignore TLEN values that are less than zero
                    if int(row['TLEN']) != 0:

                        pos = int(row['POS'])

                        # ignore positions that are in ignored positions
                        if pos not in self.SOR_ignored_positions:
                            self.SOR_pos_freq_dict[pos] += 1
                            self.SOR_read_sum += 1

                except csv.Error as e:
                    print("CSV read error occurred! Please ensure headers in the SOR file include TLEN and POS.")
                    sys.exit('file {}, line {}: {}'.format(sor_file, reader.line_num, e))

        # get rid of any positions less then the thresholding value
        killpos = list()
        for pos in self.SOR_pos_freq_dict:
            if self.SOR_pos_freq_dict[pos] < config.read_count_threshold:
                killpos.append(pos)
        for pos in killpos:
            del self.SOR_pos_freq_dict[pos]

        # create the pos array as a sorted array
        self.SOR_pos_array = np.sort(np.array(list(self.SOR_pos_freq_dict)))

        # assign minimum and maximum values
        self.SOR_pos_min = self.SOR_pos_array.min()
        self.SOR_pos_max = self.SOR_pos_array.max()

        print("{0} locations loaded with {1} reads.".format(len(self.SOR_pos_array), self.SOR_read_sum))

        return

    # method load_sCLIP loads SOR file from default directory and creates sCLIP_pos_freq_dict and sCLIP_pos_array
    def load_sCLIP(self, sCLIP_file):

        print("Loading {0} as sCLIP data...".format(sCLIP_file), end='')

        with open(sCLIP_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:

                    # ignore TLEN values that are less than zero
                    if int(row['TLEN']) != 0:

                        pos = int(row['POS'])

                        # ignore positions that are in ignored positions
                        if pos not in self.sCLIP_ignored_positions:
                            self.sCLIP_pos_freq_dict[pos] += 1
                            self.sCLIP_read_sum += 1

                except csv.Error as e:
                    print("CSV read error occurred! Please ensure headers in the sCLIP file include TLEN and POS.")
                    sys.exit('file {}, line {}: {}'.format(sCLIP_file, reader.line_num, e))

        # create the pos array as a sorted array
        self.sCLIP_pos_array = np.sort(np.array(list(self.sCLIP_pos_freq_dict)))

        # assign minimum and maximum values
        self.sCLIP_pos_min = self.sCLIP_pos_array.min()
        self.sCLIP_pos_max = self.sCLIP_pos_array.max()

        print("{0} locations loaded with {1} reads.".format(len(self.sCLIP_pos_array), self.sCLIP_read_sum))

        return

    # method make_interactive_graphical_threshold creates a histogram of the SOR data for the user to manually select
    # a density cutoff to screen potential inversion clusters using matplotlib
    # returns clusters
    # added an automatic function for laziness
    def make_SOR_interactive_graphical_threshold(self, automatic=False):

        print ("Creating thresholder for SOR data with {0}nt bins...".format(config.nbin_size))

        # class HLineBuilder allows us to define a density cutoff in the initial screen.
        class HLineBuilder:
            def __init__(self, line, x_bin):
                self.line = line
                self.xs = (0, x_bin)
                self.ys = line.get_ydata()
                self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

            def __call__(self, event):
                if event.inaxes != self.line.axes: return
                y1, y2 = event.ydata, event.ydata
                self.line.set_data(self.xs, (y1, y2))
                self.y_final = y1
                self.line.figure.canvas.draw()

        # define the bounds of the graph using data range and nbins
        seq_size = self.SOR_pos_max - self.SOR_pos_min
        nbins = int(seq_size / self.SOR_bin_size)

        # make an array of data containing elements at their frequency, ex. 1 2 2 3 3 3 ....
        data = list()
        for pos in self.SOR_pos_array:
            for i in range(0, self.SOR_pos_freq_dict[pos]):
                data.append(pos)
        data = np.array(data)

        # now make a frequency histogram using numpy
        h_densities, den_bin_edges = np.histogram(data, bins=nbins, density=True)
        self.SOR_final_bin_size = den_bin_edges[1] - den_bin_edges[0]
        self.SOR_histogram = h_densities

        if not automatic:
            # Plot the density histogram, providing visual representation of read densities
            fig, ax1 = plt.subplots()
            ax1.plot(h_densities)  # plot density histogram along axis.

            # Call a linebuilder class to allow the user to draw a cutoff visually.
            line, = ax1.plot([0], [0])  # empty line
            r = HLineBuilder(line, nbins)

            # Define plot parameters
            plt.title("Click to set read density cutoff")
            plt.xlabel("Bin")
            plt.ylabel("Bin read density")

            plt.show()

            # Our cutoff is equal to the y value of the last line drawn
            self.SOR_read_cutoff = r.y_final
            plt.close()

        # if automatic - use config threshold to define cutoff cluster read density
        if automatic:
            print("Using automatic thresholding at percentile {0}...".format(config.automatic_thresholding_perc))
            self.SOR_read_cutoff = np.percentile(h_densities, config.automatic_thresholding_perc)

        # add the left-sided bin edges to a list if they pass; these represent the left side of a potential cluster
        h_bin_left_pos_list = list()
        for i in range(0, len(h_densities)):
            if h_densities[i] >= self.SOR_read_cutoff:
                h_bin_left_pos_list.append(float(den_bin_edges[i]))

        for left_bin_edge in h_bin_left_pos_list:
            right_bin_edge = left_bin_edge + self.SOR_final_bin_size

            # Create analysis Cluster instances based on sCLIP data
            cluster = self.sCLIP_subset(left_bin_edge, right_bin_edge)
            self.clusters.append(cluster)

        print("Read Density threshold: ", self.SOR_read_cutoff)
        print("Signals detected:", len(self.clusters))

        return self.clusters

    # recreates SOR thresholding graph for saving
    def save_SOR_thresholding(self, save_path, figsize=(8, 6)):
        fig, ax1 = plt.subplots(figsize=figsize)
        ax1.plot(self.SOR_histogram)  # plot density histogram along axis.
        plt.title(self.accession_num + " SOR Density Histogram; Bins={0}".format(self.SOR_bin_size))
        plt.xlabel("Bin Number")
        plt.ylabel("Bin read density")
        plt.axhline(self.SOR_read_cutoff, color='r')
        plt.savefig(save_path)
        plt.close()

    # method returns a Cluster class based on a subset of organism SOR data
    def SOR_subset(self, start, end):
        sub_dict = dict()
        for element in self.SOR_pos_array:
            if (element >= start) and (element <= end):
                sub_dict[element] = self.SOR_pos_freq_dict[element]
        return Cluster(sub_dict,
                       cbinsize=config.cbin_size,
                       cperc=config.cbin_cutoff,
                       clustersepmin=config.cluster_min_sep,
                       clustersepmax=config.cluster_max_sep,
                       ntsepmin=config.ntpair_min_sep,
                       ntsepmax=config.ntpair_max_sep
                       )

    # Returns a Cluster instance of sCLIP data between start and end
    def sCLIP_subset(self, start, end):
        sub_dict = dict()
        for element in self.sCLIP_pos_array:
            if (element >= start) and (element <= end):
                sub_dict[element] = self.sCLIP_pos_freq_dict[element]
        return Cluster(sub_dict,
                       cbinsize=config.cbin_size,
                       cperc=config.cbin_cutoff,
                       clustersepmin=config.cluster_min_sep,
                       clustersepmax=config.cluster_max_sep,
                       ntsepmin=config.ntpair_min_sep,
                       ntsepmax=config.ntpair_max_sep)

    # method uses entrez gb data to get gene and CDS information into Bug instance using BioPython Entrez and SeqFeature
    def get_entrez_data(self, email, filename):

        Entrez.email = email  # always give the NCBI an email in case you abuse this :)

        # if we don't have the file downloaded, let's do so
        if not os.path.exists(filename):
            print("Downloading Entrez data for accession number {}....".format(self.accession_num), end='')
            net_handle = Entrez.efetch(db='nucleotide', id=self.accession_num, rettype='gb', retmode='text')
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print("Saved.")

        print("Parsing GenBank data...", end='')
        # load the file onto an entrez parser
        handle = open(filename)
        record = SeqIO.read(handle, format='gb')

        # load the bare sequence onto bug instance
        self.sequence = record.seq

        # now, let's go through the genbank record and load all CDS into our gene list in the Bug class
        elements = record.features
        for element in elements:
            # check to see if we have a CDS
            if element.type == "CDS":
                # if so, load the sequence features onto bug instance
                self.genes.append(element)

        self.name = record.annotations['organism']

        print("{0} CDS detected in organism {1}.".format(len(self.genes), self.name))

        # if we got a bad genbank file, then we may not have any CDSs. Abort!
        if len(self.genes) == 0:
            print("No CDS Detected! Accession number may be wrong or pointing to an incomplete genbank reference."
                  "Please provide a different accession number for genetic alignment.")
            quit()

        return

    # NEW! 02_04_2019
    # This method is now being added to support mapping using multiple contigs and the like.

    # returns a position-frequency dictionary containing sCLIP data between given positions
    def _sCLIP_pfd_slice(self, start, end):

        a = dict()
        for element in self.sCLIP_pos_array:
            if (element >= start) and (element <= end):
                a[element] = self.sCLIP_pos_freq_dict[element]
        return a

    # method attempts to locate organism genes around a inversion pairs given a nt tolerance and gene proximity max
    # returns a list of biopython seqfeatures (ie. genes) and is_intergenic
    def align_inversion_to_genes(self, inv_pair, ntol=10000, max_genes=1000):

        print("Aligning inversion pair {0} to genes from {1}...".format(inv_pair, self.name), end='')

        # inversion pair region in which to search for genes +/- some nts
        inv_min = inv_pair[0] - ntol
        inv_max = inv_pair[1] + ntol

        # represent middle of inversion region as average of min and max
        inv_avg = (inv_pair[0] + inv_pair[1]) / 2

        # total number of genes in organism, so we don't get an index error by going to far
        total_genes = len(self.genes)

        # hit_scores list tells us the difference of distance of the middle of the gene to the middle of the inv
        hit_scores = list()

        # vars
        gene_start = -1  # nt location of start of gene
        i = 0
        loci = list()  # genes that make the cut

        # while the beginning of the gene location does not exceed the inversion max position
        while (gene_start <= inv_max) and (i < total_genes):

            gene_start = self.genes[i].location.start
            gene_end = self.genes[i].location.end
            gene_avg = (gene_start + gene_end) / 2

            # does the end of the gene peek into the start of the inversion range?
            if (gene_end >= inv_min) and (gene_start <= inv_min):
                loci.append(self.genes[i])
                hit_scores.append(abs(inv_avg - gene_avg))

            # does the gene lie squarely within the inversion region?
            if (gene_start >= inv_min) and (gene_end <= inv_max):
                loci.append(self.genes[i])
                hit_scores.append(abs(inv_avg - gene_avg))

            # does the start of the gene lie within the inversion region?
            if (gene_start <= inv_max) and (gene_end >= inv_max):
                loci.append(self.genes[i])
                hit_scores.append(abs(inv_avg - gene_avg))

            i += 1

            # if the number of genes we got exceeded our threshold, trim off the edges
            if len(loci) > max_genes:
                excess = len(loci) - max_genes

                for x in range(0, excess):
                    sorted_scores = sorted(hit_scores, reverse=True)

                    # r is our element to remove based on the highest distance score
                    r = hit_scores.index(sorted_scores[0])
                    loci.pop(r)
                    hit_scores.pop(r)

        print("Done.")

        # Append loci to bug. In theory, should line up with inversions by index.
        self.loci.append(loci)

        # check to see if our inversion site is intergenic
        is_intergenic = self._is_intergenic(inv_pair, loci)

        return loci, is_intergenic

    # draws a gene diagram showing the cluster and nearby loci using GenomeDiagram
    def draw_gene_diagram(self, inv_pair, genes, save_path):

        # print("Drawing gene diagram for", inv_pair)

        # define the tick interval based on the start and end of the genes (should already be ordered)
        track_start = genes[0].location.start
        track_end = genes[len(genes) - 1].location.end
        s_tick_int = int((track_end - track_start) / 5)

        # create an empty genome diagram
        gdd = GenomeDiagram.Diagram(self.accession_num)
        gdt_features = gdd.new_track(1, greytrack=True, scale_smalltick_interval=s_tick_int,
                                     scale_smalltick_labels=True,
                                     scale_smallticks=0.1, scale_fontangle=0, scale_fontsize=4, name=self.accession_num)
        gds_features = gdt_features.new_set()

        # for each loci, annotate the diagram
        for orf in genes:

            # describe the orf
            loctag = orf.qualifiers['locus_tag'][0]
            product = orf.qualifiers['product'][0]

            # define orientation based on strand
            if orf.strand == 1:
                angle = 15
                pos = 'left'
            if orf.strand == -1:
                angle = -195
                pos = 'right'

            # draw the orf
            gds_features.add_feature(orf, name=loctag + ": " + product, label=True, sigil="BIGARROW",
                                     label_size=4, arrowhead_length=0.2, label_angle=angle,
                                     label_position=pos, arrowshaft_height=0.3)


        # for the cluster, annotate inversion positions

        feature = SeqFeature(FeatureLocation(int(inv_pair[0]), int(inv_pair[0]) + 1), strand=0)
        gds_features.add_feature(feature, name='   START',
                                 label=True, color="purple", label_position="left",
                                 label_angle=45, sigil='BOX', label_color='purple', label_size=6)

        feature = SeqFeature(FeatureLocation(int(inv_pair[1]), int(inv_pair[1]) + 1), strand=0)
        gds_features.add_feature(feature, name='   END',
                                 label=True, color="purple", label_position="left",
                                 label_angle=45, sigil='BOX', label_color='purple', label_size=6)

        # draw and save the graph
        gdd.draw(format='linear', pagesize=(16 * cm, 10 * cm), fragments=1,
                 start=track_start - 500, end=track_end + 500)
        gdd.write(save_path, "pdf")

        return

    # sees if the inversions are intergenic or not
    def _is_intergenic(self, inv_pair, loci):

        for orf in loci:
            start = orf.location.start
            end = orf.location.end

            # does the end of the gene peek into the start of the inversion range?
            if (end >= inv_pair[0]) and (start <= inv_pair[0]):
                return 'N'

            # does the gene lie squarely within the inversion region?
            if (start >= inv_pair[0]) and (end <= inv_pair[1]):
                return 'N'

            # does the start of the gene lie within the inversion region?
            if (start <= inv_pair[1]) and (end >= inv_pair[1]):
                return 'N'

        # if there's no overlap, return yes
        return 'Y'

    # sees if there's any inverted repeat complementation
    # one thing to do is to add another loop to look upstream for more matches
    def align_inverted_repeats(self, inv_pair, seed, start, end, tries, tol):

        # grab the sequence upstream and downstream of the inversion switch dictated by start and end
        up_seq = self.sequence[inv_pair[0] + start:]
        down_seq = self.sequence[inv_pair[1] + start:]

        for i in range(0, end-seed-start):

            # define the seed sequence
            seed_seq = up_seq[i:i + seed]

            for j in range(0, tries):

                # define the downstream segment
                d_seq = down_seq[j:j + seed]

                # get the reverse complement
                rc = d_seq.reverse_complement()

                # if it matches, start building the sequence here
                if seed_seq == rc:

                    match_seq = seed_seq
                    fails = 0
                    k = -1
                    fix_indicies = []

                    # to avoid index errors counting down, define a new down_seq where it looks backwards
                    r_seq = self.sequence[inv_pair[0]:inv_pair[1] + j + start][::-1]

                    # keep looking while we haven't exceeded failure
                    while fails <= tol:

                        # add one nt to the seed sequence
                        k += 1
                        seed_seq = seed_seq + up_seq[i + seed + k]

                        # get the downstream seq and its rc
                        # add to front, because we are looking backwards on the downstream segment
                        d_seq = r_seq[k] + d_seq
                        rc = d_seq.reverse_complement()

                        # check for matching
                        if seed_seq != rc:
                            fails += 1
                            fix_indicies.append(k + seed)
                        else:
                            fails = 0

                        # allow the d_seq to be the seed_seq to allow for continued repeat search despite mismatches
                        d_seq = seed_seq.reverse_complement()

                        # set the match_seq to the current seed_seq
                        match_seq = seed_seq

                    # turn those mismatches into n's
                    match = match_seq.tomutable()
                    for index in fix_indicies:
                        match[index] = 'n'

                    # lop off fails off the end; these didn't match
                    return match[:-fails]

        return 'Nothing found'

    # make analysis file
    def make_analysis_file(self, save_path):

        print("Creating analysis file for {0}...".format(self.accession_num), end='')

        # Organism info
        append_to_csv(['Organism:', self.name], save_path)
        append_to_csv(['Accession number:', self.accession_num], save_path)

        # Overall info
        labels = ['Number of signals detected', 'Number of inversion pairs detected']
        data = [len(self.clusters), len(self.inversions)]
        for i in range(0, len(labels)):
            d = (labels[i], data[i])
            append_to_csv(d, save_path)
        append_to_csv([''], save_path)

        # Cluster info
        header = ['Signal Start', 'Start Reads', 'Signal End', 'End Reads', 'True Pair?', 'Inversion Length',
                  'Combined Read Count', 'Percent Read to Cluster', 'Percent Read to All SOR Reads']
        append_to_csv(header, save_path)

        for cluster in self.clusters:

            if cluster.is_single_signal == 0:
                sig = 'Y'
            else:
                sig = 'N'

            start = cluster.best_nt_pair[0][0]
            sreads = cluster.best_nt_pair[0][1]
            end = cluster.best_nt_pair[1][0]
            ereads = cluster.best_nt_pair[1][1]
            length = end - start
            creads = sreads + ereads
            preads = 100 * (creads / cluster.reads)
            ptreads = 100 * (creads / self.SOR_read_sum)

            data = [start, sreads, end, ereads, sig, length, creads, preads, ptreads]
            append_to_csv(data, save_path)
        append_to_csv([''], save_path)

        # Maybe you can throw in run parameters if you'd like...just reference config
        print("Done.")

        return

    # make inversion file containing the goods
    def make_inversion_file(self, save_path):

        print("Making inversions file...")

        # Organism info
        append_to_csv(['Organism:', self.name], save_path)
        append_to_csv(['Accession number:', self.accession_num], save_path)

        # Overall info
        labels = ['Number of signals detected', 'Number of inversion pairs detected']
        data = [len(self.clusters), len(self.inversions)]
        for i in range(0, len(labels)):
            d = (labels[i], data[i])
            append_to_csv(d, save_path)
        append_to_csv([''], save_path)

        # Inversion pair info
        headers = ['Cluster Number', 'Inversion start', 'Start Reads', 'Inversion end', 'End reads', 'Inversion length',
                   'Detected inverted repeats around start',
                   'Intergenic?', 'Nearby locus_tags']
        append_to_csv(headers, save_path)

        # data
        i = 0  # for index alignment
        for inversion in self.inversions:

            inv_pair = inversion.best_nt_pair
            start = inv_pair[0][0]
            sreads = inv_pair[0][1]
            end = inv_pair[1][0]
            ereads = inv_pair[1][1]
            length = end - start
            comp = self.rcscores[i]
            is_intergenic = self.is_intergenic[i]
            tags = []
            orfs = self.loci[i]
            for orf in orfs:
                tags.append(orf.qualifiers['locus_tag'][0])
            i += 1

            data = [i, start, sreads, end, ereads, length, comp, is_intergenic, tags]

            append_to_csv(data, save_path)

        return


# Cluster instance takes a position-frequency dictionary and can perform inversion detection algorithms
class Cluster:

    def __init__(self, pos_freq_dict, cbinsize=40, cperc=98,
                 clustersepmin=0, clustersepmax=10000, ntsepmin=0, ntsepmax=10000):

        # Cluster data characteristics
        self.pos_freq_dict = pos_freq_dict                              # position:frequency dictionary of the cluster
        self.pos_array = np.array(list(pos_freq_dict))                  # unique position array of the cluster
        self.cluster_start = self.pos_array.min()                       # cluster start
        self.cluster_end = self.pos_array.max()                         # cluster end
        self.reads = self._count()                                      # total number of reads in cluster
        self.freq_min = np.array(list(pos_freq_dict.values())).min()    # minimum read frequency
        self.freq_max = np.array(list(pos_freq_dict.values())).max()    # maximum read frequency
        self.pos_freq_max = list(pos_freq_dict.keys())[list(pos_freq_dict.values()).index(self.freq_max)]

        # Clustering parameters
        self.bin_size = cbinsize                    # how many nucleotides each bin should span
        self.c_sep_min = clustersepmin              # limit to how close cluster pairs can be
        self.c_sep_max = clustersepmax              # limit to how far cluster pairs can be
        self.n_sep_min = ntsepmin                   # limit to how close nucleotide pairs can be
        self.n_sep_max = ntsepmax                   # limit to how far nt pairs can be in the end
        self.count_percentile_threshold = cperc     # initial thresholding of counts for bins
        self.bin_size_tol = 5                       # nt size tolerance of cluster binning
        self.bins = 0                               # what we eventually settled on for bins after clustering
        self.final_cbin_size = 0                    # what size we eventually got for the bins after clustering
        self.per_dif_threshold = config.per_dif_thr # if a single nt in a pair has 95% of the signal, its a spike

        # Data from binning histogram
        self.freq_histogram_counts = np.array([])   # array of pos:freq dict from histogram containing counts of bins
        self.freq_histogram_edges = np.array([])    # array of pos:freq dict from histogram containing left bin edges
        self.cluster_bin_dictionary = dict()        # dictionary of histogram bin pos:frequency

        # Results from clustering analysis
        self.all_cluster_bin_pairs = list()             # list of all unique cluster bin pairs
        self.filtered_cluster_bin_pairs = list()        # lift of completely filtered bin pairs
        self.best_bin_pair = [(-1, -1), (-1, -1)]       # bin spans that best fulfill the conditions
        self.best_nt_pair = [(-1, 0), (-1, 0)]          # best scoring nucleotide pair with counts
        self.best_nt_pair_dist = 0                      # number of nt apart the pair is
        self.best_nt_pair_sum = 0                       # sum of scores of the best nt pair
        self.graph_nt_stream = config.graph_stream      # amount of nt upstream and downstream shown when drawing graphs
        self.filtered_cluster_bin_dictionary = dict()   # filtered cluster bin:frequency dict for indexing
        self.is_single_signal = 0                       # if one, may be a useless cluster (no buddy)
        self.signal = 0                                 # if a single hit, this is the nt with all the reads

    # private method to count reads in its position:frequency array
    def _count(self):
        pfd = self.pos_freq_dict
        foo = 0
        for key in pfd:
            foo += pfd[key]
        return foo

    # private method to return frequency array for histogram creation
    def _make_freq_histogram(self):

        # enumerate an array based on position and frequency
        data = list()
        for pos in self.pos_array:
            for i in range(0, self.pos_freq_dict[pos]):
                data.append(pos)
        data = np.array(data)
        bins = int((self.cluster_end - self.cluster_start) / self.bin_size)

        # create numpy histogram
        counts, edges = np.histogram(data, bins=bins)

        # check to make sure the histogram bin size is right...sometimes, based on the pos array, it gets a bit small
        bin_size = edges[1] - edges[0]

        # adjust the bin sizes until we meet our bin size tolerance
        while abs(self.bin_size - bin_size) > self.bin_size_tol:
            if self.bin_size > bin_size:
                bins -= 1
            else:
                bins += 1
            # recreate the histogram
            counts, edges = np.histogram(data, bins=bins)
            bin_size = edges[1] - edges[0]

        # now that everything should be good, return the data array and the array of edges along with a
        # dictionary tying the two
        cbin_dict = dict()
        for i in range(0, len(counts)):
            cbin_dict[edges[i]] = counts[i]
        self.final_cbin_size = bin_size
        self.bins = bins

        return counts, edges, cbin_dict

    # private method filters out positions that exceed the given percentile of data
    def _filter_by_read_count(self, cperc):

        cluster_bin_cutoff_perc_val = np.percentile(self.freq_histogram_counts, cperc)

        for pos in self.cluster_bin_dictionary:
            read_count = self.cluster_bin_dictionary[pos]
            if read_count >= cluster_bin_cutoff_perc_val:
                self.filtered_cluster_bin_dictionary[pos] = self.cluster_bin_dictionary[pos]
        return

    # private method makes unique cluster bin pairs based on filtered cluster dictionary
    def _generate_cluster_bin_pairs(self):

        pass_pos_list = list()
        for pos in self.filtered_cluster_bin_dictionary:
            pass_pos_list.append(pos)

        for i in range(0, len(pass_pos_list) - 1):
            this_bin_pos = pass_pos_list[i]
            rest_bin_pos = pass_pos_list[i + 1:]
            for pos in rest_bin_pos:
                self.all_cluster_bin_pairs.append((this_bin_pos, pos))
        return

    # returns an array subset based on data bounds
    def _sub_array(self, pos_start, pos_end):

        a = list()
        for element in self.pos_array:
            if (element >= pos_start) and (element <= pos_end):
                a.append(element)
        return np.array(a)

    # filters bin pairs by class parameters
    def _filter_bin_pairs(self):

        self.filtered_cluster_bin_pairs = self.all_cluster_bin_pairs
        for bin_pair in self.all_cluster_bin_pairs:
            bin1 = bin_pair[0]
            bin2 = bin_pair[1]
            if abs(bin2 - bin1) >= self.c_sep_max or abs(bin2 - bin1) <= self.c_sep_min:
                self.filtered_cluster_bin_pairs.remove(bin_pair)
        return

    # finds the maximally scoring pair of cluster bins
    def _find_max_pair(self):
        bin_count_max, bin_max_pair = 0, (-1, -1)
        for bin_pair in self.filtered_cluster_bin_pairs:
            bin1 = bin_pair[0]
            bin2 = bin_pair[1]
            read_count = self.cluster_bin_dictionary[bin1] + self.cluster_bin_dictionary[bin2]

            if read_count > bin_count_max:
                bin_count_max = read_count
                bin_max_pair = bin_pair

        bin1_lb, bin1_ub = bin_max_pair[0], bin_max_pair[0]+ self.final_cbin_size
        bin2_lb, bin2_ub = bin_max_pair[1], bin_max_pair[1] + self.final_cbin_size
        self.best_bin_pair = ((bin1_lb, bin1_ub), (bin2_lb, bin2_ub))

        return

    # finds the best nucleotides in this cluster region
    def _find_best_nucleotide(self, arr):

        best_nt = -1
        read_max = 0

        for pos in arr:
            if self.pos_freq_dict[pos] > read_max:
                best_nt = pos
                read_max = self.pos_freq_dict[pos]

        return best_nt, read_max

    # returns a frequency array over a position array using pos_freq_dict
    def _make_freq_array(self, pos_array):
        data = list()
        for pos in pos_array:
            for i in range(0, self.pos_freq_dict[pos]):
                data.append(pos)
        return np.array(data)

    # looks at the nt pair and gives and idea of the legitness of the cluster based on class parameters.
    def _assess_nt_pair(self):

        pos1 = self.best_nt_pair[0][0]
        pos2 = self.best_nt_pair[1][0]

        score1 = self.best_nt_pair[0][1]
        score2 = self.best_nt_pair[1][1]

        pos_dif = abs(pos1 - pos2)
        per_dif = 100 * (abs(score1 - score2) / (score1 + score2))

        # now, is one of the nucleotides scoring nearly 100% of the data?
        if per_dif > self.per_dif_threshold:
            self.is_single_signal = 1

        # is the combined read count greater than threshold?
        if score1 + score2 <= config.cluster_count_threshold:
            self.is_single_signal = 1

        # or, are the nts waaay too close or far somehow despite our cluster distance thresholding?
        if (pos_dif >= self.n_sep_max) or (pos_dif <= self.n_sep_min):
            self.is_single_signal = 1

        # if the signal is up, find the offending nucleotide
        if self.is_single_signal == 1:

            if score2 > score1:
                self.signal = (pos2, score2)
            else:
                self.signal = (pos1, score1)
        return

    # uses matplotlib and sns to draw and save an illustration of the histogram data of the suggested inversion cluster
    def draw_inversion_site(self, save_path, figsize=(8, 6), show_fig='n'):

        # generate histogram data over this pos_array
        # first make a sub_array a little upstream and downstream of our positions
        pos1 = self.best_nt_pair[0][0]
        pos2 = self.best_nt_pair[1][0]

        start = pos1 - self.graph_nt_stream
        end = pos2 + self.graph_nt_stream
        arr = self._sub_array(start, end)

        # now make the frequency histogram
        dx = self._make_freq_array(arr)

        # use seabourn to draw our initial density histogram with a gaussian fit
        sns.distplot(dx, bins=30)

        # write vertical lines where we suspect the inversion pair to be
        #plt.axvline(pos1, color='r')
        #plt.axvline(pos2, color='r')

        # label our axes
        plt.xlabel('Nucleotide position')
        plt.ylabel('Read density')

        # draw arrows to annotate the vertical lines so you can see the exact position
        c_start = 'Cluster start: ' + str(pos1)
        c_end = 'Cluster end: ' + str(pos2)
        ymin, ymax = plt.ylim()
        plt.annotate(c_start, xy=(pos1, ymax), xycoords='data', xytext=(0.15, 0.95),
                     textcoords='figure fraction', arrowprops=dict(facecolor='black', shrink=0.05))

        plt.annotate(c_end, xy=(pos2, ymax), xycoords='data', xytext=(0.75, 0.95),
                     textcoords='figure fraction', arrowprops=dict(facecolor='black', shrink=0.05))

        # save our figure before showing
        plt.savefig(save_path)

        # show our figure if desired
        if show_fig == 'y':
            plt.show()

        # clear the figure
        plt.close()

        return

    # analyzes cluster to find the most probable pair of sites of nucleotide inversions
    def detect_inversion_pair(self):

        # make frequency histogram of data
        self.freq_histogram_counts, self.freq_histogram_edges, self.cluster_bin_dictionary = \
            self._make_freq_histogram()

        # filter out bins by percentile
        self._filter_by_read_count(self.count_percentile_threshold)

        # generate cluster bin pairs
        self._generate_cluster_bin_pairs()

        # filter the cluster bin pairs
        self._filter_bin_pairs()

        # if we are out of bin pairs, don't bother...
        if len(self.filtered_cluster_bin_pairs) == 0:
            self.is_single_signal = 1
            self.signal = (self.pos_freq_max, self.freq_max)
            pass

        else:

            # find the best cluster bin pair
            self._find_max_pair()

            # find the best nucleotides in there
            self.best_nt_pair[0] = self._find_best_nucleotide(self._sub_array(
                self.best_bin_pair[0][0], self.best_bin_pair[0][1]))
            self.best_nt_pair[1] = self._find_best_nucleotide(self._sub_array(
                self.best_bin_pair[1][0], self.best_bin_pair[1][1]))

            self.best_nt_pair_sum = self.best_nt_pair[0][1] + self.best_nt_pair[1][1]
            self.best_nt_pair_dist = abs(self.best_nt_pair[0][0] - self.best_nt_pair[1][0])

            # take a glance at the pair to see if we have a true inversion
            self._assess_nt_pair()

            return


