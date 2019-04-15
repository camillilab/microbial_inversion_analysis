#! usr/bin/python

"""
Primary executable script for cluster detection.

Navigate to the directory with all the scripts, type in:
python detect_inversion_clusters.py

Requires Python 3.6+

Please see config.py to define running parameters.

Author: Jacob Bourgeois
Email: jacob.bourgeois@tufts.edu
Organization: Tufts Sackler School of Biomedical Sciences; Camilli Lab

License: 3-clause BSD

"""

from cluster_tools import Bug
from write_tools import *
import config
import os
import csv
import shutil
import argparse
import inv_config


def main(oname):

    # Define working directory as current working directory
    cwd = os.getcwd()

    # define data paths
    sor_path = inv_config.sor_save_path
    sclip_path = inv_config.sclip_save_path
    output_path = inv_config.cluster_save_path

    output_path = os.path.join(output_path, oname)

    gbk_path = inv_config.acc_save_path

    # Make sure the paths exist
    checkdirs([output_path, gbk_path])

    # load accession nums
    acc_file = inv_config.acc_list_path
    acc_nums = list()
    with open(acc_file, 'r') as f:
        for line in f:
            acc_nums.append(line.split('\n')[0])
    print("Loaded {0} accession numbers.".format(len(acc_nums)))

    for acc_num in acc_nums:

        # create relevant output directories. If it already exists, overwrite
        working = os.path.join(output_path, acc_num)
        if os.path.exists(working):
            print("Cleaning results for {0}...".format(acc_num), end='')
            shutil.rmtree(working, ignore_errors=True)
            print("Done.")
        checkdirs([working])

        gene_diagrams_path = os.path.join(working, "Gene Diagrams")
        cluster_graph_path = os.path.join(working, "Cluster Graphs")
        checkdirs([gene_diagrams_path, cluster_graph_path])

        # create a bug instance
        bug = Bug(acc_num=acc_num)

        # get entrez data and load onto bug instance
        gbkfile = os.path.join(gbk_path, acc_num+'.gb')
        bug.get_entrez_data(filename=gbkfile, email=config.email)

        # load SOR and sCLIP data onto bug instance
        sor_data = os.path.join(sor_path, acc_num+'_sor.csv')
        sclip_data = os.path.join(sclip_path, acc_num+'_sclip.csv')

        # Sometimes, I get an error here where a non-zero array returns. I believe this may have to do with nwidth
        # being too large, so I get 0-size bins. Set up a try statement to avoid intermittent crashing.

        try:
            bug.load_SOR(sor_data)
            bug.load_sCLIP(sclip_data)

            # make graphical threshold for SOR data to determine signals
            clusters = bug.make_SOR_interactive_graphical_threshold(automatic=config.is_automatic)
            save_path = os.path.join(working, "{0} SOR Histogram".format(acc_num))
            bug.save_SOR_thresholding(save_path=save_path, figsize=config.figsize)

            # for each signal...
            c = 0
            for cluster in clusters:

                # try to detect an inversion nt pair
                print("Processing signal ({0}, {1})...".format(cluster.cluster_start, cluster.cluster_end))
                cluster.detect_inversion_pair()

                # if the signal is a clear inversion, add to the instance
                if cluster.is_single_signal != 1:
                    c += 1
                    bug.inversions.append(cluster)

                    # create a diagram of this cluster
                    save_path = os.path.join(cluster_graph_path, '{0}_{1}'.format(acc_num, c))
                    cluster.draw_inversion_site(save_path=save_path, figsize=config.figsize)

            print("{0} inversion pairs detected.".format(len(bug.inversions)))

            # for each inversion pair
            c = 0
            for inversion in bug.inversions:

                inv_pair = inversion.best_nt_pair
                locs = (inv_pair[0][0], inv_pair[1][0])

                # align entrez data to clusters
                loci, is_intergenic = bug.align_inversion_to_genes(locs, max_genes=config.max_genes)

                # flag this inversion as intergenic in the instance
                bug.is_intergenic.append(is_intergenic)

                # try to detect inverted repeats
                bug.rcscores.append(bug.align_inverted_repeats(locs, seed=config.inv_nt_seed_size, start=config.start_pos,
                                                               end=config.end_pos, tol=config.inv_nt_mismatch_tol,
                                                               tries=config.attempts))

                # create a gene diagram using GenomeDiagram
                c += 1
                save = os.path.join(gene_diagrams_path, '{0}_{1}.pdf'.format(acc_num, c))
                print("Drawing gene diagram...", end='')
                bug.draw_gene_diagram(locs, loci, save_path=save)
                print("Done.")

            # now write an analysis file
            analysis_file = os.path.join(working, "{0}_cluster_analysis.csv".format(acc_num))
            inv_file = os.path.join(working, "{0}_inversions.csv".format(acc_num))
            bug.make_analysis_file(save_path=analysis_file)
            bug.make_inversion_file(save_path=inv_file)
        except ValueError:
            print("Encountered zero-size array! Skipping...")


# make executable
if __name__ == "__main__":

    # ARGPARSE
    parser = argparse.ArgumentParser()

    # REQUIRED ARGUMENTS
    parser.add_argument('save_path', help='Path to save results to.')
    args = parser.parse_args()

    main(oname=args.save_path)


