"""
Configuration module for Cluster Detection

Author: Jacob Bourgeois
Email: jacob.bourgeois@tufts.edu
Organization: Tufts Sackler School of Biomedical Sciences; Camilli Lab

License: 3-clause BSD
"""

import inv_config

"""
Your email goes here. This is used to tell NCBI who is asking for entrez genbank data
According to NCBI guidelines, please keep queries to less than 100 files a day, and try to use on non-peak hours.
More details on NCBI's site. If you abuse this, then the NCBI denies request based on
your ip if you don't give an email!

Modified to include from inv_config.
"""
email = inv_config.email
api_key = inv_config.api_key

# Cluster Detection parameters

"""
nbin_size refers to how many nucleotides each bin is in the initial SOR read histogram
used to filter out potential clusters. Too large, and signals may become mixed together.
Too small, and signals may be split. This value should be changed according to how large
of a inversion pair you are looking for.
Be aware that changing this will alter the ideal threshold for automatic detection of clusters
"""
nbin_size = 5000

"""
"cbin_cutoff" refers to what percentile of read count an individual bin within a cluster
must have to be considered a site for potential inversion nucleotides. I'd keep this at 98

"cbin_size" refers to how many nucleotides is in each bin within a potential cluster. 10
has been working for me so far, but you can adjust is larger or smaller based on how messy
and how close together individual peaks are. For example, if I had a super messy peak at
position 10000, and another at 10002, it may be nice to have a cbin_size of at least two
so these are combined together, and not mistaken for an inversion pair of 2 nt.
"""
cbin_cutoff = 98
cbin_size = 10

"""=
"cluster_min_sep" and "cluster_max_sep" refer to how many nucleotides should be expected
between the two highest-scoring bins within a potential cluster - so, at least cluster_min_sep
and no more than cluster_max_sep. This again helps the program know what size of inversion
clusters you are looking for.
"""
cluster_min_sep = 100
cluster_max_sep = 2000

"""
"ntpair_min_sep" and "ntpair_max_sep" are like the above, but look at the distances of the
final exact nt positions. Sometimes, when two clusters are identified, the best positions
are at the respective beginning and end of a cluster. So, this makes sure you have the size
of inversion pair you want at the end.
"""
ntpair_min_sep = 50
ntpair_max_sep = 1000

"""
read_count_threshold screens initial SOR reads - if you have less than this, don't count it.
Keeping this at zero works for me - too high and and sometimes histogram binning fails
"""
read_count_threshold = 0

"""
per_dif_thr sets how different inversion pairs must be. If a nt in a read pair exceeds this percent
threshold, then the cluster is thrown out"""
per_dif_thr = 90

"""
cluster_count_threshold means there must be this many reads in the inversion pair to be counted.
"""
cluster_count_threshold = 10

"""
when using automatic SOR thresholding, only count bins containing more than the this percentile
"""
automatic_thresholding_perc = 95
is_automatic = True

# Cluster gene analysis parameters

"""
When aligning the genes, cap the nearby orfs at max_genes and only look graph_stream nts downstream
and upstream
"""
max_genes = 6
graph_stream = 1000

"""
Variables for detection of inverted repeats. Basically takes a seed sequence of size inv_nt_seed_size
and scans for a match. If it does, it continues adding to the seed sequence and checking. It'll allow
a inv_nt_mismatch_tol amount of mismatches before quitting. The seed sequence will originate at the beginning
of the switch offset by start_pos (ie. -1 means start one nt before the beginning of the detected inversion) and
stop at end_pos. Each seed will try attempts times to find homology before moving to the next seed.

"""
inv_nt_seed_size = 6
inv_nt_mismatch_tol = 3
start_pos = -5
end_pos = 15
attempts = 30


"""
figsize is a tuple that defines the size of the saved graphs of SOR and cluster images in inches
"""
figsize = (10, 8)
