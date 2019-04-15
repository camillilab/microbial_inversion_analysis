"""

File I/O tools. Helps when checking for directories and appending to csvs.

Author: Jacob Bourgeois
Email: jacob.bourgeois@tufts.edu
Organization: Tufts Sackler School of Biomedical Sciences; Camilli Lab

License: 3-clause BSD

"""

import os
import csv


# append_to_csv takes a data tuple and appends to some csv file.
def append_to_csv(data, output):

    with open(output, 'a') as o:
        writer = csv.writer(o)
        writer.writerow(data)
    return


def checkdirs(paths):
    for path in paths:
        if not os.path.exists(path):
            print("Making directory at {0}...".format(path))
            os.makedirs(path)
    return

