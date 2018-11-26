"""
run this script to download all pdb listed in the affinity benchmark set
"""
import argparse
import os

from _affinity_data import affinity_data
from _pdb import download_coord

AFFINITY_DATA_FILES = ["affinity_v1.tsv",  "affinity_v2.tsv"]

parser = argparse.ArgumentParser()
parser.add_argument( "--affinity_data_dir",          type=str,   default = "affinity" )

args = parser.parse_args()

affinity_data_files = [os.path.join(args.affinity_data_dir, file) for file in AFFINITY_DATA_FILES]

aff_data = affinity_data(affinity_data_files)
pdb_ids = aff_data.unique_pdb_ids()

for id in pdb_ids:
    print "downloading " + id
    download_coord(id)

print "%d files downloaded"%len(pdb_ids)
print "Done"
