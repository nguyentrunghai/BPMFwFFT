"""
run this script to download all pdb listed in the affinity benchmark set
"""
from __future__ import print_function

import argparse
import os
import urllib

from _affinity_data import AffinityData

AFFINITY_DATA_FILES = ["affinity_v1.tsv",  "affinity_v2.tsv"]

parser = argparse.ArgumentParser()
parser.add_argument("--affinity_data_dir", type=str, default="affinity")

args = parser.parse_args()

PDB_URL = "http://www.rcsb.org/pdb/files"
NOT_AVAIL_MESSAGE = "the requested file is not available"


def download_coord(id):
    """
    download pdb with id and save it with id.pdb (lower case)
    """
    pdb_file = id.lower() + ".pdb"
    url = os.path.join(PDB_URL, pdb_file)
    data = urllib.urlopen(url).read()
    if NOT_AVAIL_MESSAGE in data:
        raise RuntimeError(data)
    open(pdb_file, "w").write(data)
    return None


affinity_data_files = [os.path.join(args.affinity_data_dir, file) for file in AFFINITY_DATA_FILES]

aff_data = AffinityData(affinity_data_files)
pdb_ids = aff_data.unique_pdb_ids()

for pdb_id in pdb_ids:
    print("downloading " + pdb_id)
    download_coord(pdb_id)

print("%d files downloaded" % len(pdb_ids))
print("DONE")
