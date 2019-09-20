"""
to generate AMBER topology and coordinate files for bound complexes and their binding partners
"""
from __future__ import print_function

import os
import argparse

from _affinity_data import AffinityData
from _chains_combine import write_b_receptor_ligand_pdbs
from _amber_tleap import generate_prmtop

parser = argparse.ArgumentParser()
parser.add_argument("--affinity_dir", type=str, default="affinity")

parser.add_argument("--ions_cofactors_dir", type=str, default="ions_cofactors")

parser.add_argument("--cofactors_frcmod_dir", type=str, default="cofactors_frcmod")

parser.add_argument("--modeller_dir", type=str, default="modeller")

args = parser.parse_args()

AFFINITY_FILES = ["affinity_v1.tsv",  "affinity_v2.tsv"]

affinity_data_files = [os.path.join(args.affinity_dir, file) for file in AFFINITY_FILES]

TER_CUTOFF = 5
LOOP_CUTOFF = 15

if not os.path.isdir(args.modeller_dir):
    raise RuntimeError("% does not exist" % args.modeller_dir)

if not os.path.isdir(args.ions_cofactors_dir):
    raise RuntimeError("% does not exist" % args.ions_cofactors_dir)

aff = AffinityData(affinity_data_files)
b_complexes = aff.get_bound_complexes()

write_b_receptor_ligand_pdbs(b_complexes, args.modeller_dir, args.ions_cofactors_dir,
                             ter_cutoff=TER_CUTOFF, loop_cutoff=LOOP_CUTOFF)
generate_prmtop(args.cofactors_frcmod_dir)
print("Done")

