"""
to generate AMBER topology and coordinate files for bound complexes and their binding partners
"""
import os
import argparse

from _affinity_data import AffinityData
from _chains_combine import write_b_receptor_ligand_pdbs
from _amber_tleap import generate_prmtop

parser = argparse.ArgumentParser()
parser.add_argument( "--ions_cofactors_dir",     type=str, default = "ions_cofactors")
args = parser.parse_args()

AFFINITY_DATA_FILES = ["affinity_v1.tsv",  "affinity_v2.tsv"]
AFFINITY_DATA_DIR   = "/home/tnguye46/protein_binding/data_prep/affinity"
AFFINITY_DATA_FILES = [ os.path.join(AFFINITY_DATA_DIR, file) for file in AFFINITY_DATA_FILES ]

MODELLER_DIR = "/home/tnguye46/protein_binding/data_prep/coordinates/modeller"
IONS_COFACTORS_DIR = os.path.abspath(args.ions_cofactors_dir)
TER_CUTOFF   = 5
LOOP_CUTOFF  = 15

if not os.path.isdir(MODELLER_DIR):
    raise RuntimeError("% does not exist"%MODELLER_DIR)

if not os.path.isdir(IONS_COFACTORS_DIR):
    raise RuntimeError("% does not exist"%IONS_COFACTORS_DIR)

aff = AffinityData(AFFINITY_DATA_FILES)
b_complexes = aff.get_bound_complexes()

write_b_receptor_ligand_pdbs(b_complexes, MODELLER_DIR, IONS_COFACTORS_DIR, 
                                ter_cutoff=TER_CUTOFF, loop_cutoff=LOOP_CUTOFF)
generate_prmtop()
print "Done"

