"""
"""

import os
import sys
import glob
import argparse

from _loop_energy_minimize import _read_natoms
sys.path.append("../bpmfwfft")
from rotation import random_gen_rotation, systematic_gen_rotation

parser = argparse.ArgumentParser()
parser.add_argument("--coord_dir",  type=str, default="min")
parser.add_argument("--inpcrd",     type=str, default="ligand.inpcrd")
parser.add_argument("--nc",         type=str, default="rotation.nc")
parser.add_argument("--type",       type=str, default="random")
parser.add_argument("--count",      type=int, default=2000)

parser.add_argument("--submit",   action="store_true", default=False)

args = parser.parse_args()

LIGAND_INPCRD = "ligand.inpcrd"
OUTPUT_NC = "rotation.nc"
RAMDOM = "random"
SYSTEMATIC = "systematic"

if args.submit:
    this_script = os.path.abspath(sys.argv[0])
    coord_dir = os.path.abspath(args.coord_dir)

    coord_sub_dirs = glob.glob(os.path.join(coord_dir, "*"))
    coord_sub_dirs = [d for d in coord_sub_dirs if os.path.isdir(d)]
    complex_names = [os.path.basename(d) for d in coord_sub_dirs]

    natoms = {}
    for complex in complex_names:
        natoms[complex] = _read_natoms(os.path.join(coord_dir, complex, LIGAND_INPCRD))
    complex_names.sort(key=lambda name: natoms[name])

    for complex in complex_names:
        if not os.path.isdir(complex):
            os.makedirs(complex)

    for complex in complex_names:
        id = complex[:4].lower()
        ligand_inpcrd = os.path.join(coord_dir, complex, LIGAND_INPCRD)

        out_dir = os.path.abspath(complex)
        output_nc = os.path.join(out_dir, OUTPUT_NC)

        qsub_file = os.path.join(out_dir, id+"_rot.job")
        log_file  = os.path.join(out_dir, id+"_rot.log")

        qsub_script = '''#!/bin/bash
#PBS -S /bin/bash
#PBS -o %s '''%log_file + '''
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=300:00:00

source /home/tnguye46/opt/module/anaconda.sh
date

python ''' + this_script + \
        ''' --inpcrd ''' + ligand_inpcrd + \
        ''' --type '''   + args.type + \
        ''' --count %d '''%args.count + \
        ''' --nc ''' + output_nc + "\n"
        open( qsub_file, "w" ).write( qsub_script )
        if (not os.path.exists(output_nc)) or (os.path.getsize(output_nc) == 0):
            print("Submitting %s"%complex)
            os.system("qsub %s" %qsub_file)
        else:
            print("Calculation for %s is done"%complex)
else:
    if args.type == RAMDOM:
        random_gen_rotation(args.inpcrd, args.count, args.nc)
    elif args.type == SYSTEMATIC:
        systematic_gen_rotation(args.inpcrd, args.count, args.nc)
    else:
        raise RuntimeError("Unknown type")


