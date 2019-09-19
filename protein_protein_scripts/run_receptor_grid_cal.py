"""
"""

import os
import sys
import glob
import argparse

from _receptor_grid_cal import rec_grid_cal, is_nc_grid_good, get_grid_size_from_lig_rec_crd

parser = argparse.ArgumentParser()
parser.add_argument( "--max_jobs",      type=int, default = 100)
parser.add_argument( "--amber_dir",     type=str, default = "amber")
parser.add_argument( "--coord_dir",     type=str, default = "min")
parser.add_argument( "--out_dir",     type=str, default = "out")

parser.add_argument( "--lj_scale",    type=float, default = 0.6)
parser.add_argument( "--spacing",     type=float, default = 0.5)
parser.add_argument( "--buffer",      type=float, default = 1.0)

parser.add_argument( "--submit",   action="store_true", default=False)

args = parser.parse_args()

LIGAND_INPCRD   = "ligand.inpcrd"
RECEPTOR_INPCRD = "receptor.inpcrd"
RECEPTOR_PRMTOP = "receptor.prmtop"

GRID_NC = "grid.nc"
PDB_OUT = "receptor_trans.pdb"
BOX_OUT = "box.pdb"

def is_running(qsub_file, log_file, nc_file):
    if os.path.exists(qsub_file) and os.path.exists(nc_file) and (not os.path.exists(log_file)):
        return True

    if os.path.exists(qsub_file) and (not os.path.exists(nc_file)) and (not os.path.exists(log_file)):
        return True
    return False

if args.submit:
    this_script = os.path.abspath(sys.argv[0])
    amber_dir = os.path.abspath(args.amber_dir)
    coord_dir = os.path.abspath(args.coord_dir)

    amber_sub_dirs = glob.glob(os.path.join(amber_dir, "*"))
    amber_sub_dirs = [dir for dir in amber_sub_dirs if os.path.isdir(dir)]
    complex_names = [os.path.basename(dir) for dir in amber_sub_dirs]

    box_sizes = {}
    for complex in complex_names:
        rec_inpcrd = os.path.join(coord_dir, complex, RECEPTOR_INPCRD)
        lig_inpcrd = os.path.join(coord_dir, complex, LIGAND_INPCRD)
        box_sizes[complex] = get_grid_size_from_lig_rec_crd(rec_inpcrd, lig_inpcrd, args.buffer)
    complex_names.sort(key=lambda name: box_sizes[name])
    print "Complex    box size"
    for c in complex_names:
        print c, box_sizes[c]

    if args.max_jobs > 0:
        max_jobs = args.max_jobs
    else:
        max_jobs = len(complex_names)
    print "max_jobs = %d"%max_jobs

    job_count = 0
    for complex in complex_names:

        if not os.path.isdir(complex):
            os.makedirs(complex)

        id = complex[:4].lower()
        amber_sub_dir  = os.path.join(amber_dir, complex)
        coor_sub_dir   = os.path.join(coord_dir, complex) 

        out_dir = os.path.abspath(complex)

        qsub_file = os.path.join(out_dir, id+"_grid.job")
        log_file  = os.path.join(out_dir, id+"_grid.log")
        qsub_script = '''#!/bin/bash
#PBS -S /bin/bash
#PBS -o %s '''%log_file + '''
#PBS -j oe
#PBS -l nodes=1:ppn=2,walltime=300:00:00

source /home/tnguye46/opt/module/anaconda.sh
date
python ''' + this_script + \
        ''' --amber_dir ''' + amber_sub_dir + \
        ''' --coord_dir ''' + coor_sub_dir + \
        ''' --out_dir '''   + out_dir + \
        ''' --lj_scale %f'''%args.lj_scale + \
        ''' --spacing %f'''%args.spacing + \
        ''' --buffer %f'''%args.buffer + '''\n'''

        if not is_nc_grid_good(os.path.join(out_dir, GRID_NC)) and not is_running(qsub_file, log_file, 
                                                                os.path.join(out_dir, GRID_NC)): 
            print "Submitting %s"%complex
            open(qsub_file, "w").write(qsub_script)
            os.system("qsub %s" %qsub_file)
            job_count += 1
            if job_count == max_jobs:
                print "Max number of jobs %d reached."%job_count
                break
        else:
            print "Calculation for %s is done"%complex

else:
    prmtop     = os.path.join(args.amber_dir, RECEPTOR_PRMTOP)
    lj_scale   = args.lj_scale
    rec_inpcrd = os.path.join(args.coord_dir, RECEPTOR_INPCRD)
    lig_inpcrd = os.path.join(args.coord_dir, LIGAND_INPCRD)

    spacing    = args.spacing
    buffer     = args.buffer

    grid_out   = os.path.join(args.out_dir, GRID_NC)
    pdb_out    = os.path.join(args.out_dir, PDB_OUT)
    box_out    = os.path.join(args.out_dir, BOX_OUT)

    rec_grid_cal(prmtop, lj_scale, rec_inpcrd, lig_inpcrd, spacing, buffer, grid_out, pdb_out, box_out)


