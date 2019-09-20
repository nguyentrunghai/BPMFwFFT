"""
run fft sampling
"""
from __future__ import print_function

import sys
import os
import glob
import argparse

from _receptor_grid_cal import is_nc_grid_good
from _fft_sampling import sampling, is_sampling_nc_good
from _receptor_grid_cal import get_grid_size_from_nc

parser = argparse.ArgumentParser()
parser.add_argument("--max_jobs",                      type=int, default=50)
parser.add_argument("--amber_dir",                     type=str, default="amber")
parser.add_argument("--coord_dir",                     type=str, default="min")
parser.add_argument("--grid_dir",                      type=str, default="grid")
parser.add_argument("--lig_ensemble_dir",              type=str, default="rotation")

parser.add_argument("--energy_sample_size_per_ligand", type=int, default=1000)
parser.add_argument("--nr_lig_conf",                   type=int, default=500)
parser.add_argument("--start_index",                    type=int, default=0)

parser.add_argument("--out_dir",                       type=str, default="out")

parser.add_argument("--lj_scale",                      type=float, default=0.8)
parser.add_argument("--submit",   action="store_true", default=False)
args = parser.parse_args()

RECEPTOR_INPCRD = "receptor.inpcrd"
RECEPTOR_PRMTOP = "receptor.prmtop"

LIGAND_INPCRD = "ligand.inpcrd"
LIGAND_PRMTOP = "ligand.prmtop"

LIG_COOR_NC = "rotation.nc"

GRID_NC = "grid.nc"
FFT_SAMPLING_NC = "fft_sample.nc"


def is_running(qsub_file, log_file, nc_file):
    if os.path.exists(qsub_file) and (not os.path.exists(nc_file)) and (not os.path.exists(log_file)):
        return True
    if os.path.exists(qsub_file) and os.path.exists(nc_file) and (not os.path.exists(log_file)):
        return True
    return False


if args.submit:
    this_script = os.path.abspath(sys.argv[0])
    amber_dir = os.path.abspath(args.amber_dir)
    coord_dir = os.path.abspath(args.coord_dir)
    grid_dir = os.path.abspath(args.grid_dir)
    lig_ensemble_dir = os.path.abspath(args.lig_ensemble_dir)

    complex_names = glob.glob(os.path.join(grid_dir, "*"))
    complex_names = [os.path.basename(d) for d in complex_names if os.path.isdir(d)]

    complex_names = [c for c in complex_names if is_nc_grid_good(os.path.join(grid_dir, c, GRID_NC))]

    grid_sizes = {}

    for complex_name in complex_names:
        grid_sizes[complex_name] = get_grid_size_from_nc(os.path.join(grid_dir, complex_name, GRID_NC))
    complex_names.sort(key=lambda name: grid_sizes[name])
    print("Complex   grid size")

    for complex_name in complex_names:
        print(complex_name, grid_sizes[complex_name])

    pwd = os.getcwd()
    complex_names = [c for c in complex_names if not is_sampling_nc_good(
        os.path.join(pwd, c, FFT_SAMPLING_NC), args.nr_lig_conf)]
    
    if args.max_jobs > 0:
        max_jobs = args.max_jobs
    else:
        max_jobs = len(complex_names)
    print("max_jobs = %d" % max_jobs)

    job_count = 0
    for complex_name in complex_names:

        if not os.path.isdir(complex_name):
            os.makedirs(complex_name)

        idx = complex_name[:4].lower()
        amber_sub_dir = os.path.join(amber_dir, complex_name)
        coor_sub_dir = os.path.join(coord_dir, complex_name)
        grid_sub_dir = os.path.join(grid_dir, complex_name)
        lig_ensemble_sub_dir = os.path.join(lig_ensemble_dir, complex_name)

        out_dir = os.path.abspath(complex_name)
        qsub_file = os.path.join(out_dir, idx+"_fft.job")
        log_file = os.path.join(out_dir, idx+"_fft.log")
        qsub_script = '''#!/bin/bash
#PBS -S /bin/bash
#PBS -o %s '''%log_file + '''
#PBS -j oe
#PBS -l nodes=1:ppn=4,walltime=300:00:00

source /home/tnguye46/opt/module/anaconda.sh
date
python ''' + this_script + \
        ''' --amber_dir ''' + amber_sub_dir + \
        ''' --coord_dir ''' + coor_sub_dir + \
        ''' --grid_dir ''' + grid_sub_dir + \
        ''' --lig_ensemble_dir ''' + lig_ensemble_sub_dir + \
        ''' --out_dir '''   + out_dir + \
        ''' --lj_scale %f'''%args.lj_scale + \
        ''' --nr_lig_conf %d '''%args.nr_lig_conf + \
        ''' --start_index %d'''%args.start_index + \
        ''' --energy_sample_size_per_ligand %d '''%args.energy_sample_size_per_ligand + '''\n'''

        fft_sampling_nc_file = os.path.join(out_dir, FFT_SAMPLING_NC)
        if not is_running(qsub_file, log_file, fft_sampling_nc_file):

            if os.path.exists(fft_sampling_nc_file):
                print("remove file " + fft_sampling_nc_file)
                os.system("rm "+fft_sampling_nc_file)

            if os.path.exists(log_file):
                print("remove file " + log_file)
                os.system("rm "+log_file)

            print("Submitting %s" % complex_name)
            open(qsub_file, "w").write(qsub_script)
            os.system("qsub %s" % qsub_file)
            job_count += 1
            if job_count == max_jobs:
                print("Max number of jobs %d reached." % job_count)
                break
else:
    rec_prmtop = os.path.join(args.amber_dir, RECEPTOR_PRMTOP)
    lj_sigma_scal_fact = args.lj_scale
    rec_inpcrd = os.path.join(args.amber_dir, RECEPTOR_INPCRD)

    grid_nc_file = os.path.join(args.grid_dir, GRID_NC)

    lig_prmtop = os.path.join(args.amber_dir, LIGAND_PRMTOP)
    lig_inpcrd = os.path.join(args.amber_dir, LIGAND_INPCRD)

    lig_coor_nc = os.path.join(args.lig_ensemble_dir, LIG_COOR_NC)
    nr_lig_conf = args.nr_lig_conf
    start_index = args.start_index

    energy_sample_size_per_ligand = args.energy_sample_size_per_ligand
    output_nc = os.path.join(args.out_dir, FFT_SAMPLING_NC)

    sampling(rec_prmtop, lj_sigma_scal_fact,
             rec_inpcrd, grid_nc_file,
             lig_prmtop, lig_inpcrd,
             lig_coor_nc, nr_lig_conf, start_index,
             energy_sample_size_per_ligand,
             output_nc)

