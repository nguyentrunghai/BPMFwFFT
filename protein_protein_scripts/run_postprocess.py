"""
run post processing
"""

from __future__ import print_function

import sys
import os
import glob
import argparse

from _fft_sampling import is_sampling_nc_good, parse_nr_ligand_confs 
from _postprocess import post_process


parser = argparse.ArgumentParser()
parser.add_argument("--amber_dir", type=str, default="amber")
parser.add_argument("--sampling_dir", type=str, default="fft_sampling")

parser.add_argument("--nr_resample", type=int, default=100)

parser.add_argument("--out_dir", type=str, default="out")
parser.add_argument("--submit", action="store_true", default=False)
args = parser.parse_args()

RECEPTOR_PRMTOP = "receptor.prmtop"
LIGAND_PRMTOP = "ligand.prmtop"
COMPLEX_PRMTOP = "complex.prmtop"

FFT_SAMPLING_NC = "fft_sample.nc"

REC_PDB_OUT = "receptor_trans.pdb"
LIG_PDB_OUT = "ligand_resampled.pdb"
BPMF_OUT = "bpmf.pkl"


def is_sampling_good(sampling_dir):
    complex_name = sampling_dir.split("/")[-1]
    idx = complex_name[:4].lower()
    submit_file = idx + "_fft.job"
    nr_lig_confs =  parse_nr_ligand_confs(os.path.join(sampling_dir, submit_file))
    if nr_lig_confs is None:
        return False

    nc_sampling_file = os.path.join(sampling_dir, FFT_SAMPLING_NC)
    return is_sampling_nc_good(nc_sampling_file, nr_lig_confs)

if args.submit:
    this_script = os.path.abspath(sys.argv[0])
    amber_dir = os.path.abspath(args.amber_dir)
    sampling_dir = os.path.abspath(args.sampling_dir)

    complex_names = glob.glob(os.path.join(sampling_dir, "*"))
    complex_names = [os.path.basename(d) for d in complex_names if os.path.isdir(d)]
    complex_names = [c for c in complex_names if is_sampling_good(os.path.join(sampling_dir, c))]
    print(complex_names)

    for complex_name in complex_names:
        if not os.path.isdir(complex_name):
            os.makedirs(complex_name)

        idx = complex_name[:4].lower()
        amber_sub_dir = os.path.join(amber_dir, complex_name)
        sampling_sub_dir = os.path.join(sampling_dir, complex_name)
        out_dir = os.path.abspath(complex_name)

        qsub_file = os.path.join(out_dir, idx+"_post.job")
        log_file = os.path.join(out_dir, idx+"_post.log")
        qsub_script = '''#!/bin/bash
#PBS -S /bin/bash
#PBS -o %s '''%log_file + '''
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=300:00:00
module load ambertools/14
source /home/tnguye46/opt/module/anaconda.sh
date
python ''' + this_script + \
        ''' --amber_dir ''' + amber_sub_dir + \
        ''' --sampling_dir ''' + sampling_sub_dir +\
        ''' --out_dir '''   + out_dir + \
        ''' --nr_resample %d'''%args.nr_resample + '''\n'''

        bpmf_out = os.path.join(out_dir, BPMF_OUT)
        if not os.path.exists(bpmf_out):
            open(qsub_file, "w").write(qsub_script)
            print("Submiting " + qsub_file)
            os.system("qsub %s" % qsub_file)
else:
    rec_prmtop = os.path.join(args.amber_dir, RECEPTOR_PRMTOP)
    lig_prmtop = os.path.join(args.amber_dir, LIGAND_PRMTOP)
    complex_prmtop = os.path.join(args.amber_dir, COMPLEX_PRMTOP)

    sampling_nc_file = os.path.join(args.sampling_dir, FFT_SAMPLING_NC)
    nr_resampled_complexes = args.nr_resample

    sander_tmp_dir = args.out_dir

    rec_pdb_out = os.path.join(args.out_dir, REC_PDB_OUT)
    lig_pdb_out = os.path.join(args.out_dir, LIG_PDB_OUT)
    bpmf_pkl_out = os.path.join(args.out_dir, BPMF_OUT )

    post_process(rec_prmtop, lig_prmtop, complex_prmtop, sampling_nc_file, 
            nr_resampled_complexes, 
            sander_tmp_dir,
            rec_pdb_out, lig_pdb_out, bpmf_pkl_out)

