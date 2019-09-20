"""
to run minimization for loop structures
"""

from __future__ import print_function

import os
import sys
import glob
import argparse

from _loop_energy_minimize import openMM_minimize


parser = argparse.ArgumentParser()
parser.add_argument("--amber_dir", type=str, default="amber")
parser.add_argument("--minimize_all",  action="store_true", default=False)

parser.add_argument("--submit",  action="store_true", default=False)

parser.add_argument("--in_dir",  type=str, default="inp")
parser.add_argument("--out_dir", type=str, default="out")
parser.add_argument("--out_pdb", type=str, default="complex.pdb")

args = parser.parse_args()


if args.submit:

    this_script = os.path.abspath(sys.argv[0])

    amber_dir = os.path.abspath(args.amber_dir)
    amber_sub_dirs = glob.glob(os.path.join(amber_dir, "*"))
    amber_sub_dirs = [dir for dir in amber_sub_dirs if os.path.isdir(dir)]
    complex_names = [os.path.basename(d) for d in amber_sub_dirs]

    for complex in complex_names:
        if not os.path.isdir(complex):
            os.makedirs(complex)

    for complex in complex_names:
        print("Minimizing", complex)
        id = complex[:4].lower()
        in_dir  = os.path.join(amber_dir, complex)
        out_dir = os.path.abspath(complex)

        qsub_file = os.path.join(out_dir, id+"_min.job")
        log_file  = os.path.join(out_dir, id+"_min.log")
        qsub_script = '''#!/bin/bash
#PBS -S /bin/bash
#PBS -o %s '''%log_file + '''
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=300:00:00

source /home/tnguye46/opt/module/anaconda.sh
date
python ''' + this_script + \
        ''' --in_dir ''' + in_dir + \
        ''' --out_dir ''' + out_dir + \
        ''' --out_pdb ''' + args.out_pdb + '''\n'''
        open(qsub_file, "w").write( qsub_script)
        os.system("qsub %s" %qsub_file)

else:
    openMM_minimize( args.in_dir, args.out_dir, out_pdb=args.out_pdb, minimize_all=args.minimize_all)
