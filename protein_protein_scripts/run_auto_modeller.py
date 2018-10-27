"""
"""
import os
import sys
import shutil
import glob
import argparse

from _fix_pdb import AddMissing

parser = argparse.ArgumentParser()
parser.add_argument( "--pdb_dir",     type=str, default = "pdbs" )
parser.add_argument( "--pdb",     type=str, default = "xxx.pdb" )
parser.add_argument( "--submit",   action="store_true", default=False )

args = parser.parse_args()

if args.submit:
    
    this_script = os.path.abspath(sys.argv[0])
    pwd = os.getcwd()

    working_dirs = {}
    pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
    for pdb in pdb_files:
        id = os.path.basename(pdb)[:-4]

        if not os.path.isdir(id):
            os.makedirs(id)
        working_dirs[id] = os.path.join(pwd, id)
        shutil.copy(pdb, working_dirs[id])

    for id, dir in working_dirs.items():

        pdb_file_name = id+".pdb"

        qsub_file = os.path.join(dir, id+"_modeller.job")
        log_file  = os.path.join(dir,id+"_modeller.log" )
        qsub_script = '''
#!/bin/bash
#PBS -S /bin/bash
#PBS -o %s '''%log_file + '''
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=300:00:00

source /home/tnguye46/opt/module/anaconda.sh
cd ''' + dir + '''
date
python ''' + this_script + \
        ''' --pdb ''' + pdb_file_name
        open( qsub_file, "w" ).write( qsub_script )
        os.system( "qsub %s" %qsub_file )

else:
    fix = AddMissing(args.pdb)
    fix.write_dummy_alignments()
    fix.do_auto_align()
    fix.do_auto_modeller()
    fix.write_vmd_script()
