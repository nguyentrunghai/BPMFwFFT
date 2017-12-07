"""
run postprocessing FFT sampling data in implicit solvent models 
to estimate binding potential of mean force
"""

import sys
import argparse

# change this 
sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from postprocess import PostProcess_PL

parser = argparse.ArgumentParser()
parser.add_argument( "--receptor_prmtop",           type=str, default="receptor.prmtop")
parser.add_argument( "--ligand_prmtop",             type=str, default="ligand.prmtop")
parser.add_argument( "--complex_prmtop",            type=str, default="complex.prmtop")

parser.add_argument( "--fft_sampling_nc_file",      type=str, default="fft_sample.nc")

parser.add_argument( "--solvent_phases",            type=str, default="OpenMM_OBC2")
parser.add_argument( "--nr_resampled_complexes",    type=int, default=100)

parser.add_argument( "--sander_tmp_out_dir",        type=str, default="./")

parser.add_argument( "--receptor_pdb_out",          type=str, default="rec.pdb")
parser.add_argument( "--ligand_pdb_out",            type=str, default="lig.pdb")
parser.add_argument( "--results_out",               type=str, default="bpmf.pkl")

args = parser.parse_args()

solvent_phases = args.solvent_phases.split()
randomly_translate_complex = False
temperature = 300.

post_pro = PostProcess_PL(args.receptor_prmtop, args.ligand_prmtop, args.complex_prmtop,
                            args.fft_sampling_nc_file, 
                            solvent_phases,
                            args.nr_resampled_complexes,
                            randomly_translate_complex,
                            temperature,
                            args.sander_tmp_out_dir)

post_pro.write_rececptor_pdb(args.receptor_pdb_out)
post_pro.write_resampled_ligand_pdb(args.ligand_pdb_out)
post_pro.pickle_bpmf(args.results_out)

