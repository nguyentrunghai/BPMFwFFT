"""
to run FFT sampling for an ensemble of ligand configurations and a rigid receptor
"""
from __future__ import print_function

import sys
import os
import argparse

import numpy as np
import netCDF4 as nc

# change this 
sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from fft_sampling import Sampling_PL

parser = argparse.ArgumentParser()

parser.add_argument( "--receptor_prmtop",           type=str, default="receptor.prmtop")
parser.add_argument( "--receptor_inpcrd",           type=str, default="receptor.inpcrd")
parser.add_argument( "--lj_scale_factor",           type=float, default=1.0)
parser.add_argument( "--bsite",                     type=str, default="measured_binding_site.py")
parser.add_argument( "--grid_nc_file",              type=str, default="grid.nc")

parser.add_argument( "--ligand_prmtop",             type=str, default="ligand.prmtop")
parser.add_argument( "--ligand_inpcrd",             type=str, default="ligand.inpcrd")

parser.add_argument( "--ligand_traj_nc_file",       type=str, default="traj.nc")
parser.add_argument( "--number_ligand_samples",     type=int, default=1000)
parser.add_argument( "--randomly_sampling_ligand",  action="store_true", default=False)

parser.add_argument( "--number_translations_per_ligand_sample_stored", type=int, default=1000)

parser.add_argument( "--nc_out_file",   type=str, default="fft_sample.nc")

args = parser.parse_args()


def is_sampling_done(nc_file, number_ligand_samples):
    if not os.path.exists(nc_file):
        return False

    nc_handle = nc.Dataset(nc_file, "r")
    cond1 = nc_handle.variables["lig_positions"][:].shape[0] == number_ligand_samples
    if not cond1:
        return False

    cond2 = type(nc_handle.variables["lig_positions"][:]) == np.ndarray
    if not cond2:
        return False

    return True


if not is_sampling_done(args.nc_out_file, args.number_ligand_samples):

    lig_nc_handle = nc.Dataset(args.ligand_traj_nc_file, "r")
    if args.number_ligand_samples > lig_nc_handle.variables["positions"].shape[0]:
        raise Exception("The number of ligand samples stored in " + args.ligand_traj_nc_file + 
                " is less than %d"%args.number_ligand_samples )

    if not args.randomly_sampling_ligand:
        ligand_samples = lig_nc_handle.variables["positions"][0 : args.number_ligand_samples]
    else:
        sel_ind = np.random.choice(lig_nc_handle.variables["positions"].shape[0], size=args.number_ligand_samples, replace=False)
        ligand_samples = lig_nc_handle.variables["positions"][sel_ind]

    lig_nc_handle.close()

    sampler = Sampling_PL(  args.receptor_prmtop, args.lj_scale_factor, args.receptor_inpcrd,
                            args.bsite, args.grid_nc_file,
                            args.ligand_prmtop, args.ligand_inpcrd,
                            ligand_samples,
                            args.number_translations_per_ligand_sample_stored,
                            args.nc_out_file,
                            temperature=300.)
    sampler.run_sampling()

    print("Sampling Done")

else:
    print(args.nc_out_file + " is good, nothing to be done!")


