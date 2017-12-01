"""
rum temperature hamiltonian replica exchange simulation
"""

import os
import sys
import argparse

import numpy as np
import netCDF4

# change this 
sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from md_openmm import OpenMM_TREMD

parser = argparse.ArgumentParser()
parser.add_argument( "--ligand_prmtop",             type=str, default="ligand.prmtop")
parser.add_argument( "--ligand_inpcrd",             type=str, default="ligand.inpcrd")

parser.add_argument( "--phase",                     type=str, default = "OpenMM_Gas")

parser.add_argument( "--steps_per_iteration",       type=int, default = 500)
parser.add_argument( "--niterations",               type=int, default = 1000)
parser.add_argument( "--rotations_per_iteration",   type=int, default = 500)

parser.add_argument( "--low_temperature",           type=float, default = 300.)
parser.add_argument( "--high_temperature",          type=float, default = 600.)
parser.add_argument( "--ntemperatures",             type=int, default = 8)

parser.add_argument( "--nc_traj_file",              type=str, default = "traj.nc")

args = parser.parse_args()

def is_md_done(nc_file, niterations):
    if not os.path.isfile(nc_file):
        return False

    if os.path.getsize(nc_file) == 0:
        return False

    nc_handle = netCDF4.Dataset(nc_file)
    if nc_handle.variables["positions"].shape[0] < niterations:
        return False
    
    if type(nc_handle.variables["positions"][-1]) == np.ma.core.MaskedArray:
        return False
    return True

def geometric_progression(low, high, n):
    assert low > 0 and high > 0, "low and high must be positive"
    assert high > low, "high must be higher than low"
    log_scale = np.linspace(np.log(low), np.log(high), n)
    return np.exp(log_scale)

if not is_md_done(args.nc_traj_file, args.niterations):
    temperatures = geometric_progression(args.low_temperature, args.high_temperature, args.ntemperatures)

    tremd = OpenMM_TREMD(args.ligand_prmtop, args.ligand_inpcrd, args.phase, temperatures)
    tremd.run(args.nc_traj_file, args.steps_per_iteration, args.niterations, args.rotations_per_iteration)

else:
    print args.nc_traj_file + " is good, so nothing to be done!"

