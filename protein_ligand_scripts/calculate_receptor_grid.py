"""
to calculate receptor grids
"""
import sys
import os
import argparse

import netCDF4 as nc

# change this 
sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from grids import RecGrid, is_nc_grid_good


parser = argparse.ArgumentParser()
parser.add_argument( "--receptor_prmtop",       type=str, default = "receptor.prmtop")
parser.add_argument( "--receptor_inpcrd",       type=str, default = "receptor.inpcrd")
parser.add_argument( "--lj_scale_factor",       type=float, default = 1.0)
parser.add_argument( "--bsite",                 type=str, default = "measured_binding_site.py")
parser.add_argument( "--spacing",               type=float, default = 0.25)

parser.add_argument( "--grid_nc_out",           type=str, default = "grid.nc")
parser.add_argument( "--pdb_out",               type=str, default = "receptor.pdb")
parser.add_argument( "--enclosing_box_out",     type=str, default = "box.pdb")
args = parser.parse_args()

if not is_nc_grid_good(args.grid_nc_out):

    potential_grid = RecGrid(args.receptor_prmtop, args.lj_scale_factor, 
                            args.receptor_inpcrd, args.bsite,
                            args.grid_nc_out,
                            new_calculation=True,
                            spacing=args.spacing)

    potential_grid.write_pdb(args.pdb_out, "w")
    potential_grid.write_box(args.enclosing_box_out)

else:
    print args.grid_nc_out + " exists and is good. So nothing to be done!"
