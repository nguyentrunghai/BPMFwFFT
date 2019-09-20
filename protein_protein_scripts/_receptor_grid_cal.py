"""
functions that perform grid calculation
"""

from __future__ import print_function

import sys
import os
import numpy as np
import netCDF4

sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from IO import InpcrdLoad
from grids import Grid, RecGrid


def _distance(coord1, coord2):
    assert len(coord1)==len(coord1)==3, "coord must have len 3"
    d = np.array(coord1) - np.array(coord2)
    d = (d**2).sum()
    return np.sqrt(d)


def _max_inter_atom_distance(inpcrd):
    """
    inpcrd: str, name of inpcrd file
    """
    crd = InpcrdLoad(inpcrd).get_coordinates()
    max_d = 0.
    for i in range(crd.shape[0]-1):
        for j in range(i+1, crd.shape[0]):
            d = _distance(crd[i], crd[j])
            if d > max_d:
                max_d = d
    return max_d


def _max_box_edge(inpcrd):
    """
    """
    crd = InpcrdLoad(inpcrd).get_coordinates()
    dx = crd[:,0].max() - crd[:,0].min()
    dy = crd[:,1].max() - crd[:,1].min()
    dz = crd[:,2].max() - crd[:,2].min()
    return max([dx, dy, dz])


def rec_grid_cal(prmtop, lj_scale, rec_inpcrd, lig_inpcrd, 
                spacing, buffer, grid_out, pdb_out, box_out):
    """
    prmtop: str, prmtop file for receptor
    lj_scale:   float, 0 < lj_scale <=1
    rec_inpcrd: str, inpcrd file for receptor
    lig_inpcrd: str, inpcrd file for ligand, used to determine grid size
    spacing:    float
    spacing:    float
    grid_out:   str, name of output nc file
    pdb_out:    str, name of output pdb file
    box_out:    str, name of output box
    """
    #ligand_max_size = _max_inter_atom_distance(lig_inpcrd)
    #print "Ligand maximum inter-atomic distance: %f"%ligand_max_size

    ligand_max_box_edge = _max_box_edge(lig_inpcrd)
    print("Ligand maximum box edge: %f"%ligand_max_box_edge)
    total_buffer = np.ceil(ligand_max_box_edge + buffer)
    print("Total buffer for receptor gird: %f"%total_buffer)

    bsite_file = None
    potential_grid = RecGrid(prmtop, lj_scale,
                                rec_inpcrd, 
                                bsite_file,
                                grid_out,
                                new_calculation=True, 
                                spacing=spacing, buffer=total_buffer)

    potential_grid.write_pdb(pdb_out, "w")
    potential_grid.write_box(box_out)

    return None


def is_nc_grid_good(nc_grid_file):
    if not os.path.exists(nc_grid_file):
        return False

    if os.path.getsize(nc_grid_file) == 0:
        return False

    nc_handle = netCDF4.Dataset(nc_grid_file, "r")
    nc_keys = nc_handle.variables.keys()
    grid_keys = Grid().get_allowed_keys()
    for key in grid_keys:
        if key not in nc_keys:
            return False
    return True


def get_grid_size_from_nc(grid_nc_file):
    nc_handle = netCDF4.Dataset(grid_nc_file, "r")
    return nc_handle.variables["counts"][0]


def get_grid_size_from_lig_rec_crd(rec_inpcrd, lig_inpcrd, buffer):
    ligand_max_box_edge = _max_box_edge(lig_inpcrd)
    buffer_due_to_ligand = np.ceil(ligand_max_box_edge + buffer)
    box_size = _max_box_edge(rec_inpcrd) +  2.0*buffer_due_to_ligand
    return box_size


