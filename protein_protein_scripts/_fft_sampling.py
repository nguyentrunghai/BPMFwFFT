"""
functions to run FFT sampling
"""
from __future__ import print_function

import sys
import os

import netCDF4
import numpy as np

sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from fft_sampling import Sampling

BSITE_FILE = None


def sampling(rec_prmtop, lj_sigma_scal_fact, 
                rec_inpcrd, grid_nc_file,
                lig_prmtop, lig_inpcrd, 
                lig_coor_nc, nr_lig_conf, start_index,
                energy_sample_size_per_ligand,
                output_nc):
    lig_nc_handle = netCDF4.Dataset(lig_coor_nc, "r")
    lig_coord_ensemble = lig_nc_handle.variables["positions"][start_index : start_index + nr_lig_conf]
    lig_nc_handle.close()

    sampler = Sampling(rec_prmtop, lj_sigma_scal_fact, rec_inpcrd,
                        BSITE_FILE, grid_nc_file, lig_prmtop, lig_inpcrd,
                        lig_coord_ensemble,
                        energy_sample_size_per_ligand, 
                        output_nc,
                        temperature=300.)

    sampler.run_sampling()
    print("Sampling Done")
    return None


def is_sampling_nc_good(nc_file, nr_extracted_lig_conf):
    if not os.path.exists(nc_file):
        return False

    try:
        nc_handle = netCDF4.Dataset(nc_file, "r")
    except RuntimeError as e:
        print(nc_file)
        print(e)
        return True
    else:
        pass
    cond1 = nc_handle.variables["lig_positions"][:].shape[0] == nr_extracted_lig_conf
    if not cond1:
        return False

    cond2 = type(nc_handle.variables["lig_positions"][:]) == np.ndarray
    if not cond2:
        return False

    return True


def parse_nr_ligand_confs(submit_file):
    if os.path.exists(submit_file):
        with open(submit_file, "r") as F:
            for line in F:
                if "--nr_lig_conf" in line:
                    nr_confs = line.split("--nr_lig_conf")[-1]
                    nr_confs = nr_confs.split()[0]
                    nr_confs = int(nr_confs)
                    return nr_confs
    return None


