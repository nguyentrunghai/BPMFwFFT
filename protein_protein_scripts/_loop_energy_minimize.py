"""
"""
from __future__ import print_function

import os

import numpy as np
import simtk.openmm
import simtk.openmm.app
import simtk.unit
import mdtraj

REC_MIN_LIST = "receptor_minimize_list.dat"
LIG_MIN_LIST = "ligand_minimize_list.dat"

COMPLEX_INPCRD = "complex.inpcrd"
COMPLEX_PRMTOP = "complex.prmtop"

REC_INPCRD = "receptor.inpcrd"
LIG_INPCRD = "ligand.inpcrd"


def _read_res_to_minimize(file):
    """
    read which residue to minimize
    :param file: str
    :return:
    """
    residues_min = []
    with open(file, "r") as F:
        for line in F:
            if "#" in line:
                nresidues = int(line.split()[2])
            else:
                residues_min.append(int(line))
    return residues_min, nresidues


def _combine_res_to_minimize_rec_lig(rec_res_to_minimize, lig_res_to_minimize):
    """
    combine residues to minimize in both receptor and ligand
    :param rec_res_to_minimize:
    :param lig_res_to_minimize:
    :return:
    """
    rec_residues_min, rec_nresidues = _read_res_to_minimize(rec_res_to_minimize)

    lig_residues_min, lig_nresidues = _read_res_to_minimize(lig_res_to_minimize)
    lig_residues_min = [i+rec_nresidues for i in lig_residues_min]

    all_res_to_minimize = rec_residues_min + lig_residues_min

    return all_res_to_minimize


def openMM_minimize(in_dir, out_dir, out_pdb="min.pdb", minimize_all=False):
    """
    if minimize_all, the whole complex will be minimize
    res_to_minimize:    list of residue indices (starting 1) allowed to move

    :param in_dir:
    :param out_dir:
    :param out_pdb:
    :param minimize_all:
    :return:
    """
    rec_res_to_minimize = os.path.join(in_dir, REC_MIN_LIST)
    lig_res_to_minimize = os.path.join(in_dir, LIG_MIN_LIST)
    res_to_minimize = _combine_res_to_minimize_rec_lig(rec_res_to_minimize, lig_res_to_minimize)
    res_to_minimize = [i-1 for i in res_to_minimize]    # make it start 0

    prmtop = simtk.openmm.app.AmberPrmtopFile( os.path.join(in_dir, COMPLEX_PRMTOP) )
    inpcrd = simtk.openmm.app.AmberInpcrdFile( os.path.join(in_dir, COMPLEX_INPCRD) )
    system = prmtop.createSystem(nonbondedMethod = simtk.openmm.app.NoCutoff, 
            constraints = None, implicitSolvent = None)

    atoms = prmtop.topology.atoms()
    for i, atom in enumerate(atoms):
        if atom.residue.index not in res_to_minimize:
            system.setParticleMass(i, 0*simtk.unit.dalton)

    integrator = simtk.openmm.VerletIntegrator(0.001*simtk.unit.picoseconds)
    simulation = simtk.openmm.app.Simulation( prmtop.topology, system, integrator )
    simulation.context.setPositions( inpcrd.positions )
    #simulation.reporters.append(mdtraj.reporters.DCDReporter(out_prefix + ".dcd", 100))

    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    print("Energy before loop minimization ", energy.value_in_unit(simtk.unit.kilocalorie_per_mole))

    print("Loop minimizing ...")
    simulation.minimizeEnergy()

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy()
    print("Energy after loop minimization ", energy.value_in_unit(simtk.unit.kilocalorie_per_mole))

    #----
    print("Whole minimizing ...")
    positions = state.getPositions()
    prmtop = simtk.openmm.app.AmberPrmtopFile( os.path.join(in_dir, COMPLEX_PRMTOP) )
    system = prmtop.createSystem(nonbondedMethod = simtk.openmm.app.NoCutoff,
            constraints = None, implicitSolvent = None)
    integrator = simtk.openmm.VerletIntegrator(0.001*simtk.unit.picoseconds)
    simulation = simtk.openmm.app.Simulation( prmtop.topology, system, integrator )
    simulation.context.setPositions( positions )
    simulation.minimizeEnergy(tolerance=1.*simtk.unit.kilojoule_per_mole, maxIterations=2000)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy()
    print("Energy after whole minimization ", energy.value_in_unit(simtk.unit.kilocalorie_per_mole))

    positions = state.getPositions()
    crd = np.array( positions.value_in_unit(simtk.unit.angstrom), dtype=float )

    simtk.openmm.app.PDBFile.writeFile(simulation.topology, state.getPositions(), 
            open( os.path.join(out_dir, out_pdb), 'w') )

    rec_natoms = _read_natoms( os.path.join(in_dir, REC_INPCRD) )
    rec_crd = crd[:rec_natoms]
    lig_crd = crd[rec_natoms:]

    _write_inpcrd(rec_crd, os.path.join(out_dir, REC_INPCRD))
    _write_inpcrd(lig_crd, os.path.join(out_dir, LIG_INPCRD))

    print("Minimizing done!")
    return None


def _read_natoms(inpcrd):
    with open(inpcrd, "r") as F:
        F.readline()
        natoms = int(F.readline())
    return natoms


def _write_inpcrd(crd, file):
    out_str = """default_name\n"""
    out_str += """%6d\n"""%len(crd)
    for i in range(len(crd)):
        out_str += """%12.7f%12.7f%12.7f"""%( crd[i][0], crd[i][1], crd[i][2] )
        if (i+1)%2 == 0:
            out_str += """\n"""
    open( file, "w").write( out_str )
    return None


