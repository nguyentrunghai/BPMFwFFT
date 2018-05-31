"""
run MD, replica exchange (parallel tempering) and calcaulte energy using OpenMM
"""

from sys import stdout

import copy

import numpy as np
import netCDF4

import simtk.openmm
import simtk.openmm.app
import simtk.unit

from rotation import random_rotation

openmm_solvent_models = {  "OpenMM_Gas":None,
                            "OpenMM_GBn":simtk.openmm.app.GBn,
                            "OpenMM_GBn2":simtk.openmm.app.GBn2,
                            "OpenMM_HCT":simtk.openmm.app.HCT,
                            "OpenMM_OBC1":simtk.openmm.app.OBC1,
                            "OpenMM_OBC2":simtk.openmm.app.OBC2 }


KB = 0.001987204134799235  # kcal/mol/K


class OpenMM_MD(object):
    def __init__(self, prmtop, inpcrd, phase="OpenMM_Gas", temperature=300.):
        """
        """
        self._prmtop = simtk.openmm.app.AmberPrmtopFile(prmtop)
        inpcrd = simtk.openmm.app.AmberInpcrdFile(inpcrd)
        
        selected_solvent = openmm_solvent_models[phase]
        system = self._prmtop.createSystem(nonbondedMethod = simtk.openmm.app.NoCutoff, \
                constraints = simtk.openmm.app.HBonds, implicitSolvent = selected_solvent)
        integrator = simtk.openmm.LangevinIntegrator(temperature*simtk.unit.kelvin, 1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
        
        self._simulation = simtk.openmm.app.Simulation(self._prmtop.topology, system, integrator)
        self._simulation.context.setPositions(inpcrd.positions)
        print("Energy minimizing")
        self._simulation.minimizeEnergy()

    def _initialize_nc(self, nc_file_name, niterations):
        nc_handle = netCDF4.Dataset(nc_file_name, mode="w", format="NETCDF4")
        natoms = len( list( self._prmtop.topology.atoms() ) )

        nc_handle.createDimension("three", 3)
        nc_handle.createDimension("natoms", natoms)
        nc_handle.createDimension("nconfs", niterations)

        nc_handle.createVariable("positions", "f8", tuple(["nconfs", "natoms", "three"]))
        return nc_handle
    
    def run(self, nc_file_name, steps_per_iteration=500, niterations=1000):
        """
        """
        self._simulation.reporters.append(simtk.openmm.app.StateDataReporter(stdout, steps_per_iteration,
                                            step=True, potentialEnergy=True))
        
        nc_handle = self._initialize_nc(nc_file_name, niterations)

        for iteration in range(niterations):
            self._simulation.step(steps_per_iteration)
            state = self._simulation.context.getState(getPositions=True, getEnergy=True)
            
            positions = copy.deepcopy(state.getPositions())
            conf = np.array(positions.value_in_unit(simtk.unit.angstrom), dtype=float)    
            conf = random_rotation(conf)

            nc_handle.variables["positions"][iteration,:,:] = conf
        nc_handle.close()

        return None


class OpenMM_TREMD(object):
    def __init__(self, prmtop, inpcrd, phase, temperatures):
        """
        :param prmtop: str, name of AMBER prmtop file
        :param inpcrd: str, name of AMBER coordinate file
        :param phase: str
        :param temperatures: list or ndarray of float
        """
        self._temperatures = list(temperatures)
        self._simulations  = self._create_simulations(prmtop, inpcrd, phase, temperatures)
        self._acepted_exchange = np.zeros([len(self._temperatures)-1], dtype=float)
        self._nr_attempts = np.zeros([len(self._temperatures)-1], dtype=float)
        self._start = 0

    def run(self, nc_file_name, steps_per_iteration, niterations, rotations_per_iteration):
        """
        """
        nc_handle = self._initialize_nc(nc_file_name, niterations, rotations_per_iteration)

        nrotations = 0
        for iteration in range(niterations):

            self._md_evolve(steps_per_iteration)

            energies = self._get_potential_energies()
            nc_handle.variables["energies"][iteration, :] = np.array(energies, dtype=float)

            positions = self._get_positions()
            for state, p in enumerate(positions):
                crd = np.array(p.value_in_unit(simtk.unit.angstrom), dtype=float)
                nc_handle.variables["positions"][iteration, state, :, :] = crd

            crd = np.array(positions[0].value_in_unit(simtk.unit.angstrom), dtype=float)
            for rotation in range(rotations_per_iteration):
                rotated_crd = random_rotation(crd)
                nc_handle.variables["rotated_positions"][nrotations, :, :] = rotated_crd
                nrotations += 1

            self._exchange(self._start)
            self._switch_start()

            print("iteration ", iteration)
            print("energies ", energies)

        acceptance_rate = self._acepted_exchange / self._nr_attempts
        nc_handle.variables["acceptance_rate"][:] = acceptance_rate
        nc_handle.close()
        return None

    def _switch_start(self):
        if self._start == 0:
            self._start = 1
        elif self._start == 1:
            self._start = 0
        else:
            raise RuntimeError("self._start is %d"%self._start)
        return

    def _create_simulation(self, prmtop, inpcrd, phase, temperature):
        prmtop = simtk.openmm.app.AmberPrmtopFile(prmtop)
        inpcrd = simtk.openmm.app.AmberInpcrdFile(inpcrd)
        selected_solvent = openmm_solvent_models[phase]

        system = prmtop.createSystem(nonbondedMethod=simtk.openmm.app.NoCutoff, constraints=simtk.openmm.app.HBonds, implicitSolvent=selected_solvent)
        integrator = simtk.openmm.LangevinIntegrator(temperature*simtk.unit.kelvin, 1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)

        simulation = simtk.openmm.app.Simulation(prmtop.topology, system, integrator)
        simulation.context.setPositions(inpcrd.positions)
        simulation.minimizeEnergy()
        return simulation

    def _create_simulations(self, prmtop, inpcrd, phase, temperatures):
        simulations = []
        for temperature in temperatures:
            sim = self._create_simulation(prmtop, inpcrd, phase, temperature)
            simulations.append(sim)
        return simulations

    def _md_evolve(self, steps):
        for sim in self._simulations:
            sim.step(steps)
        return None

    def _get_potential_energies(self):
        energies = []
        for sim in self._simulations:
            state = sim.context.getState(getEnergy=True)
            energy = copy.deepcopy(state.getPotentialEnergy())
            e = energy.value_in_unit(simtk.unit.kilocalorie_per_mole)
            energies.append(e)
        return energies

    def _get_positions(self):
        positions = []
        for sim in self._simulations:
            state = sim.context.getState(getPositions=True)
            pos = copy.deepcopy(state.getPositions())
            positions.append(pos)
        return positions

    def _get_velocities(self):
        velocities = []
        for sim in self._simulations:
            state = sim.context.getState(getVelocities=True)
            vel = copy.deepcopy(state.getVelocities())
            velocities.append(vel)
        return velocities

    def _metropolis_prob(self, T1, E1, T2, E2):
        assert T1 > 0 and T2 > 0, "T1 and T2 must be possitive"
        deltaE = ((1./T1 - 1./T2) / KB) * (E2 - E1)
        return np.exp(-deltaE)

    def _exchange(self, start):
        assert start in [0, 1], "start must be either 0 or 1"
        simulation_pairs = zip( self._simulations[start : -1 : 2], self._simulations[start+1 : : 2] )

        energies = self._get_potential_energies()
        energy_pairs = zip( energies[start : -1 : 2], energies[start+1 : : 2] )

        temperature_pairs = zip( self._temperatures[start : -1 : 2], self._temperatures[start+1 : : 2] )

        positions = self._get_positions()
        position_pairs = zip( positions[start : -1 : 2], positions[start+1 : : 2] )

        velocities = self._get_velocities()
        velocity_pairs = zip( velocities[start : -1 : 2], velocities[start+1 : : 2] )

        pair_indices = range(start, len(self._temperatures)-1, 2) 

        for i in range(len(simulation_pairs)):
            E1 = energy_pairs[i][0]
            E2 = energy_pairs[i][1]
            T1 = temperature_pairs[i][0]
            T2 = temperature_pairs[i][1]

            x_prob = self._metropolis_prob(T1, E1, T2, E2)
            if x_prob >= 1. or x_prob > np.random.random():
                simulation_pairs[i][0].context.setPositions( position_pairs[i][1] )
                simulation_pairs[i][1].context.setPositions( position_pairs[i][0] )

                vel_unit = simtk.unit.nanometer/simtk.unit.picosecond

                v1 = np.sqrt(T1 / T2 ) * np.array( velocity_pairs[i][1].value_in_unit(vel_unit) )
                v2 = np.sqrt(T2 / T1 ) * np.array( velocity_pairs[i][0].value_in_unit(vel_unit) )

                simulation_pairs[i][0].context.setVelocities( simtk.unit.quantity.Quantity(v1, unit=vel_unit) )
                simulation_pairs[i][1].context.setVelocities( simtk.unit.quantity.Quantity(v2, unit=vel_unit) )

                self._acepted_exchange[pair_indices[i]] += 1. 

            self._nr_attempts[pair_indices[i]] += 1.
        return 

    def _initialize_nc(self, nc_file_name, niterations, rotations_per_iteration):
        """
        """
        nc_handle = netCDF4.Dataset(nc_file_name, mode="w", format="NETCDF4")
        natoms = len( list( self._simulations[0].topology.atoms() ) )

        nc_handle.createDimension("three", 3)
        nc_handle.createDimension("natoms", natoms)
        nc_handle.createDimension("nstates", len(self._simulations) )
        nc_handle.createDimension("niterations", niterations)
        nc_handle.createDimension("nrotations", niterations * rotations_per_iteration)
        nc_handle.createDimension("npairs", len(self._simulations) - 1)

        nc_handle.createVariable("energies", "f8", tuple(["niterations", "nstates"]) )
        nc_handle.createVariable("positions", "f4", tuple(["niterations", "nstates", "natoms", "three"]))
        nc_handle.createVariable("rotated_positions", "f4", tuple(["nrotations", "natoms", "three"]))
        nc_handle.createVariable("acceptance_rate", "f8", tuple(["npairs"]))
        nc_handle.createVariable("temperatures", "f8", tuple(["nstates"]))

        nc_handle.variables["temperatures"][:] = np.array(self._temperatures)
        return nc_handle


def openmm_energy(prmtop_file, crd, phase):
    """
    crd in angstroms
    crd can be an inpcrd file name, a array with shape (natom, 3) or shape (nframe, natom, 3)
    """
    implicit_solvents = openmm_solvent_models
    selected_solvent = implicit_solvents[phase]
    
    if type(crd) == str:
        inpcrd = simtk.openmm.app.AmberInpcrdFile(crd)
        position = inpcrd.positions
        crd = position.value_in_unit(simtk.unit.angstrom) 
    
    crd_ensemble = np.array(crd, dtype=float)
    if len(crd_ensemble.shape) == 2:
        crd_ensemble = np.array([crd_ensemble], dtype=float)
    
    if len(crd_ensemble.shape) != 3:
        raise RuntimeError("crd_ensemble has wrong shape")
    
    prmtop = simtk.openmm.app.AmberPrmtopFile(prmtop_file)
    system = prmtop.createSystem(nonbondedMethod=simtk.openmm.app.NoCutoff,
                                 constraints=None, implicitSolvent=selected_solvent)
    dummy_integrator = simtk.openmm.VerletIntegrator(0.002*simtk.unit.picoseconds)
    simulation = simtk.openmm.app.Simulation(prmtop.topology, system, dummy_integrator)
    
    pot_energies = []
    for conf in crd_ensemble:
        simulation.context.setPositions(conf / 10.)  # convert to nano-meter which is internal unit of OpenMM
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        pot_energies.append(energy.value_in_unit(simtk.unit.kilocalorie_per_mole))
    
    return np.array(pot_energies, dtype=float)

