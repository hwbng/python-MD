from lennard_jones import pairwise_calc,distance_calc
from verlet import integrator
from coulomb import pairwise_charges
import numpy as np
from write import output_xyz,output_energy,output_debug

def run(Atoms,start,end,cell,dt,temp_bath,cutoff,path,coul=True):
    '''
    Main MD loop.
    Args:
        Atoms (class Atoms): Array of atoms.
        start (float): Start time of simulation.
        end (float): End time of simulation.
        cell (ndarray): Simulation cell size.
        dt (float): Timestep.
        temp_bath (float): Temperature of the thermostat.
        cutoff (float): Cutoff distance for the potential.
    '''
    sigma, epsilon = pairwise_calc(Atoms)
    for t in np.arange(start,end,dt):
        r,r_unit = distance_calc(Atoms)
        charges = pairwise_charges(Atoms)

        Atoms = integrator(Atoms,sigma,epsilon,r,r_unit,dt,temp_bath,cell,charges,cutoff,coul=coul)

        output_xyz(path,Atoms,t,cell)
        output_energy(path,Atoms,t)
        output_debug(path,Atoms,t)

        print(f'Time: {t} ps, Potential Energy: {Atoms.get_potential()}, Kinetic Energy: {Atoms.get_kinetic()}, Total Energy: {Atoms.get_total()}')
