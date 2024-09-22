import numpy as np
from constants import kb

def kE_calc(Atoms):
    '''
    Calculates the kinetic energy of the system.
    Parameters:
    - Atoms (class Atoms): An object representing the atoms.
    Returns:
    - Atoms (class Atoms): An object representing the atoms.
    '''
    atom_kE = 0.5*Atoms.mass*Atoms.vel**2
    for i in range(Atoms.natoms):
        Atoms[i].kE = atom_kE[i]
    Atoms.update_arr()
    return np.sum(atom_kE)

def thermostat(Atoms,temp_bath):
    '''
    Rescales the velocities of the atoms to match the desired temperature.
    Parameters:
    - Atoms (class Atoms): An object representing the atoms.
    - temp_bath (float): The desired temperature.
    Returns:
    - Atoms (class Atoms): An object representing the atoms.
    '''
    kE = kE_calc(Atoms)
    temp_sys = (2*kE)/(3*kb*Atoms.natoms) #system temperature using stat. mech
    scale = np.sqrt(temp_bath/temp_sys)
    for i in range(Atoms.natoms):
        Atoms[i].vel *= scale
    Atoms.update_arr()
    kE = kE_calc(Atoms)
    return Atoms

