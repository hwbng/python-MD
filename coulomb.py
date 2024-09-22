from constants import k
import numpy as np

def coulomb_force(Atoms,charges,r,r_unit,cutoff):
    '''
    Calculates the Coulomb force acting on each atom.
    Parameters:
    - Atoms (class Atoms): An object representing the atoms.
    - charges (ndarray (natoms,natoms)): Array of pairwise charges (q1*q2) between atoms.
    - r (ndarray (natoms,natoms)): Array of pairwise distances between atoms.
    - r_unit (ndarray (natoms,natoms)): Array of pairwise unit vectors between atoms.
    - cutoff (float): Cutoff distance for the potential.
    Returns:
    - Atoms (class Atoms): An object representing
    '''
    r[r > cutoff] = np.inf
    coul = k*charges/(r**2) #calculate the force magnitude
    coul = coul[:,:,np.newaxis] * r_unit #calculate the force vector

    forces = np.sum(coul,axis=(1))
    for i in range(Atoms.natoms):
        Atoms[i].force += forces[i]
    Atoms.update_arr()

    return Atoms

def coulomb_potential(Atoms,charges,r,cutoff):
    '''
    Calculate pairwise Coulombic interactions between atoms.
    Parameters:
    - Atoms (class Atoms): An object representing the atoms.
    - charges (ndarray (natoms,natoms)): Array of pairwise charges (q1*q2) between atoms.
    - r (ndarray (natoms,natoms)): Array of pairwise distances between atoms.
    - cutoff (float): Cutoff distance for the potential.
    Returns:
    - Atoms (class Atoms): An object representing the atoms.
    '''
    r[r > cutoff] = np.inf
    coul = k*charges/r

    couls = np.sum(coul,axis=(1))
    for i in range(Atoms.natoms):
        Atoms[i].pE += couls[i]
    Atoms.update_arr()

    return Atoms

def pairwise_charges(Atoms):
    '''
    Calculate pairwise charges between atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
    Returns:
        charge_pairs (ndarray (natoms,natoms)): Array of pairwise charges (q1*q2) between atoms.
    '''
    charge_pairs = np.zeros((Atoms.natoms,Atoms.natoms))
    for i in range(Atoms.natoms):
        for j in range(Atoms.natoms):
            if i != j:
                charge_pairs[i,j] = Atoms.charges[i]*Atoms.charges[j]
    return charge_pairs
