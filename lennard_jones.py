#Calculates the Lennard-Jones potential and force between atoms
import numpy as np

def pairwise_calc(Atoms):
    '''
    Calculate pairwise sigma/epsilon values.
    Sigma is calculated as the average between two atoms sigma values.
    Epsilon is calculated as the geometric mean between two atoms epsilon values.
    Args:
        Atoms (class Atoms): Array of atoms.
    Returns:
        sigma (ndarray (natoms,natoms)): Array of pairwise sigma values.
        epsilon (ndarray (natoms,natoms)): Array of pairwise epsilon values.
    '''

    sigma = np.zeros((Atoms.natoms,Atoms.natoms))
    epsilon = np.zeros((Atoms.natoms,Atoms.natoms))

    for i in range(Atoms.natoms):
        for j in range(Atoms.natoms):
            if i != j:
                if Atoms[i].group != Atoms[j].group:
                    #if ateoms are in different groups, calculate averages
                    sigma[i,j] = 0.5*(Atoms[i].sigma + Atoms[j].sigma)
                    epsilon[i,j] = np.sqrt(Atoms[i].epsilon*Atoms[j].epsilon)
                else:
                    #else the values are equal
                    sigma[i,j] = Atoms[j].sigma
                    epsilon[i,j] = Atoms[j].epsilon
    return sigma, epsilon

def distance_calc(Atoms):
    '''
    Calculate pairwise distances between atoms. Uses numpy broadcasting to vectorise
    the calculation.
    Atoms.pos[:, np.newaxis, :] creates a 3D array of shape (natoms, 1, natoms)
    and Atoms.pos[np.newaxis, :, :] creates a 3D array of shape (1, natoms, natoms).
    Subtraction of these two arrays creates a 3D array of shape (natoms, natoms, 3),
    which contains the xyz difference between atoms i and j in a matrix H_ij.

    Args:
        Atoms (class Atoms): Array of atoms.
    Returns:
        r (ndarray (natoms,natoms)): Array of pairwise distances between atoms.
        r_unit (ndarray (natoms,natoms)): Array of pairwise unit vectors between atoms.
    '''
    diff = Atoms.pos[:, np.newaxis, :] - Atoms.pos[np.newaxis, :, :]
    #calculate distance using the xyz values of the difference
    r = np.linalg.norm(diff,axis=(2))
    #avoid division by zero for distance of 0
    r[r == 0] = np.inf
    r_unit = diff/r[:,:,np.newaxis]

    return r,r_unit


def potential(Atoms,sigma,epsilon,r,cutoff):
    '''
    Calculate pairwise interactions between atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
        sigma (ndarray): Array of pairwise sigma values.
        epsilon (ndarray): Array of pairwise epsilon values.
        r (ndarray (natoms,natoms)): Array of pairwise distances between atoms.
    Returns:
        Atoms (class Atoms): Atoms with updated potential energies.
    '''
    r[r > cutoff] = np.inf #set distances greater than cutoff to infinity
    pot = 4*epsilon*( (sigma/r)**12 - (sigma/r)**6 ) #LJ formula, calculates potential in xyz directions
    potentials = np.sum(pot,axis=(1)) #total potential energy for each atom
    for i in range(Atoms.natoms):
        Atoms[i].pE = potentials[i]
    Atoms.update_arr()

    return Atoms

def force(Atoms,sigma,epsilon,r,r_unit,cutoff):
    '''
    Calculate pairwise forces between atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
        sigma (ndarray): Array of pairwise sigma values.
        epsilon (ndarray (natoms,natoms)): Array of pairwise epsilon values.
        r (ndarray (natoms,natoms)): Array of pairwise distances between atoms.
    Returns:
        Atoms (class Atoms): Atoms with updated forces.
    '''
    r[r > cutoff] = np.inf #set distances greater than cutoff to infinity
    f = (24*epsilon/r)*(2*(sigma/r)**12 - (sigma/r)**6) #dLJ formula for magnitude of force
    f = f[:,:,np.newaxis] * r_unit #multiply by unit vector to get force in xyz directions
    forces = np.sum(f,axis=1) #each atom has natoms forces acting on it, sum them to get total force
    # (e.g. atom 1 is affected by atom2, atom3, etc. so we sum all forces acting on atom1)
    for i in range(Atoms.natoms):
        Atoms[i].force = forces[i]
    Atoms.update_arr()
    return Atoms

