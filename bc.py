import numpy as np
from atoms import Atoms, Atom
#boundary conditions: reflect off cell walls
#causes atoms to get stuck in corners/edges of the cell
def bc(Atoms,cell):
    '''
    Implements reflective boundary conditions
    Parameters:
    - Atoms (class Atoms): Input system
    - cell (numpy ndarray): Simulation cell size
    Returns:
    - Atoms (class Atoms): The atomic positions after applying BC.
    '''
    for atom in Atoms:
        #mask is True if atom is outside the cell
        #use mask instead of atom.pos < 0 in np.where to prevent error with
        #updating position before updating velocity
        mask = atom.pos < 0
        #if mask is True (i.e. outside cell), velocity is reversed if travelling 'out'
        atom.vel = np.where(mask, abs(atom.vel), atom.vel)
        atom.pos = np.where(mask, cell - abs(atom.pos % cell), atom.pos)

        #mask is True if atom is outside the cell
        mask = atom.pos > cell
        #invert the direction if travelling outside the cell
        atom.vel = np.where(mask, -abs(atom.vel), atom.vel)
        atom.pos = np.where(mask, cell - (atom.pos % cell), atom.pos)

    Atoms.update_arr()
    return Atoms