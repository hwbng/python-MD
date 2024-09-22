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
        for i in range(3):
            if atom.pos[i] < 0:
                atom.pos[i] = cell[i] - abs(atom.pos[i] % cell[i])
                atom.vel[i] = abs(atom.vel[i])
            elif atom.pos[i] > cell[i]:
                atom.vel[i] = -abs(atom.vel[i])
                atom.pos[i] = cell[i] - (atom.pos[i] % cell[i])
    Atoms.update_arr()
    return Atoms