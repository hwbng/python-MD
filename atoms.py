import numpy as np

class Atom:
    '''
    Represents a single atom.
    Attributes:
        id (int): The ID of the atom.
        symbol (str): The symbol of the atom.
        mass (float): The mass of the atom.
        pos (numpy.ndarray): The position of the atom.
        vel (numpy.ndarray): The velocity of the atom.
        accel (numpy.ndarray): The acceleration of the atom.
        force (numpy.ndarray): The force acting on the atom.
        kE (float): The kinetic energy of the atom.
        pE (float): The potential energy of the atom.
    '''
    def __init__(self,group,pos,vel):
        self.id = 'undefined'
        self.group = group
        self.symbol = 'undefined'
        self.mass = 'undefined'
        self.charge = 'undefined'
        self.pos = np.array([float(i) for i in pos])
        self.vel = np.array([float(i) for i in vel])
        self.accel = np.array([0.0,0.0,0.0])
        self.force = np.array([0.0,0.0,0.0])
        self.kE = 0.0
        self.pE = 0.0

    def __str__(self):
        return f'{self.id} {self.group} {self.symbol} {self.pos[0]} {self.pos[1]} {self.pos[2]}'

class Atoms:
    '''
    Represents a list of atoms
    Attributes:
        atom_list (list): A list of Atom objects.
        natoms (int): The number of atoms in the collection.
        pos (ndarray): An array of atom positions.
        vel (ndarray): An array of atom velocities.
        accel (ndarray): An array of atom accelerations.
        force (ndarray): An array of atom forces.
        mass (ndarray): An array of atom masses.
        charges (ndarray): An array of atom charges.
        pE (ndarray): An array of atom potential energies.
        kE (ndarray): An array of atom kinetic energies
    '''
    def __init__(self,atom_list=None):
        if atom_list is None:
            atom_list = [] #allows the object to be overwritten with a new list
        self.atom_list = atom_list
        self.natoms = len(atom_list)

        self.pos = np.array([atom.pos for atom in atom_list])
        self.vel = np.array([atom.vel for atom in atom_list])
        self.accel = np.array([atom.accel for atom in atom_list])
        self.force = np.array([atom.force for atom in atom_list])
        self.mass = np.array([atom.mass for atom in atom_list])[:,np.newaxis]
        self.charges = np.array([atom.charge for atom in atom_list])
        self.pE = np.array([atom.pE for atom in atom_list])
        self.kE = np.array([atom.kE for atom in atom_list])

        for ind,atom in enumerate(self.atom_list):
            atom.id = ind+1

    def __str__(self):
        return '\n'.join([str(atom) for atom in self.atom_list])

    def __getitem__(self,index):
        return self.atom_list[index]

    def update_arr(self):
        self.pos = np.array([atom.pos for atom in self.atom_list])
        self.vel = np.array([atom.vel for atom in self.atom_list])
        self.accel = np.array([atom.accel for atom in self.atom_list])
        self.force = np.array([atom.force for atom in self.atom_list])
        self.mass = np.array([atom.mass for atom in self.atom_list])[:,np.newaxis]
        self.charges = np.array([atom.charge for atom in self.atom_list])
        self.pE = np.array([atom.pE for atom in self.atom_list])
        self.kE = np.array([atom.kE for atom in self.atom_list])

    def append(self,atom):
        self.atom_list.append(atom)
        self.natoms += 1
        atom.id = self.natoms
        self.update_arr()

    def get_potential(self):
        return np.sum([atom.pE for atom in self.atom_list])
    def get_kinetic(self):
        return np.sum([atom.kE for atom in self.atom_list])
    def get_total(self):
        return self.get_kinetic() + self.get_potential()

def species(Atoms,group,symbol,mass,charge,**kwargs):
    '''
    Adds property to each group of atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
        group (int): Group number of the atoms to change.
        symbol (str): The symbol of the group.
        mass (float): The mass of the group.
        charge (float): The charge of the group.
        **kwargs: Additional properties to assign to the group.
    Returns:
        Atoms (class Atoms): Atoms with the updated group.
    '''
    for atom in Atoms:
        if atom.group == group:
            atom.symbol = symbol
            atom.mass = mass
            atom.charge = charge
            for property,value in kwargs.items():
                setattr(atom,property,value)
    Atoms.update_arr()

    return Atoms


