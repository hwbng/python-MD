#Velocity Verlet integrator
from bc import bc
from lennard_jones import force, potential
from vel_rescaling import thermostat, kE_calc
from coulomb import coulomb_force, coulomb_potential

def pos_update(Atoms,dt):
    '''
    Update the position of the atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
        dt (float): Timestep.
    Returns:
        Atoms (class Atoms): Atoms with updated positions
    '''
    for i in range(Atoms.natoms):
        Atoms[i].pos += Atoms[i].vel*dt + 0.5*Atoms[i].accel*dt*dt
    Atoms.update_arr()

def vel_update(Atoms,acc_old,dt):
    '''
    Update the velocity of the atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
        acc_old (ndarray): Array of previous accelerations.
        dt (float): Timestep.
    Returns:
        Atoms (class Atoms): Atoms with updated velocities.
    '''
    for i in range(Atoms.natoms):
        Atoms[i].vel += 0.5*(Atoms[i].accel + acc_old[i])*dt
    Atoms.update_arr()

def accel_update(Atoms):
    '''
    Update the acceleration of the atoms.
    Args:
        Atoms (class Atoms): Array of atoms.
    Returns:
        acc_old (ndarray): Array of previous accelerations.
        Atoms (class Atoms): Atoms with updated accelerations.
    '''
    acc_old = Atoms.accel
    for i in range(Atoms.natoms):
        Atoms[i].accel = Atoms[i].force/Atoms[i].mass
    Atoms.update_arr()
    return acc_old, Atoms

def integrator(Atoms,sigma,epsilon,r,r_unit,dt,temp_bath,cell,charges,cutoff,coul):
    '''
    Velocity Verlet integrator.
    Args:
        Atoms (class Atoms): Array of atoms.
        sigma (ndarray): Array of pairwise sigma values.
        epsilon (ndarray): Array of pairwise epsilon values.
        r (ndarray (natoms,natoms)): Array of pairwise distances between atoms.
        r_unit (ndarray (natoms,natoms,3)): Array of pairwise unit vectors between atoms.
        dt (float): Timestep.
        temp_bath (float): Temperature of the thermostat.
        cell (ndarray): Simulation cell size.
        charges (ndarray): Array of pairwise charges.
        cutoff (float): Cutoff distance for the potential.
    Returns:
        Atoms (class Atoms): Atoms with updated positions, velocities, and accelerations.
    '''
    pos_update(Atoms,dt)

    Atoms = force(Atoms,sigma,epsilon,r,r_unit,cutoff)
    Atoms = potential(Atoms,sigma,epsilon,r,cutoff)

    if coul == True:
        Atoms = coulomb_force(Atoms,charges,r,r_unit,cutoff)
        Atoms = coulomb_potential(Atoms,charges,r,cutoff)

    acc_old, Atoms = accel_update(Atoms)
    vel_update(Atoms,acc_old,dt)
    if temp_bath > 0:
        Atoms = thermostat(Atoms,temp_bath)
    else:
        kE_calc(Atoms)
    Atoms = bc(Atoms,cell)


    return Atoms
