import os
def output_xyz(path,atom_list,t,cell):
    '''
    Writes the current state of the system to an xyz file using the extended xyz format
    (see https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz)
    Parameters:
    - path (str): The path to the file.
    - atom_list (class Atoms): An object representing the atoms.
    - t (float): The current time.
    - cell (numpy ndarray): The cell dimensions.
    '''
    outdir = os.path.join(path, 'outdir')
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, 'output.xyz')

    if t == 0.0:
        if os.path.exists(output_file):
            os.remove(output_file)
        with open(output_file, 'x') as file:
            pass

    with open(output_file,'a') as file:
        file.write(f'{atom_list.natoms}\n')
        file.write(f'time = {t} Lattice="{cell[0]} 0 0 0 {cell[1]} 0 0 0 {cell[2]}"\n')
        for atom in atom_list:
            file.write(f'{atom.symbol} {atom.pos[0]} {atom.pos[1]} {atom.pos[2]} {atom.vel[0]} {atom.vel[1]} {atom.vel[2]}\n')

def output_energy(path,atom_list,t):
    '''
    Writes the energy of the system to a file.
    Parameters:
    - atom_list (class Atoms): An object representing the atoms.
    - t (float): The current time.
    '''
    outdir = os.path.join(path, 'outdir')
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, 'energy.txt')


    if t == 0.0:
        if os.path.exists(output_file):
            os.remove(output_file)
        with open(output_file, 'x') as file:
            file.write('Time / ps,Potential Energy / eV,Kinetic Energy / eV,Total Energy / eV\n')
            pass

    with open(output_file,'a') as file:
        potential_energy = atom_list.get_potential()
        kinetic_energy = atom_list.get_kinetic()
        total = atom_list.get_total()
        file.write(f'{t},{potential_energy},{kinetic_energy},{total}\n')


def output_debug(path,atom_list,t):
    '''
    Writes the current state of the system to a verbose file.
    Parameters:
    - atom_list (class Atoms): An object representing the atoms.
    - t (float): The current time.
    '''
    outdir = os.path.join(path, 'outdir')
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, 'debug.txt')


    if t == 0.0:
        if os.path.exists(output_file):
            os.remove(output_file)
        with open(output_file, 'x') as file:
            pass

    with open(output_file,'a') as file:
        file.write(f'{atom_list.natoms}\n')
        file.write(f'time = {t}\n')
        for atom in atom_list:
            file.write(f'{atom.symbol} {atom.pos[0]} {atom.pos[1]} {atom.pos[2]}, force = {atom.force}, velocity = {atom.vel}, accel = {atom.accel}\n')