import numpy as np
import matplotlib.pyplot as plt
from atoms import Atom,Atoms,species
import numpy as np
from MD import run


cell = np.array([20,20,20])
atoms = Atoms()
oxy = Atom('O',np.array([12.5,10,10]),np.array([0,0,0]))
atoms.append(oxy)
oxy = Atom('O',np.array([7.5,10,10]),np.array([0,0,0]))
atoms.append(oxy)
#Values from https://www.sciencedirect.com/science/article/pii/S0378381207005912, converted to eV and Angstrom
species(Atoms=atoms,group='O',symbol='O',mass=16.0,charge=-0.8476,epsilon=0.006516018,sigma=3.1169)
species(Atoms=atoms,group='H',symbol='H',mass=1.0,charge=0.4238,sigma=0.98,epsilon=0.0003358024)

run(atoms,start=0,end=2000,dt=0.005,cell=cell,temp_bath=0,cutoff=30,path='/home/weibo/Documents/MD/test',coul=False)