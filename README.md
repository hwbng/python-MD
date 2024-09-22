# Python MD
A basic molecular dynamics program written in python, mainly to practice object oriented programming (includes an atom/atoms object similar to ASE's), the underlying theory behind MD and vectorisation using NumPy broadcasting (to calculate pairwise distances, etc.)

### Features
- Lennard-Jones and coulomb interactions
- Reflective boundary conditions
- Velocity rescaling
- output .xyz file (trajectory can be read with ASE, OVITO, etc.)

### Further implementations that could be added
- Periodic boundary conditions
- Long range forces
- Input file reader
- A proper thermostat
- Bond potentials & angles

## Examples
#### Lennard-Jones forces acting on oxygen atoms (no velocity rescaling)
![](https://github.com/hwbng/python-MD/blob/main/gifs/o2_lj.gif)

![](https://github.com/hwbng/python-MD/blob/main/gifs/lj_o32.gif)
