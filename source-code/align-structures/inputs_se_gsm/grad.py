#!/usr/bin/python

# this is for parallel gsm, so VASP does not read numOfThreads from submit script


import os
import shutil
import sys

from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import read
from ase.lattice.surface import fcc111, add_adsorbate

nargs = len(sys.argv)
argv1 = sys.argv[1]
argv2 = sys.argv[2]

fname = "scratch/structure" + argv1
tmpDir = os.getenv("PBSTMPDIR")
folder = "scratch/vasp" + argv1

# lattice and initial atom setup for unit cell
slab = fcc111("Pd", size=(4, 3, 2), vacuum=15)
add_adsorbate(slab, "C2H3CHO", 2.5, "ontop")

# current position read in
slabatoms = read(fname)
slab.set_positions(slabatoms.get_positions())

mask = [atom.tag > 1 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))
calc = Vasp(
    xc="PBE",
    lreal="Auto",
    kpts=[1, 1, 1],
    ismear=1,
    sigma=0.2,
    algo="fast",
    istart=0,
    npar=2,
    encut=300,
)
slab.set_calculator(calc)

cwd = os.getcwd()

os.chdir(tmpDir)
if not os.path.exists(folder):
    os.makedirs(folder)
os.chdir(folder)
energy = -slab.get_potential_energy()
grads = -slab.get_forces()

f = open("GRAD" + argv1, "w")
f.write(str(energy))
f.write("\n")
f.write(str(grads))
f.write("\n")
f.close()

shutil.copy("GRAD" + argv1, cwd + "/scratch/")
os.chdir(cwd)
