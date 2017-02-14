#!/usr/bin/python
import sys 
sys.path.append('/export/zimmerman/paulzim/ase')
# this is for parallel gsm, so VASP does not read numOfThreads from submit script
import os
os.environ['OMP_NUM_THREADS'] = '1'

from ase.lattice.surface import surface
from ase import Atoms,Atom
from ase.visualize import view
from ase.calculators.emt import EMT
from ase.calculators.nwchem import NWChem
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import write
from ase.io import read
from ase.lattice.surface import fcc111,add_adsorbate,fcc100,fcc110,hcp0001,bcc111
import os
import sys
import shutil

nargs = len(sys.argv)
argv1 = sys.argv[1]
argv2 = sys.argv[2]

fname = 'scratch/structure'+argv1
tmpDir = os.getenv('PBSTMPDIR')
folder = 'scratch/vasp'+argv1

#lattice and initial atom setup for unit cell
slab = fcc100('Cu', size=(4,3,2), vacuum=10)
add_adsorbate(slab, 'NH3', 2.5, 'ontop')

#current position read in
slabatoms = read(fname)
slab.set_positions(slabatoms.get_positions())

mask = [atom.tag > 1 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))
calc = Vasp(xc='PBE', lreal='Auto',kpts=[1,1,1],ismear=1,sigma=0.2,algo='fast',istart=0,npar=2,encut=300)
slab.set_calculator(calc)

cwd = os.getcwd()

os.chdir(tmpDir)
if not os.path.exists(folder):
    os.makedirs(folder)
os.chdir(folder)
energy = - slab.get_potential_energy()
grads = - slab.get_forces()

f = open('GRAD'+argv1, 'w')
f.write(str(energy))
f.write('\n')
f.write(str(grads))
f.write('\n')
f.close()

shutil.copy('GRAD'+argv1, cwd+'/scratch/')
os.chdir(cwd)


