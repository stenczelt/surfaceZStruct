# set up slab using ASE
from ase import Atoms
from ase.build import fcc111

slab_in = fcc111("Pd", size=(4, 3, 2), vacuum=15)
adsorbate_in = Atoms("C2H3CHO")

# set up slab from POSCAR
# slab_POSCAR_file_name = 'scratch/POSCAR'
