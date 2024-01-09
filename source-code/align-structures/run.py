"""
import modules
"""

from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import write, read
from ase.lattice.surface import add_adsorbate, fcc111
from ase.optimize import QuasiNewton

"""
set up slab and adsorbates by reading it from input file
"""
slab = fcc111("Pd", size=(4, 3, 2), vacuum=15)
ads = Atoms("C2H3CHO")
add_adsorbate(slab, ads, 2.5, "ontop")

slabAtoms = read("${input}", format="xyz")
slab.set_positions(slabAtoms.get_positions())

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
# calc = EMT()
slab.set_calculator(calc)

dyn = QuasiNewton(slab, trajectory="trajectory.traj")
dyn.run(fmax=0.05)

write("${outputName}", slab)
