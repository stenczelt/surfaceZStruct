#include <Molecule.h>

Molecule::Molecule(
        Atom centralAtom,
        std::vector<Atom> atoms)
:
    mCentralAtom(centralAtom),
    mAtoms{atoms},
    mNumberOfAtoms(0)
{
    // Central atom + the number of atoms
    mNumberOfAtoms = 1 + mAtoms.size();
}

