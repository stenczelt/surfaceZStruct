#include "Molecule.h"

Molecule::Molecule(
        /*Atom centralAtom,*/
        std::vector<Atom> atoms, int numOfAtoms)
:
    /*mCentralAtom(centralAtom),*/
    mAtoms(atoms),
    mNumberOfAtoms(numOfAtoms)
{
    // Central atom + the number of atoms
    mNumberOfAtoms = 1 + mAtoms.size();
}

