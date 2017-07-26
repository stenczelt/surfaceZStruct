#include "Atom.h"

Atom::Atom(
        std::string name, 
        //unsigned int atomicNumber,
        double x,
        double y,
        double z)
:
    mName(name),
    //mAtomicNumber(atomicNumber),
    mCoordinates(x, y, z)
{
}

