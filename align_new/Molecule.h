#ifndef MOLECULE_H
#define MOLECULE_H

#include "Atom.h"
#include <vector>

class Molecule
{
    private:
        //Atom mCentralAtom;
        std::vector<Atom> mAtoms;
        int mNumberOfAtoms;

    public:
        Molecule(/*Atom centralAtom,*/ std::vector<Atom> atoms, int numOfAtoms);
        ~Molecule(){};
        inline int numberOfAtoms() const { return mNumberOfAtoms; }
};


#endif
