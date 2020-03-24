#ifndef ATOM_H
#define ATOM_H

#include <string>
#include "Coordinates.h"

class Atom
{
    private:
        std::string mName;
        //unsigned int mAtomicNumber;
        Coordinates mCoordinates;

    public:
        Atom(std::string name, /*unsigned int atomicNumber, */double x, double y, double z);
        ~Atom() {}

        inline const std::string name() const { return mName; }
        //inline const int atomicNumber() const { return mAtomicNumber; }

        const Coordinates& coordinates() const { return mCoordinates; }

};


#endif
