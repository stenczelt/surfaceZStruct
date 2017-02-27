// Mina Jafari, September 2016
#ifndef STRUCTURECLASS_H
#define STRUCTURECLASS_H

//#include <vector>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
#include <stdio.h>
//#include <cmath>
//#include <cstdlib>
#include <string>


class StructureClass
{
    private:
        int mNumOfAtoms = 0;
        // Atomic symbols of molecule
        std::string* mAtomicNames;
        // Atomic numbers of molecule
        int* mAtomicNumbers;
        // Cartesian coordinates of the fragments
        double* mCoordinates;

    public:
        //void setNumOfAtoms(int inNumOfAtoms);
        //void setAtomicNames(std:string* inAtomicNames);
        //void setAtomicNumbers(int* inAtomicNumbers);
        //void setAtomicCoordinates(double* inAtomicCoordinates);
        StructureClass(int inNumOfAtoms, std::string* inAtomicNames,
                       int* inAtomicNumbers, double* inCoordinates);
        int getNumOfAtoms() const;
        std::string* getAtomicNames() const;
        int* getAtomicNumbers() const;
        double* getCoordinates() const;
};

#endif

