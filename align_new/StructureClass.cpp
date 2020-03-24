// Mina jafari
// September 2016

#include "StructureClass.h"
#include <iostream>

StructureClass::StructureClass(int inNumOfAtoms, std::string* inAtomicNames,
                               int* inAtomicNumbers, double* inCoordinates)
{
    mNumOfAtoms = inNumOfAtoms;
    for(int i=0; i<mNumOfAtoms; i++)
    {
        mAtomicNames[i] = inAtomicNames[i];
        mAtomicNumbers[i] = inAtomicNumbers[i];
    }
    for(int j=0; j<3*mNumOfAtoms; j++)
    {
        mCoordinates[j] = inCoordinates[j];
    }
}

int StructureClass::getNumOfAtoms() const
{
    return(mNumOfAtoms);
}

std::string* StructureClass::getAtomicNames() const
{
    return(mAtomicNames);
}

int* StructureClass::getAtomicNumbers() const
{
    return(mAtomicNumbers);
}

double* StructureClass::getCoordinates() const
{
    return(mCoordinates);
}
