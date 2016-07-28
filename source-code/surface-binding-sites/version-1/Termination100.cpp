// Mina Jafari
// 12-09-2015

#include "Termination100.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <vector>
#include <stdio.h>
#include <cmath>
#include <iomanip>

double Termination100::m_DELTA_Z = 2.5;

bool Termination100::setAtoms(const std::vector<std::string> &xyzFile)
{
    bool isSet = false;
    mVector = xyzFile;

    // setting mNthAtom 
    std::string unwanted;
    std::istringstream iss1(xyzFile[xyzFile.size()-1]);
    iss1 >> unwanted >> mNthAtom[0] >> mNthAtom[1] >> mNthAtom[2] >> unwanted;

    // setting mNthMinusOneAtom
    std::istringstream iss2(xyzFile[xyzFile.size()-2]);
    iss2 >> unwanted >> mNthMinusOneAtom[0] >> mNthMinusOneAtom[1] >> mNthMinusOneAtom[2] >> unwanted;

    mDeltaX = mNthAtom[0] - mNthMinusOneAtom[0];
    mDeltaY = mNthAtom[1] - mNthMinusOneAtom[1];
    
    //setting mStarAtom & mStarMinusOneAtom
    // n-1 to n vector intersects w/ Y axis
    if (std::abs(mDeltaX) > std::abs(mDeltaY))
    {
        mStarAtom[0] = mNthAtom[0]; // x component
        mStarAtom[1] = mNthAtom[1] - mDeltaX; // y component
        mStarAtom[2] = mNthAtom[2]; // z component

        mStarMinusOneAtom[0] = mNthMinusOneAtom[0]; // x component
        mStarMinusOneAtom[1] = mNthMinusOneAtom[1] - mDeltaX; // y component
        mStarMinusOneAtom[2] = mNthMinusOneAtom[2]; // z component
        isSet = true;
    }   

    return (isSet);
}

void Termination100::findHollow(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = offsetX * mDeltaX + (mDeltaX/2);
    double offY = offsetY * mDeltaX + (mDeltaX/2);
    double hollowX = mNthAtom[0] - offX;
    double hollowY = mNthAtom[1] - offY;
    double hollowZ = mNthAtom[2] + m_DELTA_Z;
    if (hollowX < 0 || hollowY < 0)
    {
        std::cout << "ERROR: offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(hollowX);
    std::string val2 = std::to_string(hollowY);
    std::string val3 = std::to_string(hollowZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("Termination100-hollow.xyz", std::ofstream::out);

    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
}  //findHollow

void Termination100::findAtop(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = offsetX * mDeltaX;
    double offY = offsetY * mDeltaX;
    double atopX = mNthAtom[0] - offX;
    double atopY = mNthAtom[1] - offY;
    double atopZ = mNthAtom[2] + m_DELTA_Z;
    if (atopX < 0 || atopY < 0)
    {
        std::cout << "ERROR: offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(atopX);
    std::string val2 = std::to_string(atopY);
    std::string val3 = std::to_string(atopZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("Termination100-atop.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findAtop

void Termination100::findBridge(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = offsetX * mDeltaX + (mDeltaX/2);
    double offY = offsetY * mDeltaX;
    double brgX = mNthAtom[0] - offX;
    double brgY = mNthAtom[1] - offY;
    double brgZ = mNthAtom[2] + m_DELTA_Z;
    if (brgX < 0 || brgY < 0)
    {
        std::cout << "ERROR: offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(brgX);
    std::string val2 = std::to_string(brgY);
    std::string val3 = std::to_string(brgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("Termination100-brg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findBridge
