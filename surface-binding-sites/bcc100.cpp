// Mina Jafari
// 12-08-2015

#include "bcc100.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <vector>
#include <stdio.h>
#include <cmath>
#include <iomanip>

double bcc100::m_DELTA_Z = 2.5;

bool bcc100::setAtoms(const std::vector<std::string> &xyzFile)
{
    bool isSet = false;
    mVector = xyzFile;

    std::string unwanted;
    std::istringstream iss1(xyzFile[xyzFile.size()-1]);
    iss1 >> unwanted >> mNthAtom[0] >> mNthAtom[1] >> mNthAtom[2] >> unwanted;

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

void bcc100::findHollow()
{
    double hollowX = mNthAtom[0] - (mDeltaX/2);
    double hollowY = mNthAtom[1] - (mDeltaX/2);
    double hollowZ = mNthAtom[2] + m_DELTA_Z;

    std::string val1 = std::to_string(hollowX);
    std::string val2 = std::to_string(hollowY);
    std::string val3 = std::to_string(hollowZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc100-hollow.xyz", std::ofstream::out);

    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
}  //findHollow

void bcc100::findAtop()
{
    std::string val1 = std::to_string(mStarMinusOneAtom[0]);
    std::string val2 = std::to_string(mStarMinusOneAtom[1]);
    std::string val3 = std::to_string(mStarMinusOneAtom[2] + m_DELTA_Z);

    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc100-atop.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findAtop

void bcc100::findBridge()
{
    double brgX = mNthMinusOneAtom[0];
//    double brgY = mNthMinusOneAtom[1] - ( std::abs(mStarMinusOneAtom[1] - mNthMinusOneAtom[1]) / 2 );
    double brgY = mNthMinusOneAtom[1] - (mDeltaX/2);
    double brgZ = mNthAtom[2] + m_DELTA_Z;

    std::string val1 = std::to_string(brgX);
    std::string val2 = std::to_string(brgY);
    std::string val3 = std::to_string(brgZ);

    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc100-brg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findBridge
