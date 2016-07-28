// Mina Jafari
// 12-08-2015

#include "bcc111.h"
#include "Termination111.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <vector>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <cmath>

void bcc111::findShallowBridge(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = offsetX*mDeltaX + (mNthAtom[0]-mStarAtom[0])/2;
    double offY = offsetY*mDistance + mDistance/6;
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
    ofs.open ("bcc111-shallowBrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findShallowBridge

void bcc111::findDeepBridge(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = offsetX*mDeltaX + (mStarAtom[0]-mNthMinusOneAtom[0])/2;
    double offY = offsetY*mDistance + mDistance/2; 
    double brgX = mStarAtom[0] - offX;
    double brgY = mNthMinusOneAtom[1] - offY;
    double brgZ = mNthAtom[2] + m_DELTA_Z;

    std::string val1 = std::to_string(brgX);
    std::string val2 = std::to_string(brgY);
    std::string val3 = std::to_string(brgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc111-deepBrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findDeepBridge
