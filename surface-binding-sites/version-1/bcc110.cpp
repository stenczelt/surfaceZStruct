// Mina Jafari
// 12-08-2015

#include "bcc110.h"
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

double bcc110::m_DELTA_Z = 2.5;

bool bcc110::setAtoms(const std::vector<std::string> &xyzFile)
{
    bool isSet = true;
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

    mDistance = std::sqrt(2)/2 * mDeltaX; 

    // setting mStarAtom & mStarMinusOneAtom
    // n-1 to n vector intersects w/ Y axis
    if (std::abs(mDeltaX) > std::abs(mDeltaY))
    {   
        mStarAtom[0] = mNthAtom[0] - mDeltaX/2; // x component
        mStarAtom[1] = mNthAtom[1] - mDistance; // y component
        mStarAtom[2] = mNthAtom[2]; // z component
        if (!(isFound(mStarAtom[0], mStarAtom[1], mStarAtom[2])))
        {
            std::cout << "ERROR: setting StarAtom failed\n";
            isSet = false;
        }
        mStarMinusOneAtom[0] = mNthMinusOneAtom[0] - mDeltaX/2; // x component
        mStarMinusOneAtom[1] = mNthMinusOneAtom[1] - mDistance; // y component
        mStarMinusOneAtom[2] = mNthMinusOneAtom[2]; // z component
        if (!(isFound(mStarMinusOneAtom[0], mStarMinusOneAtom[1], mStarMinusOneAtom[2])))
        {
            std::cout << "ERROR: setting StarMinusOneAtom failed\n";
            isSet = false;
        }
    }
    else { std::cout << "ERROR: not an ASE generated input file"; }
/*
    const double tolerance = 0.20;
    for (int j=(mVector.size()-1); j>1; j--)
    {
        double temp [3];
        std::string unwanted;
        std::istringstream iss(mVector[j]);
        iss >> unwanted >> temp[0] >> temp[1] >> temp[2] >> unwanted;
        if (temp[2] < (mNthAtom[2]-tolerance))
        {
            mSecLayerZ = temp[2];
            break;
        }
    }
    for (int j=(mVector.size()-1); j>1; j--)
    {
        double temp [3];
        std::string unwanted;
        std::istringstream iss(mVector[j]);
        iss >> unwanted >> temp[0] >> temp[1] >> temp[2] >> unwanted;
        if (temp[2] < (mSecLayerZ-tolerance))
        {
            mThirLayerZ = temp[2];
            break;
        }
    }
*/

    return (isSet);
}

bool bcc110::isFound(const double &inX, const double &inY, const double &inZ) // check if the calculated * atoms exist
{
    const double tolerance = 0.15;
    for (int j=(mVector.size()-1); j>1; j--)
    {
        double temp [3];
        std::string unwanted;
        std::istringstream iss(mVector[j]);
        iss >> unwanted >> temp[0] >> temp[1] >> temp[2] >> unwanted;
        if (  ((temp[0] <= (inX + tolerance)) && (temp[0] >= (inX - tolerance)) ) &&
              ((temp[1] <= (inY + tolerance)) && (temp[1] >= (inY - tolerance)) ) &&   
              ((temp[2] <= (inZ + tolerance)) && (temp[2] >= (inZ - tolerance)) )  )   
        {
            return (true);
        }
    }
    return (false);
}

void bcc110::findHollow(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = 0;
    if (offsetY == 0)
    {
        offX = offsetX*mDeltaX + mDeltaX/2;
    }
    else if (offsetY > 0)
    {
        offX = (offsetX+offsetY) * mDeltaX;
    }
    double offY = offsetY*mDistance + mDistance/2;
    double hollowX = mNthAtom[0] - offX;
    double hollowY = mNthAtom[1] - offY;
    double hollowZ = mNthAtom[2] + m_DELTA_Z;
    if (hollowX < 0 || hollowY < 0)
    {
        std::cout << "ERROR: hollow offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(hollowX);
    std::string val2 = std::to_string(hollowY);
    std::string val3 = std::to_string(hollowZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc110-hollow.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findHollow

void bcc110::findAtop(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = 0;
    if (offsetY == 0)
    {
        offX = offsetX*mDeltaX;
    }
    else if (offsetY > 0)
    {
        offX = offsetX*mDeltaX + offsetY*mDeltaX/2;
    }
    double offY = offsetY*mDistance;
    double atopX = mNthAtom[0] - offX;
    double atopY = mNthAtom[1] - offY;
    double atopZ = mNthAtom[2] + m_DELTA_Z;
    if (atopX < 0 || atopY < 0)
    {
        std::cout << "ERROR: atop offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(atopX);
    std::string val2 = std::to_string(atopY);
    std::string val3 = std::to_string(atopZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc110-atop.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findAtop

void bcc110::findLongBridge(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = 0;
    if (offsetY == 0)
    {
        offX = offsetX*mDeltaX + mDeltaX/2;
    }
    else if (offsetY > 0)
    {
        offX = (offsetX+offsetY) * mDeltaX;
    }
    double offY = offsetY*mDistance;
    double LbrgX = mNthAtom[0] - offX;
    double LbrgY = mNthAtom[1] - offY;
    double LbrgZ = mNthAtom[2] + m_DELTA_Z;
    if (LbrgX < 0 || LbrgY < 0)
    {
         std::cout << "ERROR: long bridge offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(LbrgX);
    std::string val2 = std::to_string(LbrgY);
    std::string val3 = std::to_string(LbrgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc110-Longbrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findLongBridge

void bcc110::findShortBridge(const unsigned int offsetX, const unsigned int offsetY)
{
    double offX = 0;
    if (offsetY == 0)
    {
        offX = offsetX*mDeltaX + mDeltaX/4;
    }
    else if (offsetY > 0)
    {
        offX = offsetX*mDeltaX + (offsetY*3*mDeltaX/4);
    }
    double offY = offsetY*mDistance + mDistance/2;
    double SbrgX = mNthAtom[0] - offX;
    double SbrgY = mNthAtom[1] - offY;
    double SbrgZ = mNthAtom[2] + m_DELTA_Z;
    if (SbrgX < 0 || SbrgY < 0)
    {
        std::cout << "ERROR: short bridge offset out of the slab boundry" << std::endl;
    }

    std::string val1 = std::to_string(SbrgX);
    std::string val2 = std::to_string(SbrgY);
    std::string val3 = std::to_string(SbrgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("bcc110-Shortbrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findShortBridge
