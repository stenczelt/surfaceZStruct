// Mina Jafari
// 12-09-2015

#include "diamond100.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <vector>
#include <stdio.h>
#include <cmath>
#include <iomanip>

double diamond100::m_DELTA_Z = 2.5;

bool diamond100::setAtoms(const std::vector<std::string> &xyzFile)
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
        if (!(isFound(mStarAtom[0], mStarAtom[1], mStarAtom[2])))
        {
            std::cout << "ERROR: setting StarAtom failed\n";
            isSet = false;
        }
        mStarMinusOneAtom[0] = mNthMinusOneAtom[0]; // x component
        mStarMinusOneAtom[1] = mNthMinusOneAtom[1] - mDeltaX; // y component
        mStarMinusOneAtom[2] = mNthMinusOneAtom[2]; // z component
        if (!(isFound(mStarMinusOneAtom[0], mStarMinusOneAtom[1], mStarMinusOneAtom[2])))
        {
            std::cout << "ERROR: setting StarMinusOneAtom failed\n";
            isSet = false;
        }
    }   
    else { std::cout << "ERROR: not an ASE generated input file"; }

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

    return (isSet);
}

bool diamond100::isFound(const double &inX, const double &inY, const double &inZ)
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
} //isFound

void diamond100::findAtop()
{
    std::string val1 = std::to_string(mStarMinusOneAtom[0]);
    std::string val2 = std::to_string(mStarMinusOneAtom[1]);
    std::string val3 = std::to_string(mStarMinusOneAtom[2] + m_DELTA_Z);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("diamond100-atop.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findAtop

void diamond100::findFirstBridge()
{
    double brgX = mNthAtom[0] - (mDeltaX/2);
    double brgY = mNthAtom[1];
    double brgZ = mNthAtom[2] + m_DELTA_Z;

    std::string val1 = std::to_string(brgX);
    std::string val2 = std::to_string(brgY);
    std::string val3 = std::to_string(brgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("diamond100-firstBrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findFirstBridge

void diamond100::findSecondBridge()
{
    double brgX = mNthMinusOneAtom[0];
    double brgY = mNthMinusOneAtom[1] - (mDeltaX/2);
    double brgZ = mNthAtom[2] + m_DELTA_Z;

    std::string val1 = std::to_string(brgX);
    std::string val2 = std::to_string(brgY);
    std::string val3 = std::to_string(brgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("diamond100-secondBrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} //findSecondBridge

void diamond100::findThirdBridge()
{
    double brgX = mNthAtom[0] - (mDeltaX/2);
    double brgY = mNthAtom[1] - (mDeltaX/2);
    double brgZ = mSecLayerZ + m_DELTA_Z;

    std::string val1 = std::to_string(brgX);
    std::string val2 = std::to_string(brgY);
    std::string val3 = std::to_string(brgZ);
    std::string newElem = "C          " + val1 + "       " + val2 + "      " + val3;

    std::ofstream ofs;
    ofs.open ("diamond100-thirdBrg.xyz", std::ofstream::out);
    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
    ofs << "\n";
    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
    {   
        ofs << *i << "\n";
    }
    ofs << newElem;
    ofs.close();
} // findThirdBridge
