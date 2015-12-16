// Mina Jafari
// 12-14-2015

#include "SurfaceClass.h"
#include <iostream> //TODO for testing
#include <cmath>
#include <fstream>
#include <iomanip>

double SurfaceClass::m_DELTA_Z = 2.5;

std::string SurfaceClass::getSurfaceType() const
{
    return (mSurfaceType);
}

int SurfaceClass::getNumOfAtoms() const
{
    return (mNumOfAtoms);
}

int SurfaceClass::getSurfaceWidth() const
{
    return (mSlabSize[0]);
}

int SurfaceClass::getSurfaceLength() const
{
    return (mSlabSize[1]);
}

int SurfaceClass::getSurfaceHeight() const
{
    return (mSlabSize[2]);
}

BindingSiteClass SurfaceClass::getBindingSite(unsigned int element) const
{
    if (element <= mBindingSites.size())
    {
        return (mBindingSites[element-1]); // TODO correct??
    }
    else
    {
        std::cout << "out of range" << std::endl;
    }
}

bool SurfaceClass::setSurfaceType(std::string inSurface)
{
    mSurfaceType = inSurface;
    return (true);
}

bool SurfaceClass::setAtoms(int numOfAtoms, double* coordinates, std::string* atomicSymbols)
{
    bool isSet = true;
    mNumOfAtoms = numOfAtoms;
    int k =0;
    for (int i=0; i<numOfAtoms; ++i)
    {
        std::vector<double> temp;
        temp.push_back(coordinates[3*k]);
        temp.push_back(coordinates[3*k+1]);
        temp.push_back(coordinates[3*k+2]);
        mCoordinates.push_back(temp);
        mAtomicSymbols.push_back(atomicSymbols[k]);
        ++k;
    }
    
    setSlabSize();

    // setting mNthAtom & mNumOfAtoms
    mNthAtom[0] = mCoordinates[mNumOfAtoms-1][0];
    mNthAtom[1] = mCoordinates[mNumOfAtoms-1][1];
    mNthAtom[2] = mCoordinates[mNumOfAtoms-1][2];
    mNthMinusOneAtom[0] = mCoordinates[mNumOfAtoms-2][0];
    mNthMinusOneAtom[1] = mCoordinates[mNumOfAtoms-2][1];
    mNthMinusOneAtom[2] = mCoordinates[mNumOfAtoms-2][2];
    mDeltaX = mNthAtom[0] - mNthMinusOneAtom[0];
    mDeltaY = mNthAtom[1] - mNthMinusOneAtom[1];
    if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        mDistance = std::sqrt(3)/2 * mDeltaX;
    }
    else if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
    {
        mDistance = 0.0;
    }
    else if (mSurfaceType == "fcc110")
    {
        mDistance = mDeltaX / std::sqrt(2);
    }
    else if (mSurfaceType == "bcc110")
    {
        mDistance = std::sqrt(2)/2 * mDeltaX;
    }
    // setting mStarAtom & mStarMinusOneAtom
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
    else 
    { 
        std::cout << "ERROR: not an ASE generated input file";
        isSet = false;
    }

    const double tolerance = 0.20;
    for (int j=(mCoordinates.size()-1); j>-1; j--)
    {
        if (mCoordinates[j][2] < (mNthAtom[2]-tolerance))
        {   
            mSecondLayerZ = mCoordinates[j][2];
            break;
        }   
    } 
    for (int j=(mCoordinates.size()-1); j>-1; j--)
    {   
        if (mCoordinates[j][2] < (mSecondLayerZ-tolerance))
        {   
            mThirdLayerZ = mCoordinates[j][2];
            break;
        }   
    }   
    return (isSet);
}

//void SurfaceClass::setSlabSize(const std::vector< std::vector<double> > &coordinates)
void SurfaceClass::setSlabSize()
{
    const double tolerance = 0.15;
    int width = 0;
    int length = 0;
    int height = 0;
    for (unsigned int i=0; i<mCoordinates.size()-1; ++i)
    {
        if (mCoordinates[i][1] <= mCoordinates[0][1]+tolerance &&
            mCoordinates[i][1] >= mCoordinates[0][1]-tolerance &&
            mCoordinates[i][2] == mCoordinates[0][2])
        {
            ++width;
        }
        if (mCoordinates[i][0] <= mCoordinates[0][0]+tolerance &&
            mCoordinates[i][0] >= mCoordinates[0][0]-tolerance &&
            mCoordinates[i][2] == mCoordinates[0][2])
        {
            ++length;
        }
        height = mNumOfAtoms / (width*length);
    }

    mSlabSize[0] = width;
    mSlabSize[1] = length;
    mSlabSize[2] = height;
}

bool SurfaceClass::isFound(const double &inX, const double &inY, const double &inZ)
{
    const double tolerance = 0.15;
    for (unsigned int i=0; i<mCoordinates.size(); ++i)
    {
        if (  ((mCoordinates[i][0] <= (inX + tolerance)) && (mCoordinates[i][0] >= (inX - tolerance)) ) &&
              ((mCoordinates[i][1] <= (inY + tolerance)) && (mCoordinates[i][1] >= (inY - tolerance)) ) &&   
              ((mCoordinates[i][2] <= (inZ + tolerance)) && (mCoordinates[i][2] >= (inZ - tolerance)) )  )   
        {   
            return (true);
        }   
    }
    return (false);
}

bool SurfaceClass::writeToFile()
{
    return (true);
}

void SurfaceClass::findHollow()
{
    if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
    {
        for (int i=0; i<mSlabSize[0]-1; ++i)
        {
            for (int j=0; j<mSlabSize[1]-1; ++j) // subtract 1 bc of # of defined sites
            {
                double offX = i * mDeltaX + (mDeltaX/2);
                double offY = j * mDeltaX + (mDeltaX/2);
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
                std::string newElem = "X          " + val1 + "       " + val2 + "      " + val3;
                BindingSiteClass aSite("hollow", hollowX, hollowY, hollowZ);
                mBindingSites.push_back(aSite);
            }
        }
//        std::ofstream ofs;
//        ofs.open ("Termination100-hollow.xyz", std::ofstream::out);

//        ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
//        ofs << "\n";
//        for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
//        {   
//            ofs << *i << "\n";
//        }   
//        ofs << newElem;
//        ofs << "\n";
//        ofs.close();
    }
    else
    {
        std::cout << "ERROR: Unsupported surface type" << std::endl;
    }
}

void SurfaceClass::findHcp()
{
    if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        for (int i=0; i<mSlabSize[0]; ++i)
        {
            for (int j=0; j<mSlabSize[1]-1; ++j) // subtract 1 bc of # of defined sites
            {
                double offX = 0;
                if (j > 0)
                {
                    offX = (i+1) * mDeltaX + mDeltaX/2;
                }
                else if (j == 0)
                {
                    offX = i * mDeltaX + mDeltaX;
                }
                double offY = j * mDistance + 2*mDistance/3;
                double hcpX = mNthAtom[0] - offX;
                double hcpY = mNthAtom[1] - offY; // equilateral triangle, third layer
                double hcpZ = mNthAtom[2] + m_DELTA_Z;
                if (hcpX < 0 || hcpY < 0)
                {   
                    std::cout << "ERROR: offset out of the slab boundry" << std::endl;
                }   
                //    This generates error for hcp0001, find a better if condition
                //    if ( (!(isFound(hcpX, hcpY, mThirLayerZ)) && (mThirLayerZ == 0.0)) || (isFound(hcpX, hcpY, mThirLayerZ)) )
                if ( (!(isFound(hcpX, hcpY, mThirdLayerZ)) && (mThirdLayerZ == 0.0)) || (isFound(hcpX, hcpY, mThirdLayerZ))
                     || (!(isFound(hcpX, hcpY, mThirdLayerZ)) && (mThirdLayerZ > 0.0)) )
                {   
                    std::string val1 = std::to_string(hcpX);
                    std::string val2 = std::to_string(hcpY);
                    std::string val3 = std::to_string(hcpZ);
                    std::string newElem = "X          " + val1 + "       " + val2 + "      " + val3;
                    BindingSiteClass aSite("hcp", hcpX, hcpY, hcpZ);
                    mBindingSites.push_back(aSite);

/*                    std::ofstream ofs;
                    ofs.open ("Termination111-hcp.xyz", std::ofstream::out);
                    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
                    ofs << "\n";
                    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
                    {   
                    ofs << *i << "\n";
                    }   
                    ofs << newElem;
                    ofs.close(); */
                }
                else
                {   
                    std::cout << "ERROR: not able to find the HCP site" << std::endl;
                }  
            } // inner fo r(j)
        } // outer for (i)
    } // outer if
    else
    {
        // not defined for other surface sites
    }
} // findHcp
    
void SurfaceClass::findFcc()
{
    if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        for (int i=0; i<mSlabSize[0]; ++i) // i is the X offset
        {
            for (int j=0; j<mSlabSize[1]-1; ++j) // subtract 1 bc of # of defined sites
            {
                double offX = 0;
                if (j == 0)
                {
                    offX = (i * mDeltaX) + mDeltaX/2;
                }
                else if (j > 0)
                {
                    offX = (i+j) * mDeltaX;
                }
                double offY = j * mDistance + mDistance/3;
                double fccX = mNthAtom[0] - offX;
                double fccY = mNthAtom[1] - offY; // equilateral triangle
                double fccZ = mNthAtom[2] + m_DELTA_Z;
                if (fccX < 0 || fccY < 0)
                {   
                    std::cout << "ERROR: offset out of the slab boundry" << std::endl;
                }   
                if (isFound(fccX, fccY, mSecondLayerZ))
                {   
                    std::string val1 = std::to_string(fccX);
                    std::string val2 = std::to_string(fccY);
                    std::string val3 = std::to_string(fccZ);
                    std::string newElem = "X          " + val1 + "       " + val2 + "      " + val3;
                    BindingSiteClass aSite("fcc", fccX, fccY, fccZ);
                    mBindingSites.push_back(aSite);

/*                    std::ofstream ofs;
                    ofs.open ("Termination111-fcc.xyz", std::ofstream::out);
                    ofs << std::to_string( atoi(mVector[0].c_str()) + 1 ) << "\n";
                    ofs << "\n";
                    for (auto i = mVector.begin()+2; i != mVector.end(); ++i)
                    {   
                    ofs << *i << "\n";
                    }   
                    ofs << newElem;
                    ofs.close(); */
                }
                else
                {   
                    std::cout << "ERROR: not able to find the HCP site" << std::endl;
                }  
            } // inner fo r(j)
        } // outer for (i)
    } // outer if
} // findFcc

void SurfaceClass::SurfaceClass::findAtop()
{
}

void SurfaceClass::findLongBridge()
{
}

void SurfaceClass::findShortBridge()
{
}

void SurfaceClass::findBridge()
{
}

void SurfaceClass::findNearbySites(int atomIndex, double radius, std::string siteType)
{
    if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
    {
        findHollow();
        //findAtop();
        //findBridge();
    }
    else if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        findHcp();
        findFcc();
        //findAtop();
        //findBridge();
    }
    else if (mSurfaceType == "fcc110" || mSurfaceType == "bcc110")
    {
        findHollow();
        //findAtop();
        //findLongBridge();
        //findShortBridge();
    }
    else
    {
        std::cout << "ERROR: not a supported surface type" << std::endl;
    }
    // here, we have the mBindingSites vector with all the binding sites
    int numOfSites = mBindingSites.size();
    for (int i=0; i<numOfSites; ++i)
    {
        if (mBindingSites[i].getType() == siteType)
        {
            if (mBindingSites[i].getX() <= mCoordinates[atomIndex-1][0]+radius &&
                mBindingSites[i].getX() >= mCoordinates[atomIndex-1][0]-radius &&
                mBindingSites[i].getY() <= mCoordinates[atomIndex-1][1]+radius &&
                mBindingSites[i].getY() >= mCoordinates[atomIndex-1][1]-radius)
            {
                // the site is in the range
                mSelectedBindingSites.push_back(mBindingSites[i]);
            }
            else
            {
                std::cout << "ERROR: no site found within the specified radius/type" << std::endl;
            }
        }
    }

}

bool SurfaceClass::writeToFile(std::string &outFile) // HERE
{
    bool success = false;
    std::ofstream ofs;
    ofs.open(outFile.c_str());

    ofs << std::to_string(mNumOfAtoms+mSelectedBindingSites.size()) << "\n";
    ofs << "\n";
    ofs << std::fixed << std::setprecision(15);
    for (int i=0; i<mCoordinates.size(); ++i)
    {   
        ofs << mAtomicSymbols[i] << "            " << mCoordinates[i][0]
            << "            " << mCoordinates[i][1]
            << "            " << mCoordinates[i][2] << "\n";
    }
    for (int j=0; j<mSelectedBindingSites.size(); ++j)
    {
        ofs << "x             ";
        ofs << mSelectedBindingSites[j].getX() << "            "; // returnin double or int?
        ofs << mSelectedBindingSites[j].getY() << "            ";
        ofs << mSelectedBindingSites[j].getZ() << "\n";
    }
    ofs.close();
    success = true;
    return (success);
}
