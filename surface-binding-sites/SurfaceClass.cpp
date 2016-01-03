// Mina Jafari
// 12-14-2015

#include "SurfaceClass.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

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

const BindingSiteClass* SurfaceClass::getBindingSite(unsigned int index) const //zero indexed
{
    if (index <= mBindingSites.size() && index >= 0)
    {
        return (&(mBindingSites[index]));
    }
    else
    {
        std::cout << "ERROR: The index requested is out of range" << std::endl;
        return (NULL);
    }
}

bool SurfaceClass::setSurfaceType(std::string inSurface)
{
    bool isSet = false;
    if (inSurface == "fcc100" || inSurface == "fcc110" || inSurface == "fcc111" ||
        inSurface == "bcc100" || inSurface == "bcc110" || inSurface == "bcc111" ||
        inSurface == "hcp0001")
    {
        mSurfaceType = inSurface;
        isSet = true;;
    }
    return (isSet);
}

bool SurfaceClass::setAtoms(int numOfAtoms, double* coordinates, std::string* atomicSymbols)
{
    resetGeometry();
    bool isSet = true;
    int k = 0;
    int numOfSurfAtoms = 0;
    std::string surfaceAtom = atomicSymbols[0];
    for (int i=0; i<numOfAtoms; ++i)
    {
        if (atomicSymbols[k] == surfaceAtom)
        {
            std::vector<double> temp;
            temp.push_back(coordinates[3*k]);
            temp.push_back(coordinates[3*k+1]);
            temp.push_back(coordinates[3*k+2]);
            mCoordinates.push_back(temp);
            mAtomicSymbols.push_back(atomicSymbols[k]);
            ++numOfSurfAtoms;
        }
        ++k;
    }
    mNumOfAtoms = numOfSurfAtoms;

    if (setSlabSize())
    {
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
            mDistance = mDeltaX;
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
        if (std::abs(mDeltaX) > std::abs(mDeltaY) && 
           (mSurfaceType == "fcc100" || mSurfaceType == "bcc100" || mSurfaceType == "fcc110"))
        {
            mStarAtom[0] = mNthAtom[0]; // x component
            mStarAtom[1] = mNthAtom[1] - mDistance; // y component
            mStarAtom[2] = mNthAtom[2]; // z component
            mStarMinusOneAtom[0] = mNthMinusOneAtom[0]; // x component
            mStarMinusOneAtom[1] = mNthMinusOneAtom[1] - mDistance; // y component
            mStarMinusOneAtom[2] = mNthMinusOneAtom[2]; // z component
        }
        else if (std::abs(mDeltaX) > std::abs(mDeltaY) && 
                (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001" || mSurfaceType == "bcc110"))
        {
            mStarAtom[0] = mNthAtom[0] - mDeltaX/2; // x component
            mStarAtom[1] = mNthAtom[1] - mDistance; // y component
            mStarAtom[2] = mNthAtom[2]; // z component
            mStarMinusOneAtom[0] = mNthMinusOneAtom[0] - mDeltaX/2; // x component
            mStarMinusOneAtom[1] = mNthMinusOneAtom[1] - mDistance; // y component
            mStarMinusOneAtom[2] = mNthMinusOneAtom[2]; // z component
        }
        else 
        { 
            std::cout << "ERROR: not an ASE generated input file" << std::endl;
            isSet = false;
        }
        if (!(isFound(mStarAtom[0], mStarAtom[1], mStarAtom[2])))
        {   
            std::cout << "ERROR: setting the StarAtom failed" << std::endl;
            isSet = false;
        }   
        if (!(isFound(mStarMinusOneAtom[0], mStarMinusOneAtom[1], mStarMinusOneAtom[2])))
        {   
            std::cout << "ERROR: setting the StarMinusOneAtom failed" << std::endl;
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
    }
    else
    {
        isSet = false;
    }
    return (isSet);
}

bool SurfaceClass::setSlabSize()
{
    bool isSet = true;
    const double tolerance = 0.15;
    int width = 0;
    int length = 0;
    int height = 0;
    int layer = 0;
    for (unsigned int i=0; i<mCoordinates.size()-1; ++i)
    {
        if (mCoordinates[i][1] <= mCoordinates[0][1]+tolerance &&
            mCoordinates[i][1] >= mCoordinates[0][1]-tolerance &&
            mCoordinates[i][2] <= mCoordinates[0][2]+tolerance &&
            mCoordinates[i][2] >= mCoordinates[0][2]-tolerance)
        {
            ++width;
        }
        if (mCoordinates[i][2] <= mCoordinates[0][2]+tolerance &&
            mCoordinates[i][2] >= mCoordinates[0][2]-tolerance)
        {
            ++layer;
        }
        length = layer / width;
        height = mNumOfAtoms / layer;
        if (height < 2)
        {
            std::cout << "ERROR: The slab should have at least 2 layers" << std::endl;
            isSet = false;
        }
    }

    mSlabSize[0] = width;
    mSlabSize[1] = length;
    mSlabSize[2] = height;
    return (isSet);
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

int SurfaceClass::findHollow()
{
    double m_DELTA_Z = 1.5;
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
                if (hollowX >= -0.05 && hollowY >= -0.05)
                {   
                    BindingSiteClass aSite("hollow", hollowX, hollowY, hollowZ);
                    mBindingSites.push_back(aSite);
                }   
            }
        }
    }
    else if (mSurfaceType == "fcc110") //TODO fcc110-hollow: include the layer above mNthAtom
    {
        for (int i=0; i<mSlabSize[0]-1; ++i)
        {
            for (int j=0; j<mSlabSize[1]-1; ++j)
            {
                double offX = i * mDeltaX + (mDeltaX/2);
                double offY = j * mDistance + (mDistance/2);
                double hollowX = mNthAtom[0] - offX;
                double hollowY = mNthAtom[1] - offY;
                double hollowZ = mNthAtom[2] + m_DELTA_Z;
                if (hollowX >= -0.05 && hollowY >= -0.05)
                {   
                    BindingSiteClass aSite("hollow", hollowX, hollowY, hollowZ);
                    mBindingSites.push_back(aSite);
                }
            }
        }
    }
    else if (mSurfaceType == "bcc110")
    {
        double limit = mSlabSize[0];
        double prevI = 0.5;
        for (int j=0; j<mSlabSize[1]; ++j)
        {
            for (double i=prevI; i<limit; i=i+0.5)
            {
                double offX = i * mDeltaX;
                double offY = j * mDistance + mDistance/2;
                double hollowX = mNthAtom[0] - offX;
                double hollowY = mNthAtom[1] - offY;
                double hollowZ = mNthAtom[2] + m_DELTA_Z;
                if (hollowX >= -0.05 && hollowY >= -0.05)
                {   
                    BindingSiteClass aSite("hollow", hollowX, hollowY, hollowZ);
                    mBindingSites.push_back(aSite);
                }
            }
            limit += 0.5;
            prevI += 0.5;
        }
    }
    else
    {
        std::cout << "ERROR: The HOLLOW binding site is not defined for this surface type" << std::endl;
        return (1200);
    }
    return (0);
}

int SurfaceClass::findHcp()
{
    double m_DELTA_Z = 2.0;
    double prevI_1 = 1.0;
    double limit_1 = mSlabSize[0]-1;
    double prevI_2 = 1.5;
    double limit_2 = mSlabSize[0];
    if (mSurfaceType == "fcc111" || mSurfaceType == "hcp0001")
    {
        for (int j=0; j<mSlabSize[1]; ++j)
        {
            for (double i=prevI_1; i<=limit_1; ++i)
            {
                double offX = i * mDeltaX;
                double offY = j * mDistance + 2*mDistance/3; // equilateral triangle
                double hcpX = mNthAtom[0] - offX;
                double hcpY = mNthAtom[1] - offY;
                double hcpZ = mNthAtom[2] + m_DELTA_Z;
                if (hcpX >= -0.05 && hcpY >= -0.05)
                {
                    BindingSiteClass aSite("hcp", hcpX, hcpY, hcpZ);
                    mBindingSites.push_back(aSite);
                }
            }
            prevI_1 += 0.5;
            limit_1 += 0.5;
        }
    }
    else if (mSurfaceType == "bcc111")
    {
        if (mSlabSize[2] < 3)
        {
            std::cout << "ERROR: This surface type requires at least three layers to define the HCP sites" << std::endl;
        }
        else
        {
            for (int j=0; j<mSlabSize[1]; ++j)
            {
                for (double i=prevI_2; i<limit_2; ++i)
                {
                    double offX = i * mDeltaX;
                    double offY = j * mDistance + 2*mDistance/3;
                    double hcpX = mNthAtom[0] - offX;
                    double hcpY = mNthAtom[1] - offY;
                    double hcpZ = mNthAtom[2] + m_DELTA_Z;
                    if (hcpX >= -0.05 && hcpY >= -0.05)
                    {
                        BindingSiteClass aSite("hcp", hcpX, hcpY, hcpZ);
                        mBindingSites.push_back(aSite);
                    }
                }
                prevI_2 += 0.5;
                limit_2 += 0.5;
            }
        }
    }
    else
    {
        std::cout << "ERROR: The HCP binding site is not defined for this surface type" << std::endl;
        return (1210);
    }
    if (mSurfaceType == "bcc111")
    {
        for (int j=0; j<1; ++j)
        {
            for (double i=-1; i<j+mSlabSize[0]; i=i+0.5)
            {
                double offX = i * mDeltaX;
                double offY = j * mDistance + mDistance/3;
                double hcpX = mNthAtom[0] - offX;
                double hcpY = mNthAtom[1] + offY;
                double hcpZ = mNthAtom[2] + m_DELTA_Z;
                if (mSlabSize[2] > 2)
                {
                    if (isFound(hcpX, hcpY, mThirdLayerZ))
                    {
                        if (hcpX >= -0.05 && hcpY >= -0.05)
                        {
                            BindingSiteClass aSite("hcp", hcpX, hcpY, hcpZ);
                            mBindingSites.push_back(aSite);
                        }
                    }
                }
            }
        }
    }
    return (0);
} // findHcp

int SurfaceClass::findFcc()
{
    double m_DELTA_Z = 2.0;
    if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        double prevI = 0.5;
        double limit = mSlabSize[0]-1;
        for (int j=0; j<mSlabSize[1]-1; ++j) // j is the Y offset
        {
            for (double i=prevI; i<limit; ++i) // i is the X offset
            {
                double offX = i * mDeltaX;
                double offY = j * mDistance + mDistance/3;
                double fccX = mNthAtom[0] - offX;
                double fccY = mNthAtom[1] - offY;
                double fccZ = mNthAtom[2] + m_DELTA_Z;
                if (isFound(fccX, fccY, mSecondLayerZ))
                {
                    BindingSiteClass aSite("fcc", fccX, fccY, fccZ);
                    mBindingSites.push_back(aSite);
                }
                else
                {   
                    //std::cout << "ERROR: not able to find the FCC site" << std::endl;
                }
            } // inner for (i)
            prevI += 0.5;
            limit += 0.5;
        } // outer for (j)
    } // outer if
    else
    {
        std::cout << "ERROR: The FCC binding site is not defined for this surface type" << std::endl;
        return (1220);
    }
    if (mSurfaceType == "bcc111")
    {
        for (int j=0; j<1; ++j)
        {
            for (double i=0; i<j+mSlabSize[0]; i=i+0.5)
            {
                double offX = i * mDeltaX;
                double offY = j * mDistance + 2*mDistance/3;
                double fccX = mNthAtom[0] - offX;
                double fccY = mNthAtom[1] + offY;
                double fccZ = mNthAtom[2] + m_DELTA_Z;
                if (mSlabSize[2] > 2)
                {
                    if (isFound(fccX, fccY, mThirdLayerZ))
                    {
                        if (fccX >= -0.05 && fccY >= -0.05)
                        {
                            BindingSiteClass aSite("fcc", fccX, fccY, fccZ);
                            mBindingSites.push_back(aSite);
                        }
                    }
                }
            }
        }
    }
    return (0);
} // findFcc

int SurfaceClass::findAtop()
{
    double m_DELTA_Z = 2.5;
    for (int i=0; i<mSlabSize[0]; ++i) // i is the X offset
    {
        for (int j=0; j<mSlabSize[1]; ++j) // j is the Y offset
        {
            double offX, offY, atopX, atopY, atopZ = 0.0;
            if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
            {
                offX = i * mDeltaX;
                offY = j * mDeltaX;
            }
            else if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
            {
                offX = (i * mDeltaX) + (j * mDeltaX/2);
                offY = j * mDistance;
            }
            else if (mSurfaceType == "fcc110")
            {
                offX = i * mDeltaX;
                offY = j * mDistance;
            }
            else if (mSurfaceType == "bcc110")
            {
                if (j == 0)
                {   
                    offX = i * mDeltaX;
                }   
                else if (j > 0)
                {   
                    offX = i * mDeltaX + j * mDeltaX/2;
                }   
                offY = j * mDistance;
            }
            else
            {
                std::cout << "ERROR: The ATOP binding site is not defined for this surface type" << std::endl;
                return (1230);
            }
            atopX = mNthAtom[0] - offX;
            atopY = mNthAtom[1] - offY;
            atopZ = mNthAtom[2] + m_DELTA_Z;

            if (atopX >= -0.05 && atopY >= -0.05)
            {
                BindingSiteClass aSite("atop", atopX, atopY, atopZ);
                mBindingSites.push_back(aSite);
            }
        } // inner for (j)
    } // outer for (i)
    return (0);
} // findAtop

int SurfaceClass::findLongBridge()
{
    double m_DELTA_Z = 2.0;
    double offX, offY, LbrgX, LbrgY, LbrgZ = 0.0;
    if (mSurfaceType == "fcc110")
    {
        for (int i=0; i<mSlabSize[0]; ++i) // i is the X offset
        {
            for (int j=0; j<mSlabSize[1]; ++j) // j is the Y offset
            {
                offX = i * mDeltaX + mDeltaX/2;
                offY = j * mDistance;

                LbrgX = mNthAtom[0] - offX;
                LbrgY = mNthAtom[1] - offY;
                LbrgZ = mNthAtom[2] + m_DELTA_Z;
                if (LbrgX >= -0.05 && LbrgY >= -0.05)
                {
                    BindingSiteClass aSite("long-bridge", LbrgX, LbrgY, LbrgZ);
                    mBindingSites.push_back(aSite);
                }
            }
        }
    }
    else if (mSurfaceType == "bcc110")
    {
        double prevI_1 = 0.5;
        int prevI_2 = 1;
        double limit_1 = mSlabSize[0]-1;
        int limit_2 = mSlabSize[0];
        for (int j=0; j<mSlabSize[1]; ++j)
        {
            if (j%2 == 0)
            {
                for (double i=prevI_1; i<limit_1; ++i)
                {
                    offX = i * mDeltaX;
                    offY = j * mDistance;

                    LbrgX = mNthAtom[0] - offX;
                    LbrgY = mNthAtom[1] - offY;
                    LbrgZ = mNthAtom[2] + m_DELTA_Z;
                    if (LbrgX >= -0.05 && LbrgY >= -0.05)
                    {
                        BindingSiteClass aSite("long-bridge", LbrgX, LbrgY, LbrgZ);
                        mBindingSites.push_back(aSite);
                    }
                }
                ++limit_1;
                ++prevI_1;
            }
            else  if (j%2 != 0)
            {
                for (int i=prevI_2; i<limit_2; ++i)
                {
                    offX = i * mDeltaX;
                    offY = j * mDistance;

                    LbrgX = mNthAtom[0] - offX;
                    LbrgY = mNthAtom[1] - offY;
                    LbrgZ = mNthAtom[2] + m_DELTA_Z;
                    if (LbrgX >= -0.05 && LbrgY >= -0.05)
                    {
                        BindingSiteClass aSite("long-bridge", LbrgX, LbrgY, LbrgZ);
                        mBindingSites.push_back(aSite);
                    }
                }
                ++limit_2;
                ++prevI_2;
            }
        }
    }
    else
    {
        std::cout << "ERROR: The LONG-BRIDGE binding site is not defined for this surface type" << std::endl;
        return (1240);
    }
    return (0);
} // findLongBridge

int SurfaceClass::findShortBridge()
{
    double m_DELTA_Z = 2.0;
    double offX, offY, SbrgX, SbrgY, SbrgZ = 0.0;
    if (mSurfaceType == "fcc110")
    {
        for (int i=0; i<mSlabSize[0]; ++i) // i is the X offset
        {
            for (int j=0; j<mSlabSize[1]-1; ++j) // j is the Y offset
            {
                offX = i * mDeltaX;
                offY = j * mDistance + mDistance/2;;

                SbrgX = mNthAtom[0] - offX;
                SbrgY = mNthAtom[1] - offY;
                SbrgZ = mNthAtom[2] + m_DELTA_Z;

                if (SbrgX >= -0.05 || SbrgY > -0.05)
                {
                    BindingSiteClass aSite("short-bridge", SbrgX, SbrgY, SbrgZ);
                    mBindingSites.push_back(aSite);
                }
            }
        }
    }
    else if (mSurfaceType == "bcc110")
    {
        double prevI = 0.25;
        double limit = mSlabSize[0]-1+0.25;
        for (double j=0.5; j<mSlabSize[1]-1; ++j)
        {
            for (double i=prevI; i<=limit; i=i+0.5)
            {
                offX = i * mDeltaX;
                offY = j * mDistance;

                SbrgX = mNthAtom[0] - offX;
                SbrgY = mNthAtom[1] - offY;
                SbrgZ = mNthAtom[2] + m_DELTA_Z;

                if (SbrgX >= -0.05 || SbrgY > -0.05)
                {
                    BindingSiteClass aSite("short-bridge", SbrgX, SbrgY, SbrgZ);
                    mBindingSites.push_back(aSite);
                }
            }
            prevI += 0.5;
            limit += 0.5;
        }
    }
    else
    {
        std::cout << "ERROR: The SHORT-BRIDGE binding site is not defined for this surface type" << std::endl;
        return (1250);
    }
    return (0);
} // findShortBridge

int SurfaceClass::findBridge()
{
    double m_DELTA_Z = 2.0;
    double offX, offY, brgX, brgY = 0.0;
    double brgZ = mNthAtom[2] + m_DELTA_Z;
    if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
    {
        for (double i=0.0; i<=mSlabSize[0]; i=i+0.5) // i is the X offset
        {
            for (double j=0.0; j<=mSlabSize[1]; j=j+0.5) // j is the Y offset
            {
                if ((i==std::floor(i) && j!=std::floor(j)) || (i!=std::floor(i) && j==std::floor(j)))
                {
                    offX = i * mDeltaX;
                    offY = j * mDeltaX;
                    brgX = mNthAtom[0] - offX;
                    brgY = mNthAtom[1] - offY;
                    if (brgX >= -0.05 && brgY >= -0.05)
                    {
                        BindingSiteClass aSite("bridge", brgX, brgY, brgZ);
                        mBindingSites.push_back(aSite);
                    }
                }
            }
        }
    }
    else if (mSurfaceType == "fcc111" || mSurfaceType == "hcp0001")
    {
        double prevI_1 = 0.5;
        double prevI_2 = 0.5;
        double prevI_3 = 0.25;
        for (double j=0.0; j<=mSlabSize[1]; j=j+0.5)
        //for (int j=0.0; j<=mSlabSize[1]; ++j)
        {
            if (j == std::floor(j))
            {
                int k = j;
                if (k%2 == 0)
                {
                    for (double i=prevI_1; i<(prevI_1+mSlabSize[0]-1); ++i)
                    {
                        offX = i * mDeltaX;
                        offY = k * mDistance;
                        brgX = mNthAtom[0] - offX;
                        brgY = mNthAtom[1] - offY;
                        if ((brgX >= -0.05 && brgY >= -0.05) && (brgX <= mNthAtom[0]-k/2*mDeltaX && brgY <= mNthAtom[1]))
                        {
                            BindingSiteClass aSite("bridge", brgX, brgY, brgZ);
                            mBindingSites.push_back(aSite);
                        }
                    }
                    ++prevI_1;
                }
                else if (k%2 != 0)
                {
                    for (double i=prevI_2+0.5; i<(prevI_2+mSlabSize[0]-1); ++i)
                    {
                        offX = i * mDeltaX;
                        offY = k * mDistance;
                        brgX = mNthAtom[0] - offX;
                        brgY = mNthAtom[1] - offY;
                        if ((brgX >= -0.05 && brgY >= -0.05) && (brgX <= mNthAtom[0]-(k+1)/2*mDeltaX && brgY <= mNthAtom[1]))
                        {
                            BindingSiteClass aSite("bridge", brgX, brgY, brgZ);
                            mBindingSites.push_back(aSite);
                        }
                    }
                    ++prevI_2;
                }
            }
            else if (j != std::floor(j))
            {
                for (double i=prevI_3; i<(prevI_3+mSlabSize[0]-0.5); i=i+0.5)
                {
                    offX = i * mDeltaX;
                    offY = j * mDistance;
                    brgX = mNthAtom[0] - offX;
                    brgY = mNthAtom[1] - offY;
                    if ((brgX >= -0.05 && brgY >= -0.05) && (brgX <= mNthAtom[0]-i*mDeltaX && brgY <= mNthAtom[1]))
                    {
                        BindingSiteClass aSite("bridge", brgX, brgY, brgZ);
                        mBindingSites.push_back(aSite);
                    }
                }
                prevI_3 += 0.5;
            }
        }
    }
    else
    {
        std::cout << "ERROR: The BRIDGE binding site is not defined for this surface type" << std::endl;
        return (1260);
    }
    return (0);
} // findBridge

void SurfaceClass::findNearbySites(const int atomIndex, const double radius, 
                                   const std::string siteType)
{
    if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
    {
        // TODO if (findHollow() != 0) {ERROR}
        findHollow();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        findHcp();
        findFcc();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == "fcc110" || mSurfaceType == "bcc110")
    {
        findHollow();
        findAtop();
        findLongBridge();
        findShortBridge();
    }
    else
    {
        std::cout << "ERROR: not a supported surface type" << std::endl;
    }
    // here, we have the mBindingSites vector with all the binding sites
    int numOfSites = mBindingSites.size();
    for (int i=0; i<numOfSites; ++i)
    {
        if (mBindingSites[i].getType() == siteType || siteType == "all")
        {
            double refX = mCoordinates[atomIndex-1][0];
            double refY = mCoordinates[atomIndex-1][1];
            double testX = mBindingSites[i].getX();
            double testY = mBindingSites[i].getY();
            if ((testX <= refX+radius && testX >= refX-radius) &&
                (testY <= refY+radius && testY >= refY-radius) &&
               !(testX <= refX+0.1 && testX >= refX-0.1 &&
                 testY <= refY+0.1 && testY >= refY-0.1) )
            {
                // the site is in the range
                mSelectedBindingSites.push_back(mBindingSites[i]);
                std::cout << "This site IS within the specified radius/type" << std::endl;
            }
            else
            {
                std::cout << "ERROR: This site is NOT within the specified radius/type" << std::endl;
            }
        }
    }
}

void SurfaceClass::findAllSites()
{
    if (mSurfaceType == "fcc100" || mSurfaceType == "bcc100")
    {
        // TODO if (findHollow() != 0) {ERROR}
        findHollow();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == "fcc111" || mSurfaceType == "bcc111" || mSurfaceType == "hcp0001")
    {
        findHcp();
        findFcc();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == "fcc110" || mSurfaceType == "bcc110")
    {
        findHollow();
        findAtop();
        findLongBridge();
        findShortBridge();
    }
    else
    {
        std::cout << "ERROR: not a supported surface type" << std::endl;
    }
    // here, we have the mBindingSites vector with all the binding sites
    int numOfSites = mBindingSites.size();
    for (int i=0; i<numOfSites; ++i)
    {
        mSelectedBindingSites.push_back(mBindingSites[i]);
        std::cout << "This site IS within the specified radius/type" << std::endl;
    }
}

bool SurfaceClass::writeToFile(std::string &outFile)
{
    bool success = false;
    std::ofstream ofs;
    ofs.open(outFile.c_str());

    ofs << std::to_string(mNumOfAtoms+mSelectedBindingSites.size()) << "\n";
    ofs << "\n";
    ofs << std::fixed << std::setprecision(15);
    for (unsigned int i=0; i<mCoordinates.size(); ++i)
    {   
        ofs << mAtomicSymbols[i] << "            " << mCoordinates[i][0]
            << "            " << mCoordinates[i][1]
            << "            " << mCoordinates[i][2] << "\n";
    }
    for (unsigned int j=0; j<mSelectedBindingSites.size(); ++j)
    {
        ofs << "x             ";
        ofs << mSelectedBindingSites[j].getX() << "            ";
        ofs << mSelectedBindingSites[j].getY() << "            ";
        ofs << mSelectedBindingSites[j].getZ() << "\n";
    }
    ofs.close();
    success = true;
    return (success);
}

void SurfaceClass::resetGeometry()
{
//    mAtomicSymbols.clear();
    std::vector<BindingSiteClass> swap1;
    mBindingSites.swap(swap1);
    mSelectedBindingSites.swap(swap1);
    std::vector<std::string> swap2;
    mAtomicSymbols.swap(swap2);
    std::vector< std::vector<double> > swap3;
    mCoordinates.swap(swap3);
}
