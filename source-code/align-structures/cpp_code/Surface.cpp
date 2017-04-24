// Mina Jafari
// 12-14-2015

#include "Surface.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

Surface::Surface()
{
}

Surface::Surface(std::vector<Atom> slabAtoms, std::vector<Atom> bindingSites,
        std::vector<Molecule> adsorbateList)
{
    //TODO
    //This is a hack, Fix it
    mSlabAtoms = slabAtoms;
    mBindingSitesTemp = bindingSites;
    mAdsorbates = adsorbateList;
}

/*Surface::Surface(std::vector<Atom> slabAtoms, std::vector<BindingSite> bindingSites,
        std::vector<Molecule> adsorbateList)
{
    mSlabAtoms = slabAtoms;
    mBindingSites = bindingSites;
    mAdsorbates = adsorbateList;
}*/

SLAB_TYPE Surface::getSurfaceType() const
{
    return (mSurfaceType);
}

int Surface::getNumOfAtoms() const
{
    return (mNumOfSurfAtoms + mNumOfAdsorbateAtoms);
}

int Surface::getSurfaceWidth() const
{
    return (mSlabSize[0]);
}

int Surface::getSurfaceLength() const
{
    return (mSlabSize[1]);
}

int Surface::getSurfaceHeight() const
{
    return (mSlabSize[2]);
}

const BindingSite* Surface::getBindingSite(unsigned int index) const //zero indexed
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

void Surface::setSurfaceType(SLAB_TYPE inSurface)
{
    mSurfaceType = inSurface;
}

bool Surface::setAtoms(int numOfAtoms, double* coordinates, std::string* atomicSymbols)
{
    resetGeometry();
    bool isSet = true;
    std::string surfaceAtom = atomicSymbols[0];
    for (int i=0; i<numOfAtoms; ++i)
    {
        if (atomicSymbols[i] == surfaceAtom)
        {
            Atom anAtom(atomicSymbols[i], coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2]);
            mSlabAtoms.push_back(anAtom);
        }
        else
        {
            Atom anAtom(atomicSymbols[i], coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2]);
            mAdsorbateAtoms.push_back(anAtom);
            //std::cout << "ERROR: Reading surfaces with adsorbates not implemented" << std::endl;
            // TODO if slab is not clean and already has adsorbates attached to it, it should be taken
            // care of at this part of the code. We need that information only for finding empty binding sites.
            //Atom anAtom(atomicSymbols[i], coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2]);
            //create Molecule
            //mAdsorbates.push_back(anAtom); //push back molecules
            /*std::vector<double> temp;
            temp.push_back(coordinates[3*k]);
            temp.push_back(coordinates[3*k+1]);
            temp.push_back(coordinates[3*k+2]);
            mAdsorbateCoord.push_back(temp);
            mAdsorbateSymbols.push_back(atomicSymbols[k]);
            ++numOfAdsorbateAtoms;*/
            //isSet = false;
        }
    }
    mNumOfSurfAtoms = mSlabAtoms.size();
    //mNumOfAdsorbateAtoms = numOfAdsorbateAtoms;
    mNumOfAdsorbateAtoms = mAdsorbateAtoms.size();

    if (mSurfaceType != ANY)
    {
        if (setSlabSize())
        {
            // setting mNthAtom & mNumOfSurfAtoms
            mNthAtom[0] = mSlabAtoms[mNumOfSurfAtoms-1].coordinates().x();
            mNthAtom[1] = mSlabAtoms[mNumOfSurfAtoms-1].coordinates().y();
            mNthAtom[2] = mSlabAtoms[mNumOfSurfAtoms-1].coordinates().z();
            mNthMinusOneAtom[0] = mSlabAtoms[mNumOfSurfAtoms-2].coordinates().x();
            mNthMinusOneAtom[1] = mSlabAtoms[mNumOfSurfAtoms-2].coordinates().y();
            mNthMinusOneAtom[2] = mSlabAtoms[mNumOfSurfAtoms-2].coordinates().z();
            mDeltaX = mNthAtom[0] - mNthMinusOneAtom[0];
            mDeltaY = mNthAtom[1] - mNthMinusOneAtom[1];
            if (mSurfaceType == FCC111 || mSurfaceType == BCC111 || mSurfaceType == HCP0001)
            {
                mDistance = std::sqrt(3)/2 * mDeltaX;
            }
            else if (mSurfaceType == FCC100 || mSurfaceType == BCC100)
            {
                mDistance = mDeltaX;
            }
            else if (mSurfaceType == FCC110)
            {
                mDistance = mDeltaX / std::sqrt(2);
            }
            else if (mSurfaceType == BCC110)
            {
                mDistance = std::sqrt(2)/2 * mDeltaX;
            }
            // setting mStarAtom & mStarMinusOneAtom
            // n-1 to n vector intersects w/ Y axis
            if (std::abs(mDeltaX) > std::abs(mDeltaY) && 
                    (mSurfaceType == FCC100 || mSurfaceType == BCC100 || mSurfaceType == FCC110))
            {
                mStarAtom[0] = mNthAtom[0]; // x component
                mStarAtom[1] = mNthAtom[1] - mDistance; // y component
                mStarAtom[2] = mNthAtom[2]; // z component
                mStarMinusOneAtom[0] = mNthMinusOneAtom[0]; // x component
                mStarMinusOneAtom[1] = mNthMinusOneAtom[1] - mDistance; // y component
                mStarMinusOneAtom[2] = mNthMinusOneAtom[2]; // z component
            }
            else if (std::abs(mDeltaX) > std::abs(mDeltaY) && 
                    (mSurfaceType == FCC111 || mSurfaceType == BCC111 || mSurfaceType == HCP0001 || mSurfaceType == BCC110))
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
            for (int j=(mSlabAtoms.size()-1); j>-1; j--)
            {
                if (mSlabAtoms[j].coordinates().z() < (mNthAtom[2]-tolerance))
                {   
                    mSecondLayerZ = mSlabAtoms[j].coordinates().z();
                    break;
                }   
            } 
            for (int j=(mSlabAtoms.size()-1); j>-1; j--)
            {   
                if (mSlabAtoms[j].coordinates().z() < (mSecondLayerZ-tolerance))
                {   
                    mThirdLayerZ = mSlabAtoms[j].coordinates().z();
                    break;
                }   
            }
        }
        else
        {
            isSet = false;
        }
    }
    return (isSet);
}

bool Surface::setSlabSize()
{
    bool isSet = true;
    const double tolerance = 0.15;
    int width = 0;
    int length = 0;
    int height = 0;
    int layer = 0;
    for (unsigned int i=0; i<mSlabAtoms.size()-1; ++i)
    {
        if (mSlabAtoms[i].coordinates().y() <= mSlabAtoms[0].coordinates().y()+tolerance &&
                mSlabAtoms[i].coordinates().y() >= mSlabAtoms[0].coordinates().y()-tolerance &&
                mSlabAtoms[i].coordinates().z() <= mSlabAtoms[0].coordinates().z()+tolerance &&
                mSlabAtoms[i].coordinates().z() >= mSlabAtoms[0].coordinates().z()-tolerance)
        {
            ++width;
        }
        if (mSlabAtoms[i].coordinates().z() <= mSlabAtoms[0].coordinates().z()+tolerance &&
                mSlabAtoms[i].coordinates().z() >= mSlabAtoms[0].coordinates().z()-tolerance)
        {
            ++layer;
        }
        length = layer / width;
        height = mNumOfSurfAtoms / layer;
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

bool Surface::isFound(const double &inX, const double &inY, const double &inZ)
{
    const double tolerance = 0.15;
    for (unsigned int i=0; i<mSlabAtoms.size(); ++i)
    {
        if (  ((mSlabAtoms[i].coordinates().x() <= (inX + tolerance)) && (mSlabAtoms[i].coordinates().x() >= (inX - tolerance)) ) &&
                ((mSlabAtoms[i].coordinates().y() <= (inY + tolerance)) && (mSlabAtoms[i].coordinates().y() >= (inY - tolerance)) ) &&   
                ((mSlabAtoms[i].coordinates().z() <= (inZ + tolerance)) && (mSlabAtoms[i].coordinates().z() >= (inZ - tolerance)) )  )   
        {   
            return (true);
        }   
    }
    return (false);
}

int Surface::findHollow()
{
    double m_DELTA_Z = 1.5;
    if (mSurfaceType == FCC100 || mSurfaceType == BCC100)
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
                    BindingSite aSite(HOLLOW, hollowX, hollowY, hollowZ);
                    mBindingSites.push_back(aSite);
                }   
            }
        }
    }
    else if (mSurfaceType == FCC110)
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
                    BindingSite aSite(HOLLOW, hollowX, hollowY, hollowZ);
                    mBindingSites.push_back(aSite);
                }
            }
        }
    }
    else if (mSurfaceType == BCC110)
    {
        double limit = mSlabSize[0];
        double prevI = 0.5;
        for (int j=0; j<mSlabSize[1]; ++j)
        {
            for (double i=prevI; i<limit-0.5; i=i+0.5)
            {
                double offX = i * mDeltaX;
                double offY = j * mDistance + mDistance/3; // i whole number and j even OR i rational and j odd
                if (i == floor(i) && j%2 == 0) // i whole number and j even
                {
                    offY = j * mDistance + 2*mDistance/3;
                }
                else if (i != floor(i) && j%2 != 0)
                {
                    offY = j * mDistance + 2*mDistance/3; // i rational and j odd
                }
                double hollowX = mNthAtom[0] - offX;
                double hollowY = mNthAtom[1] - offY;
                double hollowZ = mNthAtom[2] + m_DELTA_Z;
                if (hollowX >= -0.05 && hollowY >= -0.05)
                {   
                    BindingSite aSite(HOLLOW, hollowX, hollowY, hollowZ);
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

int Surface::findHcp()
{
    double m_DELTA_Z = 1.5;
    double prevI_1 = 1.0;
    double limit_1 = mSlabSize[0]-1;
    double prevI_2 = 1.5;
    double limit_2 = mSlabSize[0];
    if (mSurfaceType == FCC111 || mSurfaceType == HCP0001)
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
                    BindingSite aSite(HCP, hcpX, hcpY, hcpZ);
                    mBindingSites.push_back(aSite);
                }
            }
            prevI_1 += 0.5;
            limit_1 += 0.5;
        }
    }
    else if (mSurfaceType == BCC111)
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
                        BindingSite aSite(HCP, hcpX, hcpY, hcpZ);
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
    if (mSurfaceType == BCC111)
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
                            BindingSite aSite(HCP, hcpX, hcpY, hcpZ);
                            mBindingSites.push_back(aSite);
                        }
                    }
                }
            }
        }
    }
    return (0);
} // findHcp

int Surface::findFcc()
{
    double m_DELTA_Z = 1.5;
    if (mSurfaceType == FCC111 || mSurfaceType == BCC111 || mSurfaceType == HCP0001)
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
                    BindingSite aSite(FCC, fccX, fccY, fccZ);
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
    if (mSurfaceType == BCC111)
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
                            BindingSite aSite(FCC, fccX, fccY, fccZ);
                            mBindingSites.push_back(aSite);
                        }
                    }
                }
            }
        }
    }
    return (0);
} // findFcc

int Surface::findAtop()
{
    double m_DELTA_Z = 1.5;
    for (int i=0; i<mSlabSize[0]; ++i) // i is the X offset
    {
        for (int j=0; j<mSlabSize[1]; ++j) // j is the Y offset
        {
            double offX, offY, atopX, atopY, atopZ = 0.0;
            if (mSurfaceType == FCC100 || mSurfaceType == BCC100)
            {
                offX = i * mDeltaX;
                offY = j * mDeltaX;
            }
            else if (mSurfaceType == FCC111 || mSurfaceType == BCC111 || mSurfaceType == HCP0001)
            {
                offX = (i * mDeltaX) + (j * mDeltaX/2);
                offY = j * mDistance;
            }
            else if (mSurfaceType == FCC110)
            {
                offX = i * mDeltaX;
                offY = j * mDistance;
            }
            else if (mSurfaceType == BCC110)
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
                BindingSite aSite(ATOP, atopX, atopY, atopZ);
                mBindingSites.push_back(aSite);
            }
        } // inner for (j)
    } // outer for (i)
    return (0);
} // findAtop

int Surface::findLongBridge()
{
    double m_DELTA_Z = 1.5;
    double offX, offY, LbrgX, LbrgY, LbrgZ = 0.0;
    if (mSurfaceType == FCC110)
    {
        for (int i=0; i<mSlabSize[0]; ++i) // i is the X offset
        {
            for (int j=0; j<mSlabSize[1]-1; ++j) // j is the Y offset
            {
                offX = i * mDeltaX + mDeltaX/2;
                offY = j * mDistance;

                LbrgX = mNthAtom[0] - offX;
                LbrgY = mNthAtom[1] - offY;
                LbrgZ = mNthAtom[2] + m_DELTA_Z;
                if (LbrgX >= -0.05 && LbrgY >= -0.05)
                {
                    BindingSite aSite(LONG_BRIDGE, LbrgX, LbrgY, LbrgZ);
                    mBindingSites.push_back(aSite);
                }
            }
        }
    }
    else if (mSurfaceType == BCC110)
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
                        BindingSite aSite(LONG_BRIDGE, LbrgX, LbrgY, LbrgZ);
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
                        BindingSite aSite(LONG_BRIDGE, LbrgX, LbrgY, LbrgZ);
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

int Surface::findShortBridge()
{
    double m_DELTA_Z = 1.5;
    double offX, offY, SbrgX, SbrgY, SbrgZ = 0.0;
    if (mSurfaceType == FCC110)
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
                    BindingSite aSite(SHORT_BRIDGE, SbrgX, SbrgY, SbrgZ);
                    mBindingSites.push_back(aSite);
                }
            }
        }
    }
    else if (mSurfaceType == BCC110)
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
                    BindingSite aSite(SHORT_BRIDGE, SbrgX, SbrgY, SbrgZ);
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

int Surface::findBridge()
{
    double m_DELTA_Z = 1.5;
    double offX, offY, brgX, brgY = 0.0;
    double brgZ = mNthAtom[2] + m_DELTA_Z;
    if (mSurfaceType == FCC100 || mSurfaceType == BCC100)
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
                        BindingSite aSite(BRIDGE, brgX, brgY, brgZ);
                        mBindingSites.push_back(aSite);
                    }
                }
            }
        }
    }
    else if (mSurfaceType == FCC111 || mSurfaceType == HCP0001)
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
                            BindingSite aSite(BRIDGE, brgX, brgY, brgZ);
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
                            BindingSite aSite(BRIDGE, brgX, brgY, brgZ);
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
                        BindingSite aSite(BRIDGE, brgX, brgY, brgZ);
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

std::vector<int> Surface::findNearbySites(const unsigned int atomIndex, const double radius, 
        const BINDING_SITE_TYPE siteType)
{
    std::vector<int> tempSelectedBSIndices;
    /*
    if (mSurfaceType == FCC100 || mSurfaceType == BCC100)
    {
            std::cout << "33333 Im here %%%%%%%%%%%%%%%\n";
        // TODO if (findHollow() != 0) {ERROR}
        findHollow();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == FCC111 || mSurfaceType == BCC111 || mSurfaceType == HCP0001)
    {
        findHcp();
        findFcc();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == FCC110 || mSurfaceType == BCC110)
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
    */
    // here, we have the mBindingSites vector with all the binding sites
    int numOfSites = mBindingSites.size();
    for (int i=0; i<numOfSites; ++i)
    {
            //std::cout << "1111 Im here %%%%%%%%%%%%%%%\n";
        if ( (mBindingSites[i].getType() == siteType || siteType == ALL) && 
                atomIndex <= (mSlabAtoms.size()+mAdsorbateAtoms.size()+mBindingSites.size()) )
        {
            //std::cout << "Im here %%%%%%%%%%%%%%%\n";
            //double refX = mCoordinates[atomIndex-1][0];
            //double refY = mCoordinates[atomIndex-1][1];
            double refX = mBindingSites[atomIndex - mSlabAtoms.size() - mAdsorbateAtoms.size() - 1].coordinates().x();
            double refY = mBindingSites[atomIndex - mSlabAtoms.size() - mAdsorbateAtoms.size() - 1].coordinates().y();
            double testX = mBindingSites[i].coordinates().x();
            double testY = mBindingSites[i].coordinates().y();
            /*if ((testX <= refX+radius && testX >= refX-radius) &&
                (testY <= refY+radius && testY >= refY-radius) &&
               !(testX <= refX+0.1 && testX >= refX-0.1 &&
                 testY <= refY+0.1 && testY >= refY-0.1) )*/
            // compare euclidean distance
            if ( sqrt( pow(refX-testX, 2) + pow(refY-testY, 2) ) < radius )
            {
                // the site is in the range
                //mSelectedBindingSites.push_back(mBindingSites[i]);
                tempSelectedBSIndices.push_back(i); // save index of binding sites
                //std::cout << "This site IS within the specified radius/type" << std::endl;
            }
            else
            {
                //std::cout << "ERROR: This site is NOT within the specified radius/type" << std::endl;
            }
        }
    }
    return (tempSelectedBSIndices);
}

void Surface::findAllSites()
{
    if (mSurfaceType == ANY)
    {
        std::cout << "Binding sites should be added using INPUT file." << std::endl;
        return;
    }
    if (mSurfaceType == FCC100 || mSurfaceType == BCC100)
    {
        // TODO if (findHollow() != 0) {ERROR}
        findHollow();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == FCC111 || mSurfaceType == BCC111 || mSurfaceType == HCP0001)
    {
        findHcp();
        findFcc();
        findAtop();
        findBridge();
    }
    else if (mSurfaceType == FCC110 || mSurfaceType == BCC110)
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
/*    int numOfSites = mBindingSites.size();
    for (int i=0; i<numOfSites; ++i)
    {
        mSelectedBindingSites.push_back(mBindingSites[i]);
        //std::cout << "This site IS within the specified radius/type" << std::endl;
    }
    */
}

bool Surface::writeToFile(std::string &outFile)
{
    bool success = false;
    std::ofstream ofs;
    ofs.open(outFile.c_str());

    ofs << std::to_string(mNumOfSurfAtoms+mNumOfAdsorbateAtoms+mBindingSites.size()) << "\n";
    ofs << "\n";
    ofs << std::fixed << std::setprecision(15);
    for (unsigned int i=0; i<mSlabAtoms.size(); ++i)
    {   
        ofs << mSlabAtoms[i].name() << "            ";
        ofs << mSlabAtoms[i].coordinates().x() << "            ";
        ofs << mSlabAtoms[i].coordinates().y() << "            ";
        ofs << mSlabAtoms[i].coordinates().z() << "\n";
    }
    for (unsigned int i=0; i<mAdsorbateAtoms.size(); ++i)
    {   
        ofs << mAdsorbateAtoms[i].name() << "            ";
        ofs << mAdsorbateAtoms[i].coordinates().x() << "            ";
        ofs << mAdsorbateAtoms[i].coordinates().y() << "            ";
        ofs << mAdsorbateAtoms[i].coordinates().z() << "\n";
    }
    /*for (unsigned int k=0; k<mAdsorbateCoord.size(); ++k)
    {
        ofs << mAdsorbateSymbols[k] << "            " << mAdsorbateCoord[k][0]
            << "            " << mAdsorbateCoord[k][1]
            << "            " << mAdsorbateCoord[k][2] << "\n";
    }*/
    for (unsigned int j=0; j<mBindingSites.size(); ++j)
    {
        ofs << "X             ";
        ofs << mBindingSites[j].coordinates().x() << "            ";
        ofs << mBindingSites[j].coordinates().y() << "            ";
        ofs << mBindingSites[j].coordinates().z() << "\n";
    }
    ofs.close();
    success = true;
    return (success);
}

bool Surface::writeBSToFile(std::string &outFile)
{
    bool success = false;
    std::ofstream ofs;
    ofs.open(outFile.c_str());

    ofs << std::to_string(mBindingSites.size()) << "\n";
    ofs << "\n";
    ofs << std::fixed << std::setprecision(15);
    for (unsigned int j=0; j<mBindingSites.size(); ++j)
    {
        ofs << "X             ";
        ofs << mBindingSites[j].coordinates().x() << "            ";
        ofs << mBindingSites[j].coordinates().y() << "            ";
        ofs << mBindingSites[j].coordinates().z() << "\n";
    }
    ofs.close();
    success = true;
    return (success);
}




void Surface::resetGeometry()
{
//    mSurfaceSymbols.clear();
    std::vector<BindingSite> swap1;
    mBindingSites.swap(swap1);
    //mSelectedBindingSites.swap(swap1);
    /*std::vector<std::string> swap2;
    mSurfaceSymbols.swap(swap2);
    mAdsorbateSymbols.swap(swap2);
    std::vector< std::vector<double> > swap3;
    mCoordinates.swap(swap3);
    mAdsorbateCoord.swap(swap3);*/
}

const SLAB_TYPE stringToSlabType(std::string in)
{
    if (in == "fcc100")
        return FCC100;
    else if (in == "fcc110")
        return FCC110;
    else if (in == "fcc111")
        return FCC111;
    else if (in == "bcc100")
        return BCC100;
    else if (in == "bcc110")
        return BCC110;
    else if (in == "bcc111")
        return BCC111;
    else if (in == "hcp0001")
        return HCP0001;
    else if (in == "ANY")
        return ANY;
    else
        throw std::exception();
}

void Surface::addBindingSites(BindingSite BS1, BindingSite BS2, BindingSite BS3)
{
    mBindingSites.push_back(BS1);
    mBindingSites.push_back(BS2);
    mBindingSites.push_back(BS3);
}
