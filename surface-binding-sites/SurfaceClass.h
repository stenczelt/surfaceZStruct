// Mina Jafari
// 12-14-2015
// This class works alongside the BindingSiteClass. This class defines a slab with
// its size and finds all the binding sites on that surface. The binding sites within a
// radius can be found using the findNearbySites() function. To use this class, make
// an object of type SurfaceClass and call setAtoms, setSurfaceType functions on it.
//
// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

#ifndef _SURFACECLASS_H_
#define _SURFACECLASS_H_
#include "BindingSiteClass.h"
#include <string>
#include <vector>

class SurfaceClass
{
    private:
        // fcc100, 110, 111, bcc100, 110, 111, hcp0001
        std::string mSurfaceType;
        // number of surface atoms
        int mNumOfSurfAtoms = 0;
        int mNumOfAdsorbateAtoms = 0;
        // number of atoms in x, y, and z direction in slab
        int mSlabSize[3] = {0, 0, 0}; // x * y * z
        // a vector to store all the binding sites of the surface
        std::vector<BindingSiteClass> mBindingSites; 
        // a vector to store only the sites within a range and atom
        std::vector<BindingSiteClass> mSelectedBindingSites; 
        // a vector to store atomic sybmols
        std::vector<std::string> mSurfaceSymbols;
        std::vector<std::string> mAdsorbateSymbols;
        // a 2D vector to store atomic coordinates
        std::vector< std::vector<double> > mCoordinates;
        std::vector< std::vector<double> > mAdsorbateCoord;
        // some parameters
        double mDeltaX = 0.0;
        double mDeltaY = 0.0;
        double mDistance = 0.0;
        // x, y, z of some important atoms used to find binding sites
        double mNthAtom [3];
        double mNthMinusOneAtom [3];
        double mStarAtom [3];
        double mStarMinusOneAtom [3];
        // 2nd and 3rd layer z components
        double mSecondLayerZ = 0.0;
        double mThirdLayerZ = 0.0;
        void resetGeometry();

    public:
        // getter functions
        std::string getSurfaceType() const;
        int getNumOfAtoms() const;
        int getSurfaceWidth() const;
        int getSurfaceLength() const;
        int getSurfaceHeight() const;
        const BindingSiteClass* getBindingSite(unsigned int index) const; // zero indexed
        //setter functions
        bool setSurfaceType(std::string inSurface);
        bool setAtoms(int numOfAtoms, double* coordinates, std::string* atomicSymbols);
        bool setSlabSize();
        // other member functions
        bool isFound(const double &inX, const double &inY, const double &inZ);
        bool writeToFile();
        int findHollow();
        int findHcp();
        int findFcc();
        int findAtop();
        int findLongBridge();
        int findShortBridge();
        int findBridge(); //TODO shallow and deep bridge(bcc111)
        // finds sites within the given radius of the specified atom(atomIndex)
        // and the binding type. Atom indexing starts from one.
        void findNearbySites(const int atomIndex, const double range, 
                             const std::string siteType);
        void findAllSites();
        bool writeToFile(std::string &outFile);
};
#endif
