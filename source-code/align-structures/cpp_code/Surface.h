// Mina Jafari
// 12-14-2015
// This class works alongside the BindingSite. This class defines a slab with
// its size and finds all the binding sites on that surface. The binding sites within a
// radius can be found using the findNearbySites() function. To use this class, make
// an object of type Surface and call setAtoms, setSurfaceType functions on it.
//
// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

// A Surface object has 3 main attributes:
//   - Atoms
//   - Adsorbates (max 2 for now) 
//   - Binding Sites

#ifndef _SURFACECLASS_H_
#define _SURFACECLASS_H_
#include "BindingSite.h"
#include <string>
#include <vector>

enum SLAB_TYPE
{
    FCC100  = 0,
    FCC110  = 1,
    FCC111  = 2,
    BCC100  = 3,
    BCC110  = 4,
    BCC111  = 5,
    HCP0001 = 6,
    ANY     = 7,
};

class Surface
{
    private:
        // fcc100, 110, 111, bcc100, 110, 111, hcp0001
        SLAB_TYPE mSurfaceType;
        // number of surface atoms
        int mNumOfSurfAtoms = 0; //TODO
        int mNumOfAdsorbateAtoms = 0; //TODO
        // number of atoms in x, y, and z direction in slab
        int mSlabSize[3] = {0, 0, 0}; // x * y * z
        // a vector to store all the binding sites of the surface
        std::vector<BindingSite> mBindingSites; 
        std::vector<Atom> mBindingSitesTemp;
        // a vector to store only the sites within a range and atom
        // a vector to store atomic sybmols
        /*
        std::vector<std::string> mSurfaceSymbols;
        std::vector<std::string> mAdsorbateSymbols;
        */

        // Vector of slab atoms (Cu's for instance)
        std::vector<Atom> mSlabAtoms;
        // Vector of adsorbates (NH3 for instance), set to 2 for now
        std::vector<Molecule> mAdsorbates;
        std::vector<Atom> mAdsorbateAtoms;
        

        // a 2D vector to store atomic coordinates
        //std::vector< std::vector<double> > mCoordinates;
        //std::vector< std::vector<double> > mAdsorbateCoord;
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
        Surface();
        // temp constructor
        Surface(std::vector<Atom> slabAtoms, std::vector<Atom> bindingSites,
                std::vector<Molecule> adsorbateList);
        /*Surface(std::vector<Atom> slabAtoms, std::vector<BindingSite> bindingSites,
                std::vector<Molecule> adsorbateList);*/
        // getter functions
        SLAB_TYPE getSurfaceType() const;
        int getNumOfAtoms() const;
        int getSurfaceWidth() const;
        int getSurfaceLength() const;
        int getSurfaceHeight() const;
        const BindingSite* getBindingSite(unsigned int index) const; // zero indexed
        //setter functions
        void setSurfaceType(SLAB_TYPE inSurface);
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
        std::vector<int> findNearbySites(const unsigned int atomIndex, const double range, 
                const BINDING_SITE_TYPE siteType);
        void findAllSites();
        bool writeToFile(std::string &outFile);
        bool writeBSToFile(std::string &outFile);
        void addBindingSites(BindingSite BS1, BindingSite BS2, BindingSite BS3);
};

const SLAB_TYPE stringToSlabType(std::string in);
#endif
