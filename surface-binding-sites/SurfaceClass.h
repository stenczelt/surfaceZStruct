// Mina Jafari
// 12-14-2015

#ifndef _SURFACECLASS_H_
#define _SURFACECLASS_H_
#include "BindingSiteClass.h"
#include <string>
#include <vector>

class SurfaceClass
{
    private:
        static double m_DELTA_Z;
        // fcc100, 110, 111, bcc100, 110, 111, hcp0001
        std::string mSurfaceType = "";
        int mNumOfAtoms = 0;
        int mSlabSize[3] = {0, 0, 0}; // x * y * z
        std::vector<BindingSiteClass> mBindingSites; 
        std::vector<BindingSiteClass> mSelectedBindingSites; 
        std::vector<std::string> mAtomicSymbols;
        // a 2D vector to store atomic coordinates
        std::vector< std::vector<double> > mCoordinates;
        double mDeltaX = 0.0;
        double mDeltaY = 0.0;
        double mDistance = 0.0;
        double mNthAtom[3] = {0, 0, 0};
        double mNthMinusOneAtom [3] = {0, 0, 0};
        double mStarAtom [3] = {0, 0, 0};
        double mStarMinusOneAtom [3] = {0, 0, 0};
        double mSecondLayerZ = 0.0;
        double mThirdLayerZ = 0.0;

    public:
        std::string getSurfaceType() const;
        int getNumOfAtoms() const;
        int getSurfaceWidth() const;
        int getSurfaceLength() const;
        int getSurfaceHeight() const;
        BindingSiteClass getBindingSite(unsigned int element) const;
        bool setSurfaceType(std::string inSurface);
        bool setAtoms(int numOfAtoms, double* coordinates, std::string* atomicSymbols);
        //void setSlabSize(const std::vector< std::vector<double> > &coordinates);
        void setSlabSize();
        bool isFound(const double &inX, const double &inY, const double &inZ);
        bool writeToFile();
        void findHollow();
        void findHcp();
        void findFcc();
        void findAtop();
        void findLongBridge();
        void findShortBridge();
        void findBridge(); //TODO shallow and deep bridge(bcc111)
        void findNearbySites(int atomIndex, double range, std::string siteType);
        bool writeToFile(std::string &outFile);
};
#endif
