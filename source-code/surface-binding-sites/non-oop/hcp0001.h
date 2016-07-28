// Mina Jafari
// 12-08-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       /        /
//      /        /
//  -(N-1)*----(N)*
//    /        /

#ifndef _HCP0001_H_
#define _HCP0001_H_
#include <vector>
#include <string>

class hcp0001
{
    private:
        static double m_DELTA_Z;
        std::vector<std::string> mVector;
        double mNthAtom[3];
        double mNthMinusOneAtom [3];
        double mStarAtom [3];
        double mStarMinusOneAtom [3];
        double mDeltaX;
        double mDeltaY;
        double mDistance;
        double mSecLayerZ = 0.0;
        double mThirLayerZ = 0.0;

    public:
        bool setAtoms(const std::vector<std::string> &xyzFile);
        bool isFound(const double &inX, const double &inY, const double &inZ);
        void findAtop();
        void findBridge();
        void findFcc();
        void findHcp();
};
#endif
