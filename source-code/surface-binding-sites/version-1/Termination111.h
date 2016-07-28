// Mina Jafari
// 12-09-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       /        /
//      /        /
//  -(N-1)*----(N)*
//    /        /

#ifndef _TERMINATION111_H_
#define _TERMINATION111_H_
#include <vector>
#include <string>

// base class for fcc111, bcc111 and hcp0001 terminations
class Termination111
{
    protected:
        // displacement from surface
        static double m_DELTA_Z;
        // vector to store input file
        std::vector<std::string> mVector;
        // atom assignment based on the top scheme
        double mNthAtom[3];
        double mNthMinusOneAtom [3];
        double mStarAtom [3];
        double mStarMinusOneAtom [3];
        // delta X & Y between the last two atoms
        double mDeltaX;
        double mDeltaY;
        // vertical distance between (N) & (N)*
        double mDistance;
        // Z value of second layer
        double mSecLayerZ = 0.0;
        // Z value of third layer
        double mThirLayerZ = 0.0;

    public:
        // sets the attributes read from an input file
        bool setAtoms(const std::vector<std::string> &xyzFile);
        // searches the input file to confirm the existence of calculated coordinates
        bool isFound(const double &inX, const double &inY, const double &inZ);
        // find binding sites
        void findAtop(const unsigned int offsetX, const unsigned int offsetY);
        void findBridge(const unsigned int offsetX, const unsigned int offsetY);
        void findFcc(const unsigned int offsetX, const unsigned int offsetY);
        void findHcp(const unsigned int offsetX, const unsigned int offsetY);
};
#endif
