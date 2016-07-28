// Mina Jafari
// 12-08-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       /        /
//      /        /
//  -(N-1)*----(N)*
//    /        /

#ifndef _BCC110_H_
#define _BCC110_H_
#include <vector>
#include <string>

// class for bcc110 termiantion
class bcc110
{
    private:
        // displacement from surface
        static double m_DELTA_Z;
        // vector to store input file
        std::vector<std::string> mVector;
        // atom assignment based on the top scheme
        double mNthAtom[3];
        double mNthMinusOneAtom [3];
        double mStarAtom [3];
        double mStarMinusOneAtom [3];
        // vertical distance between (N) & (N)*
        double mDistance;
        // delta X & Y between the last two atoms
        double mDeltaX;
        double mDeltaY;

    public:
        // sets the attributes read from an input file 
        bool setAtoms(const std::vector<std::string> &xyzFile);
        // searches the input file to confirm the existence of calculated coordinates
        bool isFound(const double &inX, const double &inY, const double &inZ);
        // find binding sites
        void findHollow(const unsigned int offsetX, const unsigned int offsetY);
        void findAtop(const unsigned int offsetX, const unsigned int offsetY);
        void findLongBridge(const unsigned int offsetX, const unsigned int offsetY);
        void findShortBridge(const unsigned int offsetX, const unsigned int offsetY);
};
#endif
