// Mina Jafari
// 12-09-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

#ifndef _DIAMOND100_H_
#define _DIAMOND100_H_
#include <vector>
#include <string>

// the class for diamond100 surface termination
class diamond100
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
        // delta X & Y between the last two atoms
        double mDeltaX;
        double mDeltaY;
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
        void findAtop();
        void findFirstBridge();
        void findSecondBridge();
        void findThirdBridge();
};
#endif
