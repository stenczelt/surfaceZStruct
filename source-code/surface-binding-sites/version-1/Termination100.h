// Mina Jafari
// 12-09-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

#ifndef _TERMINATION100_H_
#define _TERMINATION100_H_
#include <vector>
#include <string>

// base class for fcc100, bcc100
class Termination100
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

    public:
        // sets the attributes read from an input file
        bool setAtoms(const std::vector<std::string> &xyzFile);
        // TODO Is the isFound function needed here?
        // find binding sites 
        void findHollow(const unsigned int offsetX, const unsigned int offsetY);
        void findAtop(const unsigned int offsetX, const unsigned int offsetY);
        void findBridge(const unsigned int offsetX, const unsigned int offsetY);
};
#endif
