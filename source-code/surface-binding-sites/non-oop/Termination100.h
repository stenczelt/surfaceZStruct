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

class Termination100
{
    protected:
        static double m_DELTA_Z;
        std::vector<std::string> mVector;
        double mNthAtom[3];
        double mNthMinusOneAtom [3];
        double mStarAtom [3];
        double mStarMinusOneAtom [3];
        double mDeltaX;
        double mDeltaY;

    public:
        bool setAtoms(const std::vector<std::string> &xyzFile);
        void findHollow(const int offsetX, const int offsetY);
        void findAtop(const int offsetX, const int offsetY);
        void findBridge(const int offsetX, const int offsetY);
};
#endif
