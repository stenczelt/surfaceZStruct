// Mina Jafari
// 12-03-2015

// atom assignment
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

#ifndef _FCC100_H_
#define _FCC100_H_
#include <vector>
#include <string>

class fcc100
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

    public:
        bool setAtoms(const std::vector<std::string> &xyzFile);
        void findHollow();
        void findAtop();
        void findBridge();
};
#endif
