// Mina Jafari
// 12-04-2015

// atom assignment
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

#ifndef _FCC110_H_
#define _FCC110_H_
#include <vector>
#include <string>

class fcc110
{
    private:
        static double m_DELTA_Z;
        std::vector<std::string> mVector;
        double mNthAtom[3];
        double mNthMinusOneAtom [3];
        double mStarAtom [3];
        double mStarMinusOneAtom [3];
        double mDistance;
        double mDeltaX;
        double mDeltaY;

    public:
//        double calcConsDis(const std::vector<std::string> &xyzFile);
        bool setAtoms(const std::vector<std::string> &xyzFile);
        bool isFound(const double &inX, const double &inY, const double &inZ);
        void findHollow();
        void findAtop();
        void findLongBridge();
        void findShortBridge();
//        bool setVector(const std::vector<std::string> &inVector);
};
#endif
