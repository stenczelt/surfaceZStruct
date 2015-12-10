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

class diamond100
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
        double mSecLayerZ = 0.0;
        double mThirLayerZ = 0.0;

    public:
        bool setAtoms(const std::vector<std::string> &xyzFile);
        bool isFound(const double &inX, const double &inY, const double &inZ);
        void findAtop();
        void findFirstBridge();
        void findSecondBridge();
        void findThirdBridge();
};
#endif
