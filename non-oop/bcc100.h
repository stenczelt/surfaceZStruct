// Mina Jafari
// 12-08-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       |        |
//       |        |
//  ---(N-1)*----(N)*
//       |        |

#ifndef _BCC100_H_
#define _BCC100_H_
#include <vector>
#include <string>

class bcc100
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
