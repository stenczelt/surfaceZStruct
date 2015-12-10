// Mina Jafari
// 12-08-2015

// atom assignment: first layer
//
//  ---(N-1)-----(N)
//       /        /
//      /        /
//  -(N-1)*----(N)*
//    /        /

#ifndef _BCC111_H_
#define _BCC111_H_
#include "Termination111.h"
#include <vector>
#include <string>

class bcc111:public Termination111
{
    public:
        void findShallowBridge(const unsigned int offsetX, const unsigned int offsetY);
        void findDeepBridge(const unsigned int offsetX, const unsigned int offsetY);
};
#endif
