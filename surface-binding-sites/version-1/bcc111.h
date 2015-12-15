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

// this is inherited from Termination111 class
class bcc111:public Termination111
{
    private:
        // this method does not apply to bcc111, this function is replaced with the 
        // following public functions
        void findBridge(const unsigned int offsetX, const unsigned int offsetY) {};
    public:
        // find binding sites that do not exist for other 111 terminations 
        void findShallowBridge(const unsigned int offsetX, const unsigned int offsetY);
        void findDeepBridge(const unsigned int offsetX, const unsigned int offsetY);
};
#endif
