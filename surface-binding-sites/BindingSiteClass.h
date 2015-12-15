// Mina Jafari
// 12-14-2015

#ifndef _BINDINGSITECLASS_H_
#define _BINDINGSITECLASS_H_
//#include <iostream>
#include <string>

class BindingSiteClass
{
    private:
        // hollow, hcp, fcc, atop, long bridge, short bridge, bridge
        std::string mType = "";
        int mCoordinate[3] = {0, 0 , 0};

    public:
        BindingSiteClass();
        BindingSiteClass(std::string inType, double inX, double inY, double ins);
        std::string getType() const;
        int getX() const;
        int getY() const;
        int getZ() const;
        void setType(std::string inType);
        void setCoordinates(int inX, int inY, int inZ);
        //overload insersion << operator to print the binding site vector

};
#endif
