// Mina Jafari
// 12-14-2015

#ifndef _BINDINGSITECLASS_H_
#define _BINDINGSITECLASS_H_
#include <string>

class BindingSiteClass
{
    private:
        // hollow, hcp, fcc, atop, long bridge, short bridge, bridge
        std::string mType = "";
        double mCoordinate[3] = {0.0, 0.0 , 0.0};

    public:
        BindingSiteClass();
        BindingSiteClass(std::string inType, double inX, double inY, double ins);
        std::string getType() const;
        double getX() const;
        double getY() const;
        double getZ() const;
        void setType(std::string inType);
        void setCoordinates(double inX, double inY, double inZ);
};
#endif
