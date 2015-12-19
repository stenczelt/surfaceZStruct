// Mina Jafari
// 12-14-2015
// This class is a helper class for SurfaceClass. It defines binding sites of
// a surface using their type and xyz coordinates.

#ifndef _BINDINGSITECLASS_H_
#define _BINDINGSITECLASS_H_
#include <string>

class BindingSiteClass
{
    private:
        // hollow, hcp, fcc, atop, long-bridge, short-bridge, bridge
        std::string mType = "";
        // coordinates of the site
        double mCoordinate[3] = {0.0, 0.0 , 0.0};

    public:
        // Ctors
        BindingSiteClass();
        BindingSiteClass(std::string inType, double inX, double inY, double ins);
        // getter functions
        std::string getType() const;
        double getX() const;
        double getY() const;
        double getZ() const;
        // setter functions
        void setType(std::string inType);
        void setCoordinates(double inX, double inY, double inZ);
};
#endif
