// Mina Jafari
// 12-14-2015
// This class is a helper class for SurfaceClass. It defines binding sites of
// a surface using their type and xyz coordinates.

#ifndef _BINDINGSITE_H_
#define _BINDINGSITE_H_
#include <string>
#include <Coordinates.h>
#include <Molecule.h>

//TODO learn about namespaces
enum BINDING_SITE_TYPE
{
    ALL          = 0,
    HOLLOW       = 1,
    HCP          = 2,
    FCC          = 3,
    ATOP         = 4,
    LONG_BRIDGE  = 5,
    SHORT_BRIDGE = 6,
    BRIDGE       = 7,
};


class BindingSite
{
    private:
        // hollow, hcp, fcc, atop, long-bridge, short-bridge, bridge
        BINDING_SITE_TYPE mType;
        // coordinates of the site
        Coordinates mCoordinates;
        Molecule* mAdsorbate;

    public:
        // Ctors
        BindingSite(BINDING_SITE_TYPE inType, double inX, double inY, double ins);
        // getter functions
        BINDING_SITE_TYPE getType() const;
        inline const Coordinates& coordinates() const { return mCoordinates; }
        
        // setter functions
        
        int connectToAdsorbate(Molecule* adsorbate);
        //int disconnectFromAdsorbate() {mAdsorbate = NULL;}
};

#endif
