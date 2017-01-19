// Mina Jafrai
// 12-14-2015

#include "BindingSite.h"
#include <string>

/*BindingSite::BindingSite()
{
}*/

BindingSite::BindingSite(BINDING_SITE_TYPE inType, double inX, double inY, double inZ):
    mType(inType),
    mCoordinates(inX, inY, inZ),
    mAdsorbate(NULL)
{
}

BINDING_SITE_TYPE BindingSite::getType() const
{
    return (mType);
}

