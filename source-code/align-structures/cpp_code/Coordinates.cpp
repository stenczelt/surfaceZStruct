#include "Coordinates.h"


Coordinates::Coordinates()
{
    mX = mY = mZ = 0.0;
}

Coordinates::Coordinates(
        double x, 
        double y, 
        double z)
: 
    mX(x), 
    mY(y), 
    mZ(z)
{
}

