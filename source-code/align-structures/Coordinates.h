#ifndef COORDINAETS_H
#define COORDINAETS_H
class Coordinates
{
    private:
        double mX;
        double mY;
        double mZ;

    public:
        Coordinates(double x, double y, double z);
        ~Coordinates(){}
        //TODO learn about inlining
        inline const double x() const { return mX; }
        inline const double y() const { return mY; }
        inline const double z() const { return mZ; }
};



#endif
