#ifndef COORDINAETS_H
#define COORDINAETS_H
class Coordinates
{
    private:
        double mX;
        double mY;
        double mZ;

    public:
        Coordinates();
        Coordinates(double x, double y, double z);
        ~Coordinates(){}
        //TODO learn about inlining
        inline double x() const { return mX; }
        inline double y() const { return mY; }
        inline double z() const { return mZ; }
};



#endif
