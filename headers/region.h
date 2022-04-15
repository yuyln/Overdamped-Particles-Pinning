#ifndef __REGION
#define __REGION

#include <point.h>
/*
 (p2)-------(p3)
 |           |
 |           |
(p0)--------(p1)
*/

class Region
{
    Point p0, p1, p2, p3;

    Region() {}
    Region(Point p0_, Point p1_, Point p2_, Point p3_):
    p0(p0_), p1(p1_), p2(p2_), p3(p3_) {}
    Region(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3):
    p0(x0, y0), p1(x1, y1), p2(x2, y2), p3(x3, y3) {}

    bool Inside(const double &x, const double &y)
    {
        return x > p0.X() && x < p1.X() && y > p0.Y() && y < p2.Y();
    }
};

#endif