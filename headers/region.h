#ifndef __REGION
#define __REGION

#include <point.h>
/*
 (p2)-------(p3)
 |           |
 |           |
(p0)--------(p1)
*/

class Rectangle
{
    Point left_below;
    double w, h;
public:

    Rectangle() {}
    Rectangle(Point l, double width, double heigth): left_below(l), w(width), h(heigth) {}
    Rectangle(double x, double y, double width, double heigth): left_below(x, y), w(width), h(heigth) {}
    const Point &LB() const noexcept { return left_below; }
    const double &W() const noexcept { return w; }
    const double &H() const noexcept { return h; }
    void operator=(const Rectangle &o) { memcpy((void*)this, (void*)&o, sizeof(Rectangle)); }

    bool Inside(const double &x, const double &y) const noexcept { return x > left_below.X() &&
                                                           x < (left_below.X() + w) &&
                                                           y > left_below.Y() &&
                                                           y < (left_below.Y() + h); }
    // Point ClosestPoint(const double &x, const double &y) const noexcept
    // {
    //     return Point(0.0, 0.0);
    // }
};

class Circle
{
    Point center;
    double r;
public:
    Circle() {}
    Circle(Point c, double radius): center(c), r(radius) {}
    Circle(double x, double y, double radius): center(x, y), r(radius) {}
    const Point &C() const noexcept { return center; }
    const double &R() const noexcept { return r; }
    void operator=(const Circle &o) { memcpy((void*)this, (void*)&o, sizeof(Circle)); }
    bool Inside(const double &x, const double &y) const noexcept { return ((x - center.X()) * (x - center.X()) +
                                                                           (y - center.Y()) * (y - center.Y())) < 
                                                                           (r * r); }

    // Point ClosestPoint(const double &x, const double &y) const noexcept
    // {
    //     return Point(center.X(), center.Y());
    // }
};

class Triangule
{
    Point p1, p2, p3;
public:
    Triangule() {}
    Triangule(Point _1, Point _2, Point _3): p1(_1), p2(_2), p3(_3) {}
    const Point &P1() const noexcept { return p1; }
    const Point &P2() const noexcept { return p2; }
    const Point &P3() const noexcept { return p3; }
    void operator=(const Triangule &o) { memcpy((void*)this, (void*)&o, sizeof(Triangule)); }
    bool Inside(const double &x, const double &y)
    {
        double alpha = ((p2.y - p3.y)*(x - p3.x) + (p3.x - p2.x)*(y - p3.y)) /
                ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y));
        double beta = ((p3.y - p1.y)*(x - p3.x) + (p1.x - p3.x)*(y - p3.y)) /
            ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y));
        double gamma = 1.0 - alpha - beta;
        return alpha > 0.0 && beta > 0.0 && gamma > 0.0;
    }

};


#endif