#ifndef __POINT
#define __POINT
class Point
{
public:
    double x, y;
    Point(): x(0.0), y(0.0) {}
    Point(double x_, double y_): x(x_), y(y_) {}
    const double &X() const noexcept { return x; }
    const double &Y() const noexcept { return y; }
};

typedef Point Vector;
#endif