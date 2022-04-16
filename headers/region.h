#ifndef __REGION
#define __REGION

#include <point.h>
#include <cstring>
#include <functions.h>
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

int ReadRectangles(Rectangle **p)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/regions/rectangles", &nfiles, &qnt, &parsed);
    int nP = 0;

    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Data", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 4)
        {
            nP++;
        }
    }
    *p = new Rectangle[nP];

    int j_ = 0;
    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Data", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 4)
        {
            double x = strtod(parsed[i][j], NULL);
            double y = strtod(parsed[i][j + 1], NULL);
            double w = strtod(parsed[i][j + 2], NULL);
            double h = strtod(parsed[i][j + 3], NULL);
            (*p)[j_] = Rectangle(x, y, w, h);
            j_++;
        }
    }
    return nP;
}

int ReadTriangules(Triangule **p)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/regions/triangules", &nfiles, &qnt, &parsed);
    int nP = 0;

    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Vertex", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 6)
        {
            nP++;
        }
    }
    *p = new Triangule[nP];

    int j_ = 0;
    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Vertex", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 6)
        {
            double x0 = strtod(parsed[i][j], NULL);
            double y0 = strtod(parsed[i][j + 1], NULL);
            double x1 = strtod(parsed[i][j + 2], NULL);
            double y1 = strtod(parsed[i][j + 3], NULL);
            double x2 = strtod(parsed[i][j + 4], NULL);
            double y2 = strtod(parsed[i][j + 5], NULL);
            Point p0(x0, y0), p1(x1, y1), p2(x2, y2);
            (*p)[j_] = Triangule(p0, p1, p2);
            j_++;
        }
    }
    return nP;
}

int ReadCircles(Circle **p)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/regions/circles", &nfiles, &qnt, &parsed);
    int nP = 0;

    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Centers", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 2)
        {
            nP++;
        }
    }
    *p = new Circle[nP];

    int j_ = 0;
    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Centers", parsed[i], qnt[i]);
        int start = 0;
        double R0 = GetValueDouble("Radius", parsed[i], qnt[i]);
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 2)
        {
            double x = strtod(parsed[i][j], NULL);
            double y = strtod(parsed[i][j + 1], NULL);
            (*p)[j_] = Circle(x, y, R0);
            j_++;
        }
    }
    return nP;
}

#endif