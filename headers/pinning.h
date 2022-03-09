#ifndef __PINNING
#define __PINNING
#include <cmath>
#include <table.h>

class Pinning
{
public:
    double x, y, U0, R0, F0;
    double R02;
    Pinning() {}
    Pinning(double x_, double y_, double U0_, double R0_) : x(x_), y(y_), U0(U0_), R0(R0_)
    {
        R02 = R0 * R0;
        F0 = 2.0 * U0 / R02;
    }

    static double Potential(const Pinning *p, const double &x, const double &y, const Table &table) noexcept
    {
        double dx = x - p->x;
        double dy = y - p->y;
        double d2 = dx * dx + dy * dy;
        return p->U0 * table(d2 / p->R02);
    }

    static void Force(const Pinning *p, const double &x, const double &y, const Table &table, double *fx, double *fy) noexcept
    {
        double dx = x - p->x;
        double dy = y - p->y;
        double d2 = dx * dx + dy * dy;
        if (d2 / p->R02 > table.getMaxRange())
        {
            *fx = 0.0;
            *fy = 0.0;
            return;
        }
        double e = table(d2 / p->R02);
        *fx = p->F0 * e * dx;
        *fy = p->F0 * e * dy;
    }
};

#endif