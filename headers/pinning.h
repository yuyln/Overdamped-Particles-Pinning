#ifndef __PINNING
#define __PINNING
#ifndef OPENCLCOMP
#include <cmath>
#endif
#include <exptable.h>

typedef struct Pinning
{
    double U0, R0, F0;
    double R02;
    double x, y;
} Pinning;

Pinning InitPinning(double x, double y, double U0, double R0)
{
    Pinning p;
    p.x = x;
    p.y = y;
    p.U0 = U0;
    p.R0 = R0;
    p.R02 = R0 * R0;
    p.F0 = 2.0 * U0 / p.R02;
    return p;
}

double PotentialPinning(const Pinning *p, const double &x, const double &y)
{
    double dx = x - p->x;
    double dy = y - p->y;
    double d2 = dx * dx + dy * dy;
    return p->U0 * GetExpTable(d2 / p->R02);
}

void ForcePinning(const Pinning *p, const double &x, const double &y, double &fx, double &fy)
{
    double dx = x - p->x;
    double dy = y - p->y;
    double d2 = dx * dx + dy * dy;
    if (d2 / p->R02 > ExpT.params.TLimit)
    {
        fx = 0.0;
        fy = 0.0;
        return;
    }
    double e = GetExpTable(d2 / p->R02);
    fx = p->F0 * e * dx;
    fy = p->F0 * e * dy;
}


#endif