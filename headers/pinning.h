#ifndef __PINNING
#define __PINNING
#ifndef OPENCLCOMP
#include <cmath>
#endif
#include <exptable.h>

typedef struct Pinning
{
    double U0, R0;
    double x, y;
} Pinning;

Pinning InitPinning(double x, double y, double U0, double R0)
{
    Pinning p;
    p.x = x;
    p.y = y;
    p.U0 = U0;
    p.R0 = R0;
    return p;
}


#endif