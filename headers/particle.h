#ifndef __PARTICLE
#define __PARTICLE
#include <cmath>
#include <bessel0table.h>
#include <bessel1table.h>


typedef struct Particle
{
    double beta, damp, U0;
    double x, y;
} Particle;

Particle InitParticle(double betadamp_, double U0_, double x_, double y_)
{
    Particle p;
    p.damp = 1.0 / sqrt(betadamp_ * betadamp_ + 1.0);
    p.beta = betadamp_ * p.damp;
    p.U0 = U0_;
    p.x = x_;
    p.y = y_;
    return p;
}

double PotentialParticle(const Particle *cur, const double x, const double y)
{
    double dx, dy, d;
    dx = cur->x - x;
    dy = cur->y - y;
    if (dx > Bess0T.params.TCut || dy > Bess0T.params.TCut)
    {
        return 0.0;
    }

    if (dx == 0.0 && dy == 0.0)
    {
        return 0.0;
    }

    d = sqrt(dx * dx + dy * dy);
    return cur->U0 * GetBess0Table(d);
}

void ForceParticle(const Particle *interact, const double x, double y, 
                   double *fx, double *fy)
{
    double dx, dy, d;
    dx = interact->x - x;
    dy = interact->y - y;
    if (dx > Bess1T.params.TCut || dy > Bess1T.params.TCut)
    {
        *fx = 0.0;
        *fy = 0.0;
        return;
    }

    if (dx == 0.0 && dy == 0.0)
    {
        *fx = 0.0;
        *fy = 0.0;
        return;
    }

    d = sqrt(dx * dx + dy * dy);
    double d1 = 1.0 / d;
    double Bes = interact->U0 * GetBess1Table(d);
    *fx = Bes * dx * d1;
    *fy = Bes * dy * d1;
}
#endif