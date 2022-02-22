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

Particle InitParticle(double betadamp_, double U0_, double x_, double y_);

double PotentialParticle(const Particle *cur, const Particle *other);
void ForceParticle(const Particle *cur, const Particle *other, double &fx, double &fy);


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

double PotentialParticle(const Particle *cur, const Particle *other)
{
    double dx, dy, d;
    dx = cur->x - other->x;
    dy = cur->y - other->y;
    if (dx > Bess0T.TCut || dy > Bess0T.TCut)
    {
        return 0.0;
    }

    if (dx == 0.0 && dy == 0.0)
    {
        return 0.0;
    }

    d = sqrt(dx * dx + dy * dy);
    return cur->U0 * other->U0 * GetBess0Table(d);
}

void ForceParticle(const Particle *cur, const Particle *other, double &fx, double &fy)
{
    double dx, dy, d;
    dx = cur->x - other->x;
    dy = cur->y - other->y;
    if (dx > Bess1T.TCut || dy > Bess1T.TCut)
    {
        fx = 0.0;
        fy = 0.0;
        return;
    }

    if (dx == 0.0 && dy == 0.0)
    {
        fx = 0.0;
        fy = 0.0;
        return;
    }

    d = sqrt(dx * dx + dy * dy);
    double d1 = 1.0 / d;
    double Bes = cur->U0 * other->U0 * GetBess1Table(d);
    fx = Bes * dx * d1;
    fy = Bes * dy * d1;
}
#endif