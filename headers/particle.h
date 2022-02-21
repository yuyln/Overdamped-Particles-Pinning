#ifndef __PARTICLE
#define __PARTICLE
#include <cmath>

class Particle
{
public:
    double beta, damp, U0;
    double x, y;

    Particle(double betadamp_, double U0_, double x_, double y_);
    Particle InitParticle(double betadamp_, double U0_, double x_, double y_);
    double PotentialParticle(const Particle &cur, const Particle &other);
    double PotentialParticle(const Particle &other);
};

Particle::Particle(double betadamp_, double U0_, double x_, double y_)
{
    damp = 1.0 / sqrt(betadamp_ * betadamp_ + 1.0);
    beta = betadamp_ * damp;
    U0 = U0_;
    x = x_;
    y = y_;
}

Particle Particle::InitParticle(double betadamp_, double U0_, double x_, double y_)
{
    return Particle(betadamp_, U0_, x_, y_);
}

#endif