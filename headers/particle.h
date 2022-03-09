#ifndef __PARTICLE
#define __PARTICLE
#include <cmath>
#include <table.h>

class Particle
{
public:
    double x, y, beta, damp, U0, Vx, Vy;
    Particle() {}
    Particle(double betadamp_, double U0_, double x_, double y_) : x(x_), y(y_), U0(U0_)
    {
        damp = 1.0 / sqrt(betadamp_ * betadamp_ + 1.0);
        beta = betadamp_ * damp;
    }

    static double Potential(const Particle *cur, const double &x, const double &y, const Table &table) noexcept
    {
        double dx, dy, d;
        dx = cur->x - x;
        dy = cur->y - y;
        if (dx > table.getMaxRange() || dy > table.getMaxRange())
        {
            return 0.0;
        }

        if (dx == 0.0 && dy == 0.0)
        {
            return 0.0;
        }

        d = sqrt(dx * dx + dy * dy);
        return cur->U0 * table(d);
    }

    static void Force(const Particle *interact, const double &x, const double &y, const Table &table,
                       double *fx, double *fy) noexcept
    {
        double dx, dy, d;
        dx = interact->x - x;
        dy = interact->y - y;
        if (dx > table.getMaxRange() || dy > table.getMaxRange())
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
        double Bes = interact->U0 * table(d);
        *fx = Bes * dx * d1;
        *fy = Bes * dy * d1;
    }
};

#endif