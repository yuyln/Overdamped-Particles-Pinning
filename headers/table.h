#ifndef __TABLE
#define __TABLE
#include <functional>
#include <cmath>

class Table
{
    const double zero = 0.0;
    double minrange, maxrange, *values, stepsize;
    size_t n;

public:
    std::function<double(double)> f;
    Table(){}

    Table(size_t n_, double minValue, double minrange_, std::function<double(double)> f): minrange(minrange_), n(n_)
    {
        values = new double[n + 1];

        double h = 0.0001;
        double x = 0.0;
        double y = f(x);
        while (y >= minValue)
        {
            x += h;
            y = f(x);
        }
        maxrange = x;
        stepsize = x / (double)n;
        x = 0.0;
        for (size_t i = 0; i <= n; ++i)
        {
            x = (double)i * stepsize;
            if (x < minrange)
            {
                x = minrange;
            }
            values[i] = f(x);
        }
    }

    const double &operator() (double x) const
    {
        if (x >= maxrange)
        {
            return zero;
        }
        else if (x <= minrange)
        {
            return values[0];
        }
        else
        {
            size_t i = (size_t)(x / stepsize);
            return values[i];
        }
    }

    const double &getMaxRange() const { return maxrange; }
    const double &getStepSize() const { return stepsize; }

    void operator=(Table l2)
    {
        minrange = l2.minrange;
        maxrange = l2.maxrange;
        values = l2.values;
        stepsize = l2.stepsize;
        n = l2.n;
    }
};

#endif