#ifndef __TABLE
#define __TABLE
#include <functional>
#include <cmath>

double FindRange(const double h, const double cutValue, const double start, const double end, std::function<double(double)> f)
{
    double x = start;
    double y = f(x);
    while (y >= cutValue && x < end)
    {
        x += h;
        y = f(x);
    }
    if (x >= end)
    {
        printf("CAUTION: X CLOSE TO END\n");
    }
    return x;
}

class Table
{
    double startValue, endValue;
    double *values, stepsize;
    double valueBelowStart, valueAfterEnd;
    size_t n;

public:
    std::function<double(double)> f;
    Table(){}

    Table(size_t n_, double start, double end,
    double valueBStart, double valueAEnd, std::function<double(double)> f_): 
    startValue(start), endValue(end), valueBelowStart(valueBStart), valueAfterEnd(valueAEnd), n(n_), f(f_)
    {
        values = new double[n + 1];
        stepsize = (end - start) / (double)n;
        double x;
        for (size_t i = 0; i <= n; ++i)
        {
            x = start + (double)i * stepsize;
            values[i] = f(x);
        }
    }

    const double &operator() (double x) const
    {
        if (x >= endValue)
        {
            return valueAfterEnd;
        }
        else if (x <= startValue)
        {
            return valueBelowStart;
        }
        else
        {
            size_t i = (size_t)((x - startValue) / stepsize);
            return values[i];
        }
    }

    const double &getMaxRange() const { return endValue; }
    const double &getStepSize() const { return stepsize; }

    void operator=(Table l2)
    {
        startValue = l2.startValue;
        endValue = l2.endValue;
        values = l2.values;
        stepsize = l2.stepsize;
        valueBelowStart = l2.valueBelowStart;
        valueAfterEnd = l2.valueAfterEnd;
        n = l2.n;
        f = l2.f;
    }
};

#endif