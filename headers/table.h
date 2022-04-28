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
public:
    double startValue, endValue;
    double *values, stepsize;
    double valueBelowStart, valueAfterEnd;
    size_t n;

    std::function<double(double)> f;
    Table(): values(nullptr) {}

    Table(const Table &o) = delete;
/*    {
        startValue = o.startValue; endValue = o.endValue;
        stepsize = o.stepsize;
        valueBelowStart = o.valueBelowStart;
        valueAfterEnd = o.valueAfterEnd;
        n = o.n;
        if (values)
            delete[] values;
        values = new double[n + 1];
        memcpy((void*)values, (void*)o.values, sizeof(double) * (n + 1));
    }*/

    void operator=(const Table &o)
    {
        startValue = o.startValue; endValue = o.endValue;
        stepsize = o.stepsize;
        valueBelowStart = o.valueBelowStart;
        valueAfterEnd = o.valueAfterEnd;
        n = o.n;
        if (values)
            delete[] values;
        values = new double[n + 1];
        memcpy((void*)values, (void*)o.values, sizeof(double) * (n + 1));
    }

    ~Table()
    {
        if (values)
            delete[] values;
    }

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
};

#endif