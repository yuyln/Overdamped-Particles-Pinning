#ifndef __BESS0TABLE
#define __BESS0TABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#include <cmath>
#include <besselfunc.h>
#endif

#include <table.h>

Table Bess0T;

void InitBess0Table(double factor, int n, double limit)
{
    Bess0T.TN = n;
    Bess0T.TLimit = limit;
    Bess0T.TAtLimit = BESSK0(limit);
    double x = 0.0;
    double y = BESSK0(x);
    double h = 0.0001;
    while (y >= factor)
    {
        x += h;
        y = BESSK0(x);
    }
    Bess0T.TCut = x;
    Bess0T.TH = x / (double)Bess0T.TN;
    Bess0T.table = new double[Bess0T.TN + 1];
    x = 0.0;
    for (size_t i = 0; i < Bess0T.TN + 1; i++)
    {
        Bess0T.table[i] = BESSK0(x);
        x += Bess0T.TH;
    }
}

double GetBess0Table(double x)
{
    if (x >= Bess0T.TCut)
    {
        return 0.0;
    }

    if (x <= Bess0T.TLimit)
    {
        return Bess0T.TAtLimit;
    }

    size_t i = x / Bess0T.TH;
    return Bess0T.table[i];
}


#endif