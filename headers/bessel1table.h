#ifndef __BESS1TABLE
#define __BESS1TABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#include <cmath>
#include <besselfunc.h>
#endif

#include <table.h>

Table Bess1T;

void InitBess1Table(double factor, int n, double limit)
{
    Bess1T.TN = n;
    Bess1T.TLimit = limit;
    Bess1T.TAtLimit = BESSK1(limit);
    double x = 0.0;
    double y = BESSK1(x);
    double h = 0.0001;
    while (y >= factor)
    {
        x += h;
        y = BESSK1(x);
    }
    Bess1T.TCut = x;
    Bess1T.TH = x / (double)Bess1T.TN;
    Bess1T.table = new double[Bess1T.TN + 1];
    x = 0.0;
    for (size_t i = 0; i < Bess1T.TN + 1; i++)
    {
        Bess1T.table[i] = BESSK1(x);
        x += Bess1T.TH;
    }
}

double GetBess1Table(double x)
{
    if (x >= Bess1T.TCut)
    {
        return 0.0;
    }

    if (x <= Bess1T.TLimit)
    {
        return Bess1T.TAtLimit;
    }

    size_t i = x / Bess1T.TH;
    return Bess1T.table[i];
}


#endif