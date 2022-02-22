#ifndef __BESS0TABLE
#define __BESS0TABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#include <cmath>
#include <besselfunc.h>
#endif

#include <table.h>

typedef struct Bess0Table
{
    TableParams params;
    double *table;
} Bess0Table;

Bess0Table Bess0T;

void InitBess0Table(double factor, int n, double limit)
{
    Bess0T.params.TN = n;
    Bess0T.params.TLimit = limit;
    Bess0T.params.TAtLimit = BESSK0(limit);
    double x = 0.0;
    double y = BESSK0(x);
    double h = 0.0001;
    while (y >= factor)
    {
        x += h;
        y = BESSK0(x);
    }
    Bess0T.params.TCut = x;
    Bess0T.params.TH = x / (double)Bess0T.params.TN;
    Bess0T.table = new double[Bess0T.params.TN + 1];
    x = 0.0;
    for (size_t i = 0; i < Bess0T.params.TN + 1; i++)
    {
        Bess0T.table[i] = BESSK0(x);
        x += Bess0T.params.TH;
    }
}

double GetBess0Table(double x)
{
    if (x >= Bess0T.params.TCut)
    {
        return 0.0;
    }

    if (x <= Bess0T.params.TLimit)
    {
        return Bess0T.params.TAtLimit;
    }

    size_t i = x / Bess0T.params.TH;
    return Bess0T.table[i];
}


#endif