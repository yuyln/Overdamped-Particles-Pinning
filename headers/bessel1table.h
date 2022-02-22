#ifndef __BESS1TABLE
#define __BESS1TABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#include <cmath>
#include <besselfunc.h>
#endif

#include <table.h>

typedef struct Bess1Table
{
    TableParams params;
    double *table;
} Bess1Table;

Bess1Table Bess1T;

void InitBess1Table(double factor, int n, double limit)
{
    Bess1T.params.TN = n;
    Bess1T.params.TLimit = limit;
    Bess1T.params.TAtLimit = BESSK0(limit);
    double x = 0.0;
    double y = BESSK0(x);
    double h = 0.0001;
    while (y >= factor)
    {
        x += h;
        y = BESSK0(x);
    }
    Bess1T.params.TCut = x;
    Bess1T.params.TH = x / (double)Bess1T.params.TN;
    Bess1T.table = new double[Bess1T.params.TN + 1];
    x = 0.0;
    for (size_t i = 0; i < Bess1T.params.TN + 1; i++)
    {
        Bess1T.table[i] = BESSK0(x);
        x += Bess1T.params.TH;
    }
}

double GetBess1Table(double x)
{
    if (x >= Bess1T.params.TCut)
    {
        return 0.0;
    }

    if (x <= Bess1T.params.TLimit)
    {
        return Bess1T.params.TAtLimit;
    }

    size_t i = x / Bess1T.params.TH;
    return Bess1T.table[i];
}


#endif