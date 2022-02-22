#ifndef __EXPTABLE
#define __EXPTABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#include <cmath>
#include <besselfunc.h>
#endif

#include <table.h>

Table ExpT;

void InitExpTable(double factor, int n, double limit)
{
    ExpT.TN = n;
    ExpT.TLimit = limit;
    ExpT.TAtLimit = exp(-limit);
    double x = 0.0;
    double y = exp(-x);
    double h = 0.0001;
    while (y >= factor)
    {
        x += h;
        y = exp(-x);
    }
    ExpT.TCut = x;
    ExpT.TH = x / (double)ExpT.TN;
    ExpT.table = new double[ExpT.TN + 1];
    x = 0.0;
    for (size_t i = 0; i < ExpT.TN + 1; i++)
    {
        ExpT.table[i] = exp(-x);
        x += ExpT.TH;
    }
}

double GetExpTable(double x)
{
    if (x >= ExpT.TCut)
    {
        return 0.0;
    }

    if (x <= ExpT.TLimit)
    {
        return ExpT.TAtLimit;
    }

    size_t i = x / ExpT.TH;
    return ExpT.table[i];
}


#endif