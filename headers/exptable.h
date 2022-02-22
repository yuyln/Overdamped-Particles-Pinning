#ifndef __EXPTABLE
#define __EXPTABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#include <cmath>
#include <besselfunc.h>
#endif

#include <table.h>

typedef struct ExpTable
{
    TableParams params;
    double *table;
} ExpTable;

ExpTable ExpT;

void InitExpTable(double factor, int n, double limit)
{
    ExpT.params.TN = n;
    ExpT.params.TLimit = limit;
    ExpT.params.TAtLimit = exp(-limit);
    double x = 0.0;
    double y = exp(-x);
    double h = 0.0001;
    while (y >= factor)
    {
        x += h;
        y = exp(-x);
    }
    ExpT.params.TCut = x;
    ExpT.params.TH = x / (double)ExpT.params.TN;
    ExpT.table = new double[ExpT.params.TN + 1];
    x = 0.0;
    for (size_t i = 0; i < ExpT.params.TN + 1; i++)
    {
        ExpT.table[i] = exp(-x);
        x += ExpT.params.TH;
    }
}

double GetExpTable(double x)
{
    if (x >= ExpT.params.TCut)
    {
        return 0.0;
    }

    if (x <= ExpT.params.TLimit)
    {
        return ExpT.params.TAtLimit;
    }

    size_t i = x / ExpT.params.TH;
    return ExpT.table[i];
}


#endif