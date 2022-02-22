#ifndef __TABLE
#define __TABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#endif

typedef struct
{
    double *table, TCut, TH, TLimit, TAtLimit;
    size_t TN;
} Table;

#endif