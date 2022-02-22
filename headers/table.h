#ifndef __TABLE
#define __TABLE

#ifndef OPENCLCOMP
#include <cstdlib>
#endif

typedef struct
{
    double TCut, TH, TLimit, TAtLimit;
    size_t TN;
} TableParams;

#endif