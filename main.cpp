#define PARSER_IMPLEMENTATION
#define OPENCLWRAPPER_IMPLEMTATION
#define PROFILER_IMPLEMENTATION
#include <cstdlib>
#include <cstdio>
#include <profiler.h>
#include <OpenCLWrapper.h>
#include <parser.h>
#include <functions.h>
#include <besselfunc.h>
#include <particle.h>
#include <exptable.h>
#include <bessel0table.h>
#include <bessel1table.h>
#include <pinning.h>

int main()
{   
    Pinning *p;
    Particle *part;
    InitBess0Table(1.0e-4, 1e7, 1.0);
    InitBess1Table(1.0e-4, 1e7, 1.0);
    InitExpTable(1.0e-4, 1e7, 0.0);
    int n = InitPinnings(&p);
    printf("%d\n", n);
    int n_ = InitParticles(&part);
    printf("%d\n", n_);

    for (int i = 0; i < n; i++)
    {
        printf("%.3f %.3f\n", p[i].x, p[i].y);
    }
    printf("------------------\n");
    for (int i = 0; i < n_; i++)
    {
        printf("%.3f %.3f\n", part[i].x, part[i].y);
    }
    return 0;
}