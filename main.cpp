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
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/pinnings", &nfiles, &qnt, &parsed);
    double R01 = GetValueDouble("R0", parsed[0], qnt[0]);
    double R02 = GetValueDouble("R0", parsed[1], qnt[1]);
    printf("%f\n", R01);
    printf("%f\n", R02);

    int ip = GetIndexOfTag("Positions", parsed[0], qnt[0]);
    for (int i = ip + 1; i < qnt[0]; i += 2)
    {
        printf("%f %f\n", strtod(parsed[0][i], NULL), strtod(parsed[0][i + 1], NULL));
    }

    InitBess0Table(1.0e-4, 1e7, 1.0);
    InitBess1Table(1.0e-4, 1e7, 1.0);

    FILE *f = fopen("bess0.out", "wb");

    double h = 0.01;

    for (double x = -1.0; x <= 10.0; x += h)
    {
        fprintf(f, "%.3f\t%.3f\n", x, GetBess0Table(x));
    }

    fclose(f);

    f = fopen("bess1.out", "wb");

    for (double x = -1.0; x <= 10.0; x += h)
    {
        fprintf(f, "%.3f\t%.3f\n", x, GetBess1Table(x));
    }

    fclose(f);

    InitExpTable(1.0e-4, 1e7, 0.0);

    f = fopen("exp.out", "wb");

    for (double x = -1.0; x <= 10.0; x += h)
    {
        fprintf(f, "%.3f\t%.3f\n", x, GetExpTable(x * x / (0.3 * 0.3)));
    }

    fclose(f);

    f = fopen("pin.out", "wb");
    Pinning p = InitPinning(0.0, 0.0, 1.0, 0.3);
    Pinning p2 = InitPinning(1.0, 0.0, -1.0, 0.3);

    for (double x = -5.0; x <= 5.0; x += 0.05)
    {
        for (double y = -5.0; y <= 5.0; y += 0.05)
        {
            fprintf(f, "%.3f\t%.3f\t%.3f\n", x, y, PotentialPinning(&p, x, y) + PotentialPinning(&p2, x, y));
        }
    }
    fclose(f);

    return 0;
}