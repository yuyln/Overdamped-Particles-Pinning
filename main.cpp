#define PARSER_IMPLEMENTATION
#define PROFILER_IMPLEMENTATION
#define M_PIl          3.141592653589793238462643383279502884L
#include <cstdlib>
#include <cstdio>
#include <profiler.h>
#include <parser.h>
#include <table.h>
#include <besselfunc.h>
#include <particle.h>
#include <pinning.h>
#include <functions.h>
#include <box.h>
#include <matrix.h>
#include <control.h>
#include <simulator.h>

int main()
{
    FILE *f = fopen("./out/velocity.out", "w");
    fclose(f);
    Simulator s;
    //TODO: save/load system
    //      output simulator object

    StartMeasure("ALL");
    for (int t = 0; t < 1000; ++t)
    {
        double FC = t * 0.003;
        s.FixCurrent(FC);
        Integration(s);
        printf("%.2f\n", (double)t / (double)1000 * 100.0);
    }
    EndMeasure("ALL");
    PrintAll(stdout);

    return 0;
}