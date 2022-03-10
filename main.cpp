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
    Simulator s(true);
    //TODO: save/load system           (V)
    //      output simulator object    (X)
    //      GSA                        (X) 

    if (s.Recovery)
    {
        s.LoadSystem("", "");
    }
    // s.Export("");

    StartMeasure("ALL");
    for (double FC = s.FC; FC <= s.FCMax; FC += s.hFC)
    {
        s.FixCurrent(FC);
        Integration(s, "", "");
        printf("%.4f %.4f\n", FC, s.FCMax);
    }
    EndMeasure("ALL");
    PrintAll(stdout);

    return 0;
}