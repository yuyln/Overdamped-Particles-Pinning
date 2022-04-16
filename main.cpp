#define PARSER_IMPLEMENTATION
#define PROFILER_IMPLEMENTATION
#define M_PIl          3.141592653589793238462643383279502884L
#include <b.h>
#include <cstdlib>
#include <cstdio>
#include <profiler.h>
#include <parser.h>
#include <table.h>
#include <besselfunc.h>
#include <particle.h>
#include <functions.h>
#include <box.h>
#include <matrix.h>
#include <control.h>
#include <simulator.h>
#include <GSA.h>
#include <line.h>
#include <region.h>

int main()
{
    Simulator s(true);
    s.Export("./out/simulator_data.out");
    //TODO: save/load system           (V)
    //      output simulator object    (V)
    //      GSA                        (V?) 
    //      Fix Bug no pinning         (V)
    //      Multithread GSA            (X)
    //      Multithread Integration    (V)
    //      Test generalized forces/potential (V)  <- <- Compare to paper
    //      Check every snprintf, i think i forgot a +1 on every call         (X)
    //      More general tables        (V)

    // GSAParams gsap;
    // gsap.innerLoop = 100000;
    // gsap.outerLoop = 5;
    // gsap.printParam = gsap.innerLoop / 5;
    // gsap.qA = 2.8;
    // gsap.qT = 2.2;
    // gsap.qV = 2.6;

    // gsap.T0 = 10.0;
    // GSA(gsap, s);

    // gsap.T0 = 5.0;
    // GSA(gsap, s);

    // gsap.T0 = 1.0;
    // GSA(gsap, s);

    if (s.Recovery)
    {
        s.LoadSystem("", "");
    }

    StartMeasure("ALL");
    for (double FC = s.FC; FC <= s.FCMax + s.hFC; FC += s.hFC)
    {
        s.FixCurrent(FC);
        if (s.NThreads > 1)
        {
            IntegrationMult(s, "", "");
        }
        else
        {
            Integration(s, "", "");
        }
        printf("%.4f %.4f\n", FC, s.FCMax);
    }
    EndMeasure("ALL");
    PrintAll(stdout);

    return 0;
}
