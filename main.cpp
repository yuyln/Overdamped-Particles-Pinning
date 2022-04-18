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
#include <omp.h>

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

    GSAParams *gsap;
    int nparams = ReadGSAParams(&gsap);
    double *values = new double[nparams];
    Particle **p = new Particle*[nparams];
    Simulator *sims = new Simulator[nparams];
    printf("Before GSA: %.15f\n", Potential(s, s.parts));

    for (size_t I = 0; I < s.rGSAP; ++I)
    {
        #pragma omp parallel for num_threads(s.nGSAT)
        for(int i = 0; i < nparams; ++i)
        {
            p[i] = GSA(gsap[i], sims[i], &values[i]);
        }

        double min = values[0];
        for (int i = 0; i < nparams; ++i)
        {
            if (values[i] <= min)
            {
                min = values[i];

                memcpy(sims[i].parts, p[i], sims[i].nParticles * sizeof(Particle));
                memcpy(sims[i].parts1, p[i], sims[i].nParticles * sizeof(Particle));
                AttBoxes(sims[i].nParticles, sims[i].parts, &sims[i].PartForceBoxes);
                AttBoxes(sims[i].nParticles, sims[i].parts, &sims[i].PartPotentialBoxes);

                memcpy(s.parts, p[i], s.nParticles * sizeof(Particle));
                memcpy(s.parts1, p[i], s.nParticles * sizeof(Particle));
                AttBoxes(s.nParticles, s.parts, &s.PartForceBoxes);
                AttBoxes(s.nParticles, s.parts, &s.PartPotentialBoxes);

            }
        }
    }
    delete[] sims;
    delete[] p;
    delete[] values;
    printf("After GSA: %.15f\n", Potential(s, s.parts));

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