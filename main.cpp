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
    //      Multithread GSA            (V)
    //      Multithread Integration    (V)
    //      Test generalized forces/potential (V)  <- <- Compare to paper
    //      Check every snprintf, i think i forgot a +1 on every call         (X)
    //      More general tables        (V)
    //      Fix Regions Boxes          (X)

    GSAParams *gsap;
    int nparams = ReadGSAParams(&gsap);
    double *values = new double[nparams];
    Particle **p = new Particle*[nparams];
    Simulator *sims = new Simulator[nparams];
    printf("Before GSA: %.15f\n", Potential(s, s.parts));

    for (size_t I = 0; I < s.rGSAP * (s.Recovery == 0); ++I)
    {
        #pragma omp parallel for num_threads(s.nGSAT)
        for(int i = 0; i < nparams; ++i)
        {
            p[i] = GSA(gsap[i], sims[i], &values[i]);
        }

        size_t minI = 0;
        double min = values[0];
        for (int i = 0; i < nparams; ++i)
        {
            if (values[i] <= min)
            {
                min = values[i];
                minI = i;
            }
        }

        for (int i = 0; i < nparams; ++i)
        {
            memcpy(sims[i].parts, p[minI], sims[i].nParticles * sizeof(Particle));
            memcpy(sims[i].parts1, p[minI], sims[i].nParticles * sizeof(Particle));
            AttBoxes(sims[i].nParticles, sims[i].parts, &sims[i].PartForceBoxes);
            AttBoxes(sims[i].nParticles, sims[i].parts, &sims[i].PartPotentialBoxes);

            memcpy(s.parts, p[minI], s.nParticles * sizeof(Particle));
            memcpy(s.parts1, p[minI], s.nParticles * sizeof(Particle));
            AttBoxes(s.nParticles, s.parts, &s.PartForceBoxes);
            AttBoxes(s.nParticles, s.parts, &s.PartPotentialBoxes);
        }

        FILE *f = fopen("./out/GSAParticles.out", "wb");
        if (f == NULL)
        {
            fprintf(stderr, "NOT POSSIBLE: %s\n", strerror(errno));
            exit(1);
        }
        for (size_t i = 0; i < s.nParticles - 1; ++i)
        {
            fprintf(f, "%.15f\t%.15f\n", s.parts[i].x, s.parts[i].y);
        }
        size_t i = s.nParticles - 1;
        fprintf(f, "%.15f\t%.15f", s.parts[i].x, s.parts[i].y);
        fclose(f);


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