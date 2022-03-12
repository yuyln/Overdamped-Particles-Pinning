#ifndef __GSA
#define __GSA
#include <simulator.h>
#include <control.h>
#include <particle.h>

class GSAParams
{
public:
    double qA, qT, qV, T0;
    size_t outerLoop, innerLoop, printParam;
};

void GSA(const GSAParams &param, Simulator &s)
{
    double qA1, qT1, qV1, OneqA1, coef, coef1, r, pqa, df, tmp, exp1, exp2;
    Particle *pmin = new Particle[s.nParticles];
    double func0, func1, funcmin;
    double T, Tup, t, gammaUp, gammaDown, Tqt;
    memcpy(pmin, s.parts, sizeof(Particle) * s.nParticles);
    AttBoxes(s.nParticles, s.parts, &s.PartPotentialBoxes);
    funcmin = Potential(s, s.parts);
    func0 = funcmin;
    func1 = funcmin;
    printf("Starting Energy: %.15f\n", funcmin);

    qA1 = param.qA - 1.0;
    qT1 = param.qT - 1.0;
    qV1 = param.qV - 1.0;
    Tqt = param.T0 * (pow(2.0, qT1) - 1.0);
    tmp = 1.0 / qV1 - 0.5;
    gammaDown = tgamma(tmp);
    OneqA1 = 1.0 / qA1;
    exp1 = 2.0 / (3.0 - param.qV);
    coef1 = 1.0;
    exp2 = 1.0 / qV1 - 0.5;
    gammaUp = gammaDown;
    coef = coef1 * gammaUp / gammaDown;

    for (size_t i = 0; i < param.outerLoop; ++i)
    {
        t = 0.0;
        srand((unsigned)time(NULL));
        for (size_t ic = 0; ic <= param.innerLoop; ++ic)
        {
            t = t + 1.0;
            if (ic % param.printParam == 0) { printf("Inner: %zu  Outer: %zu   Minimun: %.15f\n", ic, i, funcmin); }
            T = Tqt / (pow(t + 1.0, qT1) - 1.0);
            Tup = 1.0;

            for (size_t ip = 0; ip < s.nParticles; ++ip)
            {
                double R = myrandom();
                double S = myrandom();
                double delta = coef * Tup / pow(1.0 + qV1 * R * R / pow(T, exp1), exp2);
                if(S <= 0.5){delta = -delta;}
                s.parts1[ip].x = s.parts[ip].x + delta;

                R = myrandom();
                S = myrandom();
                delta = coef * Tup / pow(1.0 + qV1 * R * R / pow(T, exp1), exp2);
                if(S <= 0.5){delta = -delta;}
                s.parts1[ip].y = s.parts[ip].y + delta;
            }
            Boundary(s);
            AttBoxes(s.nParticles, s.parts1, &s.PartPotentialBoxes);
            func1 = Potential(s, s.parts1);

            if (func1 < funcmin)
            {
                funcmin = func1;
                memcpy(pmin, s.parts1, sizeof(Particle) * s.nParticles);
                memcpy(s.parts, s.parts1, sizeof(Particle) * s.nParticles);
            }

            if (func1 <= func0)
            {
                memcpy(s.parts, s.parts1, sizeof(Particle) * s.nParticles);
                func0 = func1;
            }
            else
            {
                df = func1 - func0;
                pqa = 1.0 / pow(1.0 + qA1 * df / T, OneqA1);
                r = myrandom();
                if (r < pqa)
                {
                    memcpy(s.parts, s.parts1, sizeof(Particle) * s.nParticles);
                    func0 = func1;
                }
            }
        }
    }
    memcpy(s.parts, pmin, sizeof(Particle) * s.nParticles);
    memcpy(s.parts1, pmin, sizeof(Particle) * s.nParticles);
    delete[] pmin;

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

#endif