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

int ReadGSAParams(GSAParams **gsap)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/GSA", &nfiles, &qnt, &parsed);
    *gsap = new GSAParams[nfiles];

    for (int i = 0; i < nfiles; i++)
    {
        double qA = GetValueDouble("qA", parsed[i], qnt[i]);
        double qV = GetValueDouble("qV", parsed[i], qnt[i]);
        double qT = GetValueDouble("qT", parsed[i], qnt[i]);
        double T0 = GetValueDouble("T0", parsed[i], qnt[i]);
        size_t outer = (size_t)GetValueUInt("OL", parsed[i], qnt[i]);
        size_t inner = (size_t)GetValueUInt("IL", parsed[i], qnt[i]);
        size_t print = (size_t)GetValueUInt("PP", parsed[i], qnt[i]);
        (*gsap)[i].qA = qA;
        (*gsap)[i].qV = qV;
        (*gsap)[i].qT = qT;
        (*gsap)[i].T0 = T0;
        (*gsap)[i].outerLoop = outer;
        (*gsap)[i].innerLoop = inner;
        (*gsap)[i].printParam = inner / print;
    }
    return nfiles;
}

Particle* GSA(const GSAParams &param, Simulator &s, double *outEnergy)
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
        for (size_t ic = 1; ic <= param.innerLoop; ++ic)
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
            Boundary(s.parts1, s.nParticles, s.Lx, s.Ly);
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
    // delete[] pmin;
    *outEnergy = funcmin;
    return pmin;
}

#endif