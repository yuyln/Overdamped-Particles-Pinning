#ifndef __CONTROL
#define __CONTROL

#include <simulator.h>
#include <functions.h>
#include <string.h>
#include <errno.h>
#include <omp.h>
static const double pi = acos(-1.0);

int FindBox(const double &p, const double &BoxSize, const size_t &nMax)
{
    int Box = p / BoxSize;

    if (Box >= (int)nMax)
    {
        Box = 0;
    }
    else if (Box < 0)
    {
        Box = nMax - 1;
    }
    return Box;
}

template <typename T, bool ig=false>
void InteractWithBox(const size_t &ic, const double &x, const double &y, const size_t iX, const size_t iY, const double &Lx, const double &Ly, 
                     const Matrix<Box> &Box, const T *p, const Table &table, double &fxO, double &fyO)
{
    double fx_ = 0.0, fy_ = 0.0;
    for (int i = -1; i <= 1; ++i)
    {
        int iBx = (int)iX + i;
        bool left = iBx < 0;
        bool right = iBx >= (int)Box.nCols;
        if (right)
        {
            iBx = 0;
        }
        else if (left)
        {
            iBx = Box.nCols - 1;
        }
        int factorX = right - left;
        for (int j = -1; j <= 1; ++j)
        {
            int iBy = (int)iY + j;
            bool down = iBy < 0;
            bool up = iBy >= (int)Box.nRows;
            if (up)
            {
                iBy = 0;
            }
            else if (down)
            {
                iBy = Box.nRows - 1;
            }
            int factorY = up - down;

            for (size_t iB = 0; iB < Box(iBy, iBx).GetIn(); ++iB)
            {
                double fx, fy;
                size_t I = Box(iBy, iBx).GetIndex(iB);
                if (ig) { if (I == ic) { continue; } }
                T::Force(&p[I], x - factorX * Lx, y - factorY * Ly, table, &fx, &fy);
                fx_ += fx;
                fy_ += fy;
            }
        }
    }
    fxO = fx_;
    fyO = fy_;
}


template <typename T, bool ig=false>
double InteractWithBoxPotential(const size_t &ic, const double &x, const double &y, const size_t iX, const size_t iY, const double &Lx, const double &Ly, 
                     const Matrix<Box> &Box, const T *p, const Table &table)
{
    double ret = 0.0;
    for (int i = -1; i <= 1; ++i)
    {
        int iBx = (int)iX + i;
        bool left = iBx < 0;
        bool right = iBx >= (int)Box.nCols;
        if (right)
        {
            iBx = 0;
        }
        else if (left)
        {
            iBx = Box.nCols - 1;
        }
        int factorX = right - left;
        for (int j = -1; j <= 1; ++j)
        {
            int iBy = (int)iY + j;
            bool down = iBy < 0;
            bool up = iBy >= (int)Box.nRows;
            if (up)
            {
                iBy = 0;
            }
            else if (down)
            {
                iBy = Box.nRows - 1;
            }
            int factorY = up - down;

            for (size_t iB = 0; iB < Box(iBy, iBx).GetIn(); ++iB)
            {
                size_t I = Box(iBy, iBx).GetIndex(iB);
                if (ig) { if (I == ic) { continue; } }
                ret += T::Potential(&p[I], x - factorX * Lx, y - factorY * Ly, table);
            }
        }
    }
    return ret;
}

double Potential(const Simulator &s, const Particle *p)
{
    double retPin = 0.0, retPart = 0.0;
    for (size_t i = 0; i < s.nParticles; ++i)
    {
        int BoxXPin = FindBox(p[i].x, s.PinPotentialBoxes(0, 0).GetLx(), s.PinPotentialBoxes.nCols);
        int BoxYPin = FindBox(p[i].y, s.PinPotentialBoxes(0, 0).GetLy(), s.PinPotentialBoxes.nRows);

        retPin += InteractWithBoxPotential<Pinning, false>(i, p[i].x, p[i].y, BoxXPin, BoxYPin, s.Lx, s.Ly, s.PinPotentialBoxes, s.pins, s.PinPotentialTable);
    }

    for (size_t i = 0; i < s.nParticles; ++i)
    {
        int BoxXPart = FindBox(p[i].x, s.PartPotentialBoxes(0, 0).GetLx(), s.PartPotentialBoxes.nCols);
        int BoxYPart = FindBox(p[i].y, s.PartPotentialBoxes(0, 0).GetLy(), s.PartPotentialBoxes.nRows);

        retPart += InteractWithBoxPotential<Particle, true>(i, p[i].x, p[i].y, BoxXPart, BoxYPart, s.Lx, s.Ly, s.PartPotentialBoxes, p, s.PartPotentialTable);
    }

    return retPin + retPart / 2.0;
}

/*void InteractWithBoxParticle(const size_t &ic, const double &x, const double &y, const size_t iX, const size_t iY, const double &Lx, const double &Ly, 
                     const Matrix<Box> &Box, const Particle *p, const Table &table, double &fxO, double &fyO)
{
    double fx_ = 0.0, fy_ = 0.0;
    for (int i = -1; i <= 1; ++i)
    {
        int iBx = (int)iX + i;
        bool left = iBx < 0;
        bool right = iBx >= (int)Box.nCols;
        if (right)
        {
            iBx = 0;
        }
        else if (left)
        {
            iBx = Box.nCols - 1;
        }
        int factorX = right - left;
        for (int j = -1; j <= 1; ++j)
        {
            int iBy = (int)iY + j;
            bool down = iBy < 0;
            bool up = iBy >= (int)Box.nRows;
            if (up)
            {
                iBy = 0;
            }
            else if (down)
            {
                iBy = Box.nRows - 1;
            }
            int factorY = up - down;

            for (size_t iB = 0; iB < Box(iBy, iBx).GetIn(); ++iB)
            {
                double fx, fy;
                size_t I = Box(iBy, iBx).GetIndex(iB);
                // if (ignore) { if (I == ic) { continue; } }
                Particle::Force(&p[I], x - factorX * Lx, y - factorY * Ly, table, &fx, &fy);
                fx_ += fx;
                fy_ += fy;
            }
        }
    }
    fxO = fx_;
    fyO = fy_;
}*/

void Force(const double &x, const double &y, const size_t &i, const double &t, const Simulator &sim, double &fxO, double &fyO)
{
    double fxPart = 0.0, fyPart = 0.0;
    double fxPin = 0.0, fyPin = 0.0;
    double fx = 0.0, fy = 0.0;
    int BoxXPin = FindBox(x, sim.PinForceBoxes(0, 0).GetLx(), sim.PinForceBoxes.nCols);
    int BoxYPin = FindBox(y, sim.PinForceBoxes(0, 0).GetLy(), sim.PinForceBoxes.nRows);

    InteractWithBox<Pinning, false>(i, x, y, BoxXPin, BoxYPin, sim.Lx, sim.Ly, sim.PinForceBoxes, sim.pins, sim.PinForceTable, fxPin, fyPin);

    int BoxXPart = FindBox(x, sim.PartForceBoxes(0, 0).GetLx(), sim.PartForceBoxes.nCols);
    int BoxYPart = FindBox(y, sim.PartForceBoxes(0, 0).GetLy(), sim.PartForceBoxes.nRows);

    InteractWithBox<Particle, true>(i, x, y, BoxXPart, BoxYPart, sim.Lx, sim.Ly, sim.PartForceBoxes, sim.parts, sim.PartForceTable, fxPart, fyPart);

    fx = fxPart + fxPin;
    fy = fyPart + fyPin;

    fx += sim.DCX + sim.DCFixedX;
    fy += sim.DCY + sim.DCFixedY;

    fx += (sim.A + sim.A_F) * sin(2.0 * pi * sim.omegaX * t);
    fy += (sim.B + sim.B_F) * cos(2.0 * pi * sim.omegaY * t);

    fx += (2.0 * myrandom() - 1.0) * sim.sqrtTemp;
    fy += (2.0 * myrandom() - 1.0) * sim.sqrtTemp;

    double damp = sim.parts[i].damp;
    double beta = sim.parts[i].beta;

    fxO = (damp * fx + beta * fy) / (beta * beta + damp * damp);
    fyO = (damp * fy - beta * fx) / (beta * beta + damp * damp);
}

void Boundary(Simulator &s)
{
    for (size_t i = 0; i < s.nParticles; ++i)
    {
        if (s.parts1[i].x > s.Lx)
        {
            s.parts1[i].x -= s.Lx;
        }
        else if (s.parts1[i].x < 0.0)
        {
            s.parts1[i].x += s.Lx;
        }

        if (s.parts1[i].y > s.Ly)
        {
            s.parts1[i].y -= s.Ly;
        }
        else if (s.parts1[i].y < 0.0)
        {
            s.parts1[i].y += s.Ly;
        }
    }
}

void Step(const size_t &i, const double &t, Simulator &s, double &fxO, double &fyO)
{
    #if defined(RK4)
    double r1x, r1y, r2x, r2y, r3x, r3y, r4x, r4y;
    Force(s.parts[i].x, s.parts[i].y, i, t, s, r1x, r1y);
    Force(s.parts[i].x + r1x * s.h / 2.0, s.parts[i].y + r1y * s.h / 2.0, i, t + s.h / 2.0, s, r2x, r2y);
    Force(s.parts[i].x + r2x * s.h / 2.0, s.parts[i].y + r2y * s.h / 2.0, i, t + s.h / 2.0, s, r3x, r3y);
    Force(s.parts[i].x + r3x * s.h, s.parts[i].y + r3y * s.h, i, t + s.h, s, r4x, r4y);
    fxO = (r1x + 2.0 * r2x + 2.0 * r3x + r4x) * s.h / 6.0;
    fyO = (r1y + 2.0 * r2y + 2.0 * r3y + r4y) * s.h / 6.0;
    #elif defined(RK2)
    double r1x, r1y, r2x, r2y;
    Force(s.parts[i].x, s.parts[i].y, i, t, s, r1x, r1y);
    Force(s.parts[i].x + r1x * s.h, s.parts[i].y + r1y * s.h, i, t + s.h, s, r2x, r2y);
    fxO = (r1x + r2x) * s.h / 2.0;
    fyO = (r1y + r2y) * s.h / 2.0;
    #elif defined(EULER)
    double r1x, r1y;
    Force(s.parts[i].x, s.parts[i].y, i, t, s, r1x, r1y);
    fxO = (r1x) * s.h;
    fyO = (r1y) * s.h;
    #else
    fprintf(stderr, "CHOOSE ONE INTEGRATION ON COMPILATION -DRK4 -DRK2 -DEULER");
    exit(1);
    #endif
}

void Att(const size_t &i, const double &t, Simulator &s)
{
    double fx, fy;
    Step(i, t, s, fx, fy);
    s.parts1[i].x = s.parts[i].x + fx;
    s.parts1[i].y = s.parts[i].y + fy;

    s.parts[i].Vx = (s.parts1[i].x - s.parts[i].x) / s.h;
    s.parts[i].Vy = (s.parts1[i].y - s.parts[i].y) / s.h;
    s.parts[i].Vxm += s.parts[i].Vx / (double)s.N;
    s.parts[i].Vym += s.parts[i].Vy / (double)s.N;

    s.parts1[i].Vxm = s.parts[i].Vxm;
    s.parts1[i].Vym = s.parts[i].Vym;
}

void Reset(Simulator &s)
{
    memcpy(s.parts, s.parts1, sizeof(Particle) * s.nParticles);
}

void AttWrite(Simulator &s, const size_t &i)
{
    if (s.Write == 0 || (i % s.NCut) != 0)
    {
        return;
    }
    size_t i_ = i / s.NCut;
    for (size_t ip = 0; ip < s.nParticles; ++ip)
    {
        s.WriteX[i_ * s.nParticles + ip] = s.parts[ip].x;
        s.WriteY[i_ * s.nParticles + ip] = s.parts[ip].y;
    }
}

void WriteToFile(Simulator &s, const char *prefix, const char *suffix)
{
    if (s.Write == 0) { return; }
    int size = snprintf(NULL, 0, "./out/positions/%s_A_%.4f_AF_%.4f_B_%.4f_BF_%.4f_DC_%.4f_DCF_%.4f_%s.position", 
                           prefix, s.A, s.A_F, s.B, s.B_F, s.DCFactor * s.FC, s.DCFixed, suffix);
    char *name = new char[size];
    snprintf(name, size + 1, "./out/positions/%s_A_%.4f_AF_%.4f_B_%.4f_BF_%.4f_DC_%.4f_DCF_%.4f_%s.position", 
                           prefix, s.A, s.A_F, s.B, s.B_F, s.DCFactor * s.FC, s.DCFixed, suffix);
    FILE *f = fopen(name, "wb");
    if (f == NULL)
    {
        fprintf(stderr, "NOT POSSIBLE TO OPEN FILE %s: %s", name, strerror(errno));
        exit(1);
    }

    for (size_t t = 0; t < s.N; ++t)
    {
        if (t % s.NCut != 0) { continue; }
        size_t t_ = t / s.NCut;
        for (size_t i = 0; i < s.nParticles - 1; ++i)
        {
            fprintf(f, "%.5f\t%.5f\t", s.WriteX[t_ * s.nParticles + i], s.WriteY[t_ * s.nParticles + i]);
        }
        size_t i = s.nParticles - 1;
        fprintf(f, "%.5f\t%.5f\n", s.WriteX[t_ * s.nParticles + i], s.WriteY[t_ * s.nParticles + i]);
    }
    fclose(f);
}

void Integration(Simulator &s, const char *prefixSave, const char *suffixSave)
{
    s.ZeroVelocity();
    s.SaveSystem(prefixSave, suffixSave);
    for (size_t i = 0; i < s.N; ++i)
    {
        AttBoxes(s.nParticles, s.parts, &s.PartForceBoxes);
        double t = i * s.h;
        #pragma omp parallel for num_threads(s.NThreads)
        for (size_t ip = 0; ip < s.nParticles; ++ip)
        {
            Att(ip, t, s);
        }
        Boundary(s);
        Reset(s);
        AttWrite(s, i);
        if (i % s.NCut == 0) { printf("%.2f\n", (double)i / (double)s.N * 100.0); }
    }
    WriteToFile(s, prefixSave, suffixSave);
    FILE *vel = fopen("./out/velocity.out", "a");
    s.CalculateAverageVelocity();
    fprintf(vel, "%.15f\t%.15f\t%.15f\n", s.FC, s.VXm, s.VYm);
    fclose(vel);

    s.CalculateAverageVelocityPerBeta();
    for (std::map<double, size_t>::iterator i = s.betaQnt.begin(); i != s.betaQnt.end(); ++i)
    {
        double betadamp = i->first;
        char *name;
        size_t size = snprintf(NULL, 0, "./out/velocity_per_beta/velocity_%.5f.out", betadamp);
        name = new char[size + 1];
        snprintf(name, size, "./out/velocity_per_beta/velocity_%.5f.out", betadamp);
        FILE *f = fopen(name, "a");
        delete[] name;
        if (f == NULL)
        {
            fprintf(stderr, "NOT POSSIBLE TO OPEN FILE %s: %s\n", name, strerror(errno));
            exit(1);
        }
        fprintf(f, "%.15f\t%.15f\t%.15f\n", s.FC, s.VmxBeta[betadamp], s.VmyBeta[betadamp]);
        fclose(f);
    }

}
#endif