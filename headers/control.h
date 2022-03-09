#ifndef __CONTROL
#define __CONTROL

#include <simulator.h>
#include <functions.h>
#include <string.h>
static const double pi = acos(-1.0);

template <typename T>
void InteractWithBox(const double &x, const double &y, const size_t iX, const size_t iY, const double &Lx, const double &Ly, 
                     const Matrix<Box> &Box, const T *p, const Table &table, double &fxO, double &fyO)
{
    double fx_ = 0.0, fy_ = 0.0;
    for (int i = -1; i <= 1; ++i)
    {
        int iBx = iX + i;
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
            int iBy = iY + i;
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

            for (size_t iB = 0; iB < Box(iBx, iBy).GetIn(); ++iB)
            {
                double fx, fy;
                int I = Box(iBx, iBy).GetIndex(iB);
                T::Force(&p[I], x - factorX * Lx, y - factorY * Ly, table, &fx, &fy);
                fx_ += fx;
                fy_ += fy;
            }
        }
    }
    fxO = fx_;
    fyO = fy_;
}

void Force(const double &x, const double &y, const size_t &i, const double &t, const Simulator &sim, double &fxO, double &fyO)
{
    double fxPart = 0.0, fyPart = 0.0;
    double fxPin = 0.0, fyPin = 0.0;
    double fx = 0.0, fy = 0.0;
    int BoxXPin = x / sim.PinningBoxes(0, 0).GetLx();
    int BoxYPin = y / sim.PinningBoxes(0, 0).GetLy();

    if (BoxXPin >= (int)sim.PinningBoxes.nCols)
    {
        BoxXPin = 0;
    }
    else if (BoxXPin < 0)
    {
        BoxXPin = sim.PinningBoxes.nCols - 1;
    }

    if (BoxYPin >= (int)sim.PinningBoxes.nRows)
    {
        BoxYPin = 0;
    }
    else if (BoxYPin < 0)
    {
        BoxYPin = sim.PinningBoxes.nRows - 1;
    }

    InteractWithBox(x, y, BoxXPin, BoxYPin, sim.Lx, sim.Ly, sim.PinningBoxes, sim.pins, sim.ExpTable, fxPin, fyPin);

    int BoxXPart = x / sim.ParticleBoxes(0, 0).GetLx();
    int BoxYPart = y / sim.ParticleBoxes(0, 0).GetLy();

    if (BoxXPart >= (int)sim.ParticleBoxes.nCols)
    {
        BoxXPart = 0;
    }
    else if (BoxXPart < 0)
    {
        BoxXPart = sim.ParticleBoxes.nCols - 1;
    }

    if (BoxYPart >= (int)sim.ParticleBoxes.nRows)
    {
        BoxYPart = 0;
    }
    else if (BoxYPart < 0)
    {
        BoxYPart = sim.ParticleBoxes.nRows - 1;
    }

    InteractWithBox(x, y, BoxXPart, BoxYPart, sim.Lx, sim.Ly, sim.ParticleBoxes, sim.parts, sim.BK0Table, fxPart, fyPart);

    fx = fxPart + fxPin;
    fy = fyPart + fyPin;

    fx += sim.DCX + sim.DCX_F;
    fy += sim.DCY + sim.DCY_F;

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
            s.parts1[i].x = 0.0;
        }
        else if (s.parts1[i].x < 0.0)
        {
            s.parts1[i].x = s.Lx;
        }

        if (s.parts1[i].y > s.Ly)
        {
            s.parts1[i].y = 0.0;
        }
        else if (s.parts1[i].y < 0.0)
        {
            s.parts1[i].y = s.Ly;
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
    double r1x, r1y
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
}

void Reset(Simulator &s)
{
    memcpy(s.parts, s.parts1, sizeof(Particle) * s.nParticles);
}
#endif