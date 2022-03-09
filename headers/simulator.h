#ifndef __SIM
#define __SIM

#include <matrix.h>
#include <box.h>
#include <table.h>
#include <particle.h>
#include <pinning.h>
#include <functions.h>

typedef struct Simulator
{
    Table ExpTable, BK0Table, BK1Table;
    Matrix<Box> PinningBoxes, ParticleBoxes;
    Pinning *pins;
    Particle *parts;
    Particle *parts1;

    double VXm, VYm;

    double FC;

    double omegaX, omegaY;
    double ACXFactor, ACYFactor;
    double A, B;
    double A_F, B_F;
    
    double DCFactor;
    double DCX, DCY;
    double DCAng;

    double DCFixed;
    double DCFixedX, DCFixedY;
    double DCFixedAng;

    double sqrtTemp, Lx, Ly;
    size_t nParticles, nPinnings;
    double h, tmax;
    size_t N, NCut;
    double *WriteX, *WriteY;
    bool Write;

    Simulator()
    {
        nPinnings = InitPinnings(&pins);
        nParticles = InitParticles(&parts);
        InitParticles(&parts1);
        FILE *f = fopen64("./input/input.in", "rb");
        char *data = ReadFile(f);
        fclose(f);
        int nData;
        char **d = Parse(data, &nData);
        N = (size_t)GetValueInt("N", d, nData);
        NCut = (size_t)GetValueInt("NCut", d, nData);
        Write = (size_t)GetValueInt("WRITE", d, nData);
        WriteX = new double[N * nParticles * Write / NCut];
        WriteY = new double[N * nParticles * Write / NCut];

        tmax = GetValueDouble("tmax", d, nData);
        h = tmax / (double)N;
        Lx = GetValueDouble("LX", d, nData);
        Ly = GetValueDouble("LY", d, nData);
        ACXFactor = GetValueDouble("ACXFA", d, nData);
        ACYFactor = GetValueDouble("ACYFA", d, nData);

        A_F = GetValueDouble("ACXFI", d, nData);
        B_F = GetValueDouble("ACYFI", d, nData);

        DCFactor = GetValueDouble("DCFA", d, nData);
        DCAng = GetValueDouble("DCANG", d, nData);

        DCFixed = GetValueDouble("DCFI", d, nData);
        DCFixedAng = GetValueDouble("DCANGFI", d, nData);
        DCFixedX = DCFixed * cos(DCFixedAng);
        DCFixedY = DCFixed * sin(DCFixedAng);
        sqrtTemp = 0.0;

        omegaX = tmax / GetValueDouble("NACX", d, nData);
        omegaY = tmax / GetValueDouble("NACY", d, nData);

        ExpTable = Table(100000, 1.0e-4, 0.0, [](double x)
                           { return exp(-x); });
        BK0Table = Table(100000, 1.0e-4, 1.0, [](double x)
                           { return BESSK0(x); });
        BK1Table = Table(100000, 1.0e-4, 1.0, [](double x)
                           { return BESSK1(x); });

        double R0Max = -1.0;
        for (size_t i = 0; i < nPinnings; ++i)
        {
            if (pins[i].R0 >= R0Max)
            {
                R0Max = pins[i].R0;
            }
        }

        PinningBoxes = CreateBoxes(R0Max * sqrt(ExpTable.getMaxRange()), nPinnings, Lx, Ly, pins);
        ParticleBoxes = CreateBoxes(BK1Table.getMaxRange(), nParticles, Lx, Ly, parts);
        AttBoxes(nPinnings, pins, &PinningBoxes);
    }

    void FixCurrent(double FC_)
    {
        FC = FC_;
        A = FC_ * ACXFactor;
        B = FC_ * ACYFactor;
        DCX = FC_ * DCFactor * cos(DCAng);
        DCY = FC_ * DCFactor * sin(DCAng);
    }
} Simulator;
#endif