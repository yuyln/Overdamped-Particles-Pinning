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
    double *Vxm_ind, *Vym_ind;

    double FC, FCMax, hFC;
    size_t NCurrents;

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
    bool Recovery;

    Simulator(bool CreateFoldersEtc)
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
        Write = (bool)GetValueInt("WRITE", d, nData);
        Recovery = (bool)GetValueInt("RECOVERY", d, nData);
        WriteX = new double[N * nParticles * Write / NCut];
        WriteY = new double[N * nParticles * Write / NCut];

        tmax = GetValueDouble("tmax", d, nData);
        h = tmax / (double)N;
        Lx = GetValueDouble("LX", d, nData);
        Ly = GetValueDouble("LY", d, nData);
        ACXFactor = GetValueDouble("ACXFA", d, nData);
        ACYFactor = GetValueDouble("ACYFA", d, nData);

        FCMax = GetValueDouble("FCMAX", d, nData);
        NCurrents = (size_t)GetValueInt("CURRSTEPS", d, nData);
        hFC = FCMax / (double)NCurrents;

        A_F = GetValueDouble("ACXFI", d, nData);
        B_F = GetValueDouble("ACYFI", d, nData);

        DCFactor = GetValueDouble("DCFA", d, nData);
        DCAng = GetValueDouble("DCANG", d, nData) * M_PIl / 180.0;

        DCFixed = GetValueDouble("DCFI", d, nData);
        DCFixedAng = GetValueDouble("DCANGFI", d, nData) * M_PIl / 180.0;
        DCFixedX = DCFixed * cos(DCFixedAng);
        DCFixedY = DCFixed * sin(DCFixedAng);
        sqrtTemp = sqrt(GetValueDouble("TEMP", d, nData));

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
        FC = 0.0;
        PinningBoxes = CreateBoxes(R0Max * sqrt(ExpTable.getMaxRange()), nPinnings, Lx, Ly, pins);
        ParticleBoxes = CreateBoxes(BK1Table.getMaxRange(), nParticles, Lx, Ly, parts);
        AttBoxes(nPinnings, pins, &PinningBoxes);

        if (CreateFoldersEtc)
        {
            system("mkdir out");
            system("mkdir \"out/positions\"");
            system("mkdir saves");
            system("mkdir \"saves/snaps\"");
            if (!Recovery)
            {
                FILE *f = fopen("./out/velocity.out", "wb");
                if (f == NULL)
                {
                    fprintf(stderr, "NOT POSSIBLE TO OPEN VELOCITY FILE: %s\n", strerror(errno));
                    exit(1);
                }
                fclose(f);

                f = fopen("./out/velocities.out", "wb");
                if (f == NULL)
                {
                    fprintf(stderr, "NOT POSSIBLE TO OPEN VELOCITIES FILE: %s\n", strerror(errno));
                    exit(1);
                }
                fclose(f);
            }
        }
    }

    void FixCurrent(double FC_)
    {
        FC = FC_;
        A = FC_ * ACXFactor;
        B = FC_ * ACYFactor;
        DCX = FC_ * DCFactor * cos(DCAng);
        DCY = FC_ * DCFactor * sin(DCAng);
    }

    void SaveSystem(const char *prefix, const char *suffix)
    {
        size_t size = snprintf(NULL, 0, "./saves/%s_save_%s.save", prefix, suffix);
        char *name = new char[size + 1];
        snprintf(name, size + 1, "./saves/%s_save_%s.save", prefix, suffix);
        FILE *f = fopen(name, "wb");
        delete[] name;

        if (f == NULL)
        {
            fprintf(stderr, "NOT POSSIBLE TO OPEN SAVE FILE: %s\n", strerror(errno));
            exit(1);
        }
        fprintf(f, "Current to do STARTINGCURRENT: %.15f\n", FC);
        fprintf(f, "Positions to start STARTINGPOSITIONS:\n");

        for (size_t i = 0; i < nParticles - 1; ++i)
        {
            fprintf(f, "%.17f\t%.17f\n", parts[i].x, parts[i].y);
        }
        size_t i = nParticles - 1;
        fprintf(f, "%.17f\t%.17f\n", parts[i].x, parts[i].y);
        fclose(f);


        size = snprintf(NULL, 0, "./saves/snaps/%s_snap_%.5f_%s.snap", prefix, FC, suffix);
        name = new char[size + 1];
        snprintf(name, size + 1, "./saves/snaps/%s_snap_%.5f_%s.snap", prefix, FC, suffix);
        f = fopen(name, "wb");
        delete[] name;
        if (f == NULL)
        {
            fprintf(stderr, "NOT POSSIBLE TO OPEN SAVE FILE: %s\n", strerror(errno));
            exit(1);
        }
        fprintf(f, "Current to do STARTINGCURRENT: %.15f\n", FC);
        fprintf(f, "Positions to start STARTINGPOSITIONS:\n");

        for (size_t i = 0; i < nParticles - 1; ++i)
        {
            fprintf(f, "%.17f\t%.17f\n", parts[i].x, parts[i].y);
        }
        i = nParticles - 1;
        fprintf(f, "%.17f\t%.17f\n", parts[i].x, parts[i].y);
        fclose(f);
    }

    void LoadSystem(const char *prefix, const char *suffix)
    {
        size_t size = snprintf(NULL, 0, "./saves/%s_save_%s.save", prefix, suffix);
        char *name = new char[size + 1];
        snprintf(name, size + 1, "./saves/%s_save_%s.save", prefix, suffix);
        FILE *f = fopen(name, "rb");
        delete[] name;

        if (f == NULL)
        {
            fprintf(stderr, "NOT POSSIBLE TO OPEN SAVE FILE: %s\n", strerror(errno));
            exit(1);
        }

        char *data = ReadFile(f);
        fclose(f);
        int nD;
        char **d = Parse(data, &nD);
        FC = GetValueDouble("STARTINGCURRENT", d, nD);
        size_t iStart = GetIndexOfTag("STARTINGPOSITIONS", d, nD) + 1;
        size_t i_ = 0;
        for (int i = iStart; i < nD; i += 2)
        {
            parts[i_].x = strtod(d[i], NULL);
            parts[i_].y = strtod(d[i + 1], NULL);
            parts1[i_].x = strtod(d[i], NULL);
            parts1[i_].y = strtod(d[i + 1], NULL);
            ++i_;
        }
        delete[] d;
        AttBoxes(nParticles, parts, &ParticleBoxes);
    }
} Simulator;
#endif