#ifndef __SIM
#define __SIM

#include <matrix.h>
#include <box.h>
#include <table.h>
#include <particle.h>
#include <pinning.h>
#include <functions.h>
#include <map>

typedef struct Simulator
{
    Table PinPotentialTable, PinForceTable,
          PartPotentialTable, PartForceTable;

    Matrix<Box> PinPotentialBoxes, PinForceBoxes,
                PartPotentialBoxes, PartForceBoxes;
    Pinning *pins;
    Particle *parts;
    Particle *parts1;
    std::map<double, size_t> betaQnt;
    std::map<double, double> VmxBeta, VmyBeta;

    double VXm, VYm;

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
    size_t N, NCut, NThreads;
    double *WriteX, *WriteY;
    bool Write;
    bool Recovery;

    Simulator(bool CreateFoldersEtc)
    {
        nPinnings = InitPinnings(&pins);
        nParticles = InitParticles(&parts);
        InitParticles(&parts1);
        for (size_t i = 0; i < nParticles; ++i)
        {
            betaQnt[parts[i].betadamp] = 0;
        }

        for (size_t i = 0; i < nParticles; ++i)
        {
            betaQnt[parts[i].betadamp]++;
        }

        FILE *f = fopen("./input/input.in", "rb");
        char *data = ReadFile(f);
        fclose(f);
        int nData;
        char **d = Parse(data, &nData);
        N = (size_t)GetValueInt("N", d, nData);
        NThreads = (size_t)GetValueInt("NTHREADS", d, nData);
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

        omegaX = GetValueDouble("NACX", d, nData) / tmax;
        omegaY = GetValueDouble("NACY", d, nData) / tmax;

        double range = FindRange(0.0001, 0.21e-3, 0.0, sqrt(Lx * Lx + Ly * Ly), [](double x) { return exp(-x); });
        PinPotentialTable = Table(1e6, 0.0, range, 1.0, 0.0, [](double x){ return exp(-x); });

        range = FindRange(0.0001, 0.21e-3, 0.0, sqrt(Lx * Lx + Ly * Ly), [](double x) { return exp(-x); });
        PinForceTable = Table(1e6, 0.0, range, 1.0, 0.0, [](double x){ return exp(-x); });

        range = FindRange(0.0001, 0.21e-3, 0.0, sqrt(Lx * Lx + Ly * Ly), [](double x) { return BESSK0(x); });
        PartPotentialTable = Table(1e6, 0.0001, range, BESSK0(0.0001), 0.0, [](double x){ return BESSK0(x); });

        range = FindRange(0.0001, 0.21e-3, 0.0, sqrt(Lx * Lx + Ly * Ly), [](double x) { return BESSK1(x) / x; });
        PartForceTable = Table(1e6, 0.0001, range, BESSK1(0.0001) / 0.0001, 0.0, [](double x){ return BESSK1(x) / x; });

        double R0Max = 0.0;
        for (size_t i = 0; i < nPinnings; ++i)
        {
            if (pins[i].R0 >= R0Max)
            {
                R0Max = pins[i].R0;
            }
        }
        if (nPinnings == 0)
        {
            R0Max = 1.0;
        }
        FC = 0.0;

        PinPotentialBoxes = CreateBoxes(R0Max * sqrt(PinPotentialTable.getMaxRange()), nPinnings, Lx, Ly, pins);
        PinForceBoxes = CreateBoxes(R0Max * sqrt(PinForceTable.getMaxRange()), nPinnings, Lx, Ly, pins);
        PartPotentialBoxes = CreateBoxes(PartPotentialTable.getMaxRange(), nParticles, Lx, Ly, parts);
        PartForceBoxes = CreateBoxes(PartForceTable.getMaxRange(), nParticles, Lx, Ly, parts);

        AttBoxes(nPinnings, pins, &PinPotentialBoxes);
        AttBoxes(nPinnings, pins, &PinForceBoxes);

        if (CreateFoldersEtc)
        {
            int i = system("mkdir out");
            i = system("mkdir \"out/positions\"");
            i = system("mkdir saves");
            i = system("mkdir \"saves/snaps\"");
            i = system("mkdir \"./out/velocity_per_beta\"");
            (void)i;
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

                for (std::map<double, size_t>::iterator i = betaQnt.begin(); i != betaQnt.end(); ++i)
                {
                    double betadamp = i->first;
                    char *name;
                    size_t size = snprintf(NULL, 0, "./out/velocity_per_beta/velocity_%.5f.out", betadamp);
                    name = new char[size + 1];
                    snprintf(name, size, "./out/velocity_per_beta/velocity_%.5f.out", betadamp);
                    f = fopen(name, "wb");
                    delete[] name;
                    if (f == NULL)
                    {
                        fprintf(stderr, "NOT POSSIBLE TO OPEN FILE %s: %s\n", name, strerror(errno));
                        exit(1);
                    }
                    fclose(f);
                }
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

        size = snprintf(NULL, 0, "./saves/snaps/%s_snap_%.5f_%d_%s.snap", prefix, FC, (int)(FC / hFC), suffix);
        name = new char[size + 1];
        snprintf(name, size + 1, "./saves/snaps/%s_snap_%.5f_%d_%s.snap", prefix, FC, (int)(FC / hFC), suffix);
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
        AttBoxes(nParticles, parts, &PartPotentialBoxes);
        AttBoxes(nParticles, parts, &PartForceBoxes);
    }

    void Export(const char *name)
    {
        FILE *f = fopen(name, "wb");
        if (f == NULL)
        {
            fprintf(stderr, "NOT POSSIBLE TO OPEN FILE %s: %s", name, strerror(errno));
            exit(1);
        }
        fprintf(f, "Number of Threads: %zu\n", NThreads);
        fprintf(f, "Integration Steps: %zu\n", N);
        fprintf(f, "Integration Time: %.6f\n", tmax);
        fprintf(f, "Integration Step: %.6f\n", h);
        fprintf(f, "Current Steps: %zu\n", NCurrents);
        fprintf(f, "Max Current: %.6f\n", FCMax);
        fprintf(f, "Current Step: %.6f\n", hFC);
        fprintf(f, "Number of Particles: %zu\n", nParticles);
        fprintf(f, "Number of Pinnings: %zu\n", nPinnings);
        fprintf(f, "AC Factor in X: %.6f\n", ACXFactor);
        fprintf(f, "AC Factor in Y: %.6f\n", ACYFactor);
        fprintf(f, "Omega in X: %.6f\n", omegaX);
        fprintf(f, "Omega in Y: %.6f\n", omegaY);
        fprintf(f, "DC Factor: %.6f\n", DCFactor);
        fprintf(f, "DC Angle (degrees): %.6Lf\n", DCAng * 180.0 / M_PIl);
        fprintf(f, "AC in X Fixed: %.6f\n", A_F);
        fprintf(f, "AC in Y Fixed: %.6f\n", B_F);
        fprintf(f, "DC Fixed: %.6f\n", DCFixed);
        fprintf(f, "DC Fixed Angle: %.6Lf\n", DCFixedAng * 180.0 / M_PIl);
        fprintf(f, "Recovery Mode: %d\n", Recovery);
        fprintf(f, "Write Positions: %d\n", Write);
        fprintf(f, "Cut Factor for Writing: %zu\n", NCut);
        fprintf(f, "---------------------------\n");
        fprintf(f, "Pinning Boxes:\n");
        fprintf(f, "Potential Cutoff: %.6f\n", PinPotentialTable.getMaxRange());
        fprintf(f, "Force Cutoff: %.6f\n", PinForceTable.getMaxRange());
        fprintf(f, "Size of Boxes in X for Pinning Potential: %.6f\n", PinPotentialBoxes(0, 0).GetLx());
        fprintf(f, "Size of Boxes in Y for Pinning Potential: %.6f\n", PinPotentialBoxes(0, 0).GetLy());
        fprintf(f, "Size of Boxes in X for Pinning Force: %.6f\n", PinForceBoxes(0, 0).GetLx());
        fprintf(f, "Size of Boxes in Y for Pinning Force: %.6f\n", PinForceBoxes(0, 0).GetLy());
        fprintf(f, "Number of Boxes in X for Pinning Potential: %zu\n", PinPotentialBoxes.nCols);
        fprintf(f, "Number of Boxes in Y for Pinning Potential: %zu\n", PinPotentialBoxes.nRows);
        fprintf(f, "Number of Boxes in X for Pinning Force: %zu\n", PinForceBoxes.nCols);
        fprintf(f, "Number of Boxes in Y for Pinning Force: %zu\n", PinForceBoxes.nRows);
        fprintf(f, "---------------------------\n");
        fprintf(f, "Particle Boxes:\n");
        fprintf(f, "Potential Cutoff: %.6f\n", PartPotentialTable.getMaxRange());
        fprintf(f, "Force Cutoff: %.6f\n", PartForceTable.getMaxRange());
        fprintf(f, "Size of Boxes in X for Particle Potential: %.6f\n", PartPotentialBoxes(0, 0).GetLx());
        fprintf(f, "Size of Boxes in Y for Particle Potential: %.6f\n", PartPotentialBoxes(0, 0).GetLy());
        fprintf(f, "Size of Boxes in X for Particle Force: %.6f\n", PartForceBoxes(0, 0).GetLx());
        fprintf(f, "Size of Boxes in Y for Particle Force: %.6f\n", PartForceBoxes(0, 0).GetLy());
        fprintf(f, "Number of Boxes in X for Particle Potential: %zu\n", PartPotentialBoxes.nCols);
        fprintf(f, "Number of Boxes in Y for Particle Potential: %zu\n", PartPotentialBoxes.nRows);
        fprintf(f, "Number of Boxes in X for Particle Force: %zu\n", PartForceBoxes.nCols);
        fprintf(f, "Number of Boxes in Y for Particle Force: %zu\n", PartForceBoxes.nRows);
#if defined(RK4)
        fprintf(f, "Integration Method: Runge-Kutta 4\n");
#elif defined(RK2)
        fprintf(f, "Integration Method: Runge-Kutta 2\n");
#elif defined(EULER)
        fprintf(f, "Integration Method: Euler\n");
#endif
        fclose(f);
    }

    void ZeroVelocity()
    {
        VXm = 0.0;
        VYm = 0.0;
        for (size_t i = 0; i < nParticles; ++i)
        {
            parts[i].Vxm = 0.0;
            parts[i].Vym = 0.0;
        }

        for (std::map<double, size_t>::iterator i = betaQnt.begin(); i != betaQnt.end(); ++i)
        {
            double beta = i->first;
            VmxBeta[beta] = 0.0;
            VmyBeta[beta] = 0.0;
        }
    }

    void CalculateAverageVelocity()
    {
        for (size_t i = 0; i < nParticles; ++i)
        {
            VXm += parts[i].Vxm / (double)nParticles;
            VYm += parts[i].Vym / (double)nParticles;
        }
    }

    // TODO: change how this is implemented
    void CalculateAverageVelocityPerBeta()
    {
        for (std::map<double, size_t>::iterator i = betaQnt.begin(); i != betaQnt.end(); ++i)
        {
            double beta = i->first;
            size_t qnt = i->second;
            for (size_t j = 0; j < nParticles; ++j)
            {
                if (beta == parts[j].betadamp)
                {
                    VmxBeta[beta] += parts[j].Vxm / (double)qnt;
                    VmyBeta[beta] += parts[j].Vym / (double)qnt;
                }
            }
        }
    }
} Simulator;
#endif