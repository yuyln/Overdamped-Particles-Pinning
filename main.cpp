#define PARSER_IMPLEMENTATION
#define PROFILER_IMPLEMENTATION
#include <cstdlib>
#include <cstdio>
#include <profiler.h>
#include <parser.h>
#include <table.h>
#include <besselfunc.h>
#include <particle.h>
#include <pinning.h>
#include <functions.h>
#include <box.h>
#include <matrix.h>
#include <control.h>

int main()
{   
    double Lx = 36.0, Ly = Lx;

    int nParticles, nPinning;
    Particle *particles;
    Pinning *pinnings;
    nParticles = InitParticles(&particles);
    nPinning = InitPinnings(&pinnings);

    Table EXPTABLE(100000, 1.0e-4, 0.0, [](double x) { return exp(-x); });
    Table BESS0TABLE(100000, 1.0e-4, 1.0, [](double x){ return BESSK0(x); });
    Table BESS1TABLE(100000, 1.0e-4, 1.0, [](double x){ return BESSK1(x); });

    Matrix<Box> EBoxes = CreateBoxes(EXPTABLE.getMaxRange(), nPinning, Lx, Ly, pinnings);

    for (size_t j = EBoxes.nRows; j --> 0;)
    {
        for (size_t i = 0; i < EBoxes.nCols; ++i)
        {
            printf("%zu ", EBoxes(j, i).GetIn());
        }
        printf("\n");
    }

    return 0;
}