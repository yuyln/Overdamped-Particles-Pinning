#ifndef __SIM
#define __SIM

#include <matrix.h>
#include <box.h>
#include <table.h>
#include <particle.h>
#include <pinning.h>

typedef struct Simulator
{
    Table ExpTable, BK0Table, BK1Table;
    Matrix<Box> PinningBoxes, ParticleBoxes;
    Pinning *pins;
    Particle *parts;
    Particle *parts1;
    double DCX, DCY, DCX_F, DCY_F;
    double A, B, A_F, B_F;
    double sqrtTemp, Lx, Ly;
    size_t nParticles, nPinnings;
    double omegaX, omegaY;
    double h;
} Simulator;
#endif