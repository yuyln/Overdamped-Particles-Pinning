#define PARSER_IMPLEMENTATION
#define PROFILER_IMPLEMENTATION
#include <cstdlib>
#include <cstdio>
#include <profiler.h>
#include <parser.h>
#include <table.h>
#include <besselfunc.h>

int main()
{   
    Table EXPTABLE(100000, 1.0e-4, 0.0, [](double x) { return exp(-x); });
    Table BESS0TABLE(100000, 1.0e-4, 1.0, [](double x){ return BESSK0(x); });
    Table BESS1TABLE(100000, 1.0e-4, 1.0, [](double x){ return BESSK1(x); });
    FILE *f = fopen("exp.out", "w");
    for (double x = 0; x < 5.0; x += 0.01)
    {
        fprintf(f, "%.3f\t%.3f\n", x, EXPTABLE(x * x / (0.3 * 0.3)));
    }
    fclose(f);
    printf("%f\n", EXPTABLE.getMaxRange());
    printf("%f\n", EXPTABLE.getStepSize());

    FILE *f_ = fopen("bess0.out", "w");
    for (double x = 0; x < 5.0; x += 0.01)
    {
        fprintf(f_, "%.15f\t%.15f\n", x, BESS0TABLE(x));
    }
    fclose(f_);
    printf("%f\n", BESS0TABLE.getMaxRange());
    printf("%f\n", BESS0TABLE.getStepSize());

    FILE *f__ = fopen("bess1.out", "w");
    for (double x = 0; x < 5.0; x += 0.01)
    {
        fprintf(f__, "%.15f\t%.15f\n", x, BESS1TABLE(x));
    }
    fclose(f__);
    printf("%f\n", BESS1TABLE.getMaxRange());
    printf("%f\n", BESS1TABLE.getStepSize());
    // Pinning *p;
    // Particle *part;
    // int n = InitPinnings(&p);
    // int n_ = InitParticles(&part);
    return 0;
}