#define PARSER_IMPLEMENTATION
#define OPENCLWRAPPER_IMPLEMTATION
#define PROFILER_IMPLEMENTATION
#include <cstdlib>
#include <cstdio>
#include <profiler.h>
#include <OpenCLWrapper.h>
#include <parser.h>
#include <dirent.h>
#include <functions.h>

int main()
{   
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/pinnings", &nfiles, &qnt, &parsed);
    double R01 = GetValueDouble("R0", parsed[0], qnt[0]);
    double R02 = GetValueDouble("R0", parsed[1], qnt[1]);
    printf("%f\n", R01);
    printf("%f\n", R02);

    int ip = GetIndexOfTag("Positions", parsed[0], qnt[0]);
    for (int i = ip + 1; i < qnt[0]; i += 2)
    {
        printf("%f %f\n", strtod(parsed[0][i], NULL), strtod(parsed[0][i + 1], NULL));
    }
    return 0;
}