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
    char **files = FindFilesInsideDir("./input/particles", &nfiles);
    for (int i = 0; i < nfiles; i++)
    {
        printf("%s\n", files[i]);
    }
    return 0;
}