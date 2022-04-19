#ifndef __FUNCS
#define __FUNCS
#define PARSER_IMPLEMENTATION
#define PROFILER_IMPLEMENTATION
#ifndef _MSC_VER
#include <dirent.h>
#else
#include <dirent_.h>
#endif
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <parser.h>
#include <particle.h>
#include <matrix.h>
#include <box.h>
#include <line.h>

int FindNumberFilesInsideDir(const char *dir);
char **FindFilesInsideDir(const char *dir, int *nf);
void ParseFilesInsideDir(const char *dir, int *nfiles, int **qnt, char ****parsed);
int InitParticles(Particle **p);
int InitLines(LineSegment **l);

inline double myrandom() noexcept
{
    return (double)rand() / RAND_MAX;
}

int FindNumberFilesInsideDir(const char *dir_)
{
    DIR *d;
    struct dirent *dir;
    d = opendir(dir_);
    int i = 0;
    if (d != NULL)
    {
        while ((dir = readdir(d)) != NULL)
        {
            int r1 = strcmp(".", dir->d_name);
            int r2 = strcmp("..", dir->d_name);
            if (r1 != 0 && r2 != 0)
            {
                i++;
            }
        }
        closedir(d);
    }
    else
    {
        fprintf(stderr, "COULD NOT OPEN DIR: %s: %s\n", dir_, strerror(errno));
    }
    return i;
}

char **FindFilesInsideDir(const char *dir_, int *nf_)
{
    int nf = FindNumberFilesInsideDir(dir_);
    *nf_ = nf;
    char **files = new char*[nf];

    DIR *d;
    struct dirent *dir;
    d = opendir(dir_);
    int i = 0;
    if (d != NULL)
    {
        while ((dir = readdir(d)) != NULL)
        {
            int r1 = strcmp(".", dir->d_name);
            int r2 = strcmp("..", dir->d_name);
            if (r1 != 0 && r2 != 0)
            {
                int n = snprintf(NULL, 0, "%s/%s", dir_, dir->d_name) + 1;
                files[i] = new char[n];
                snprintf(files[i], n, "%s/%s", dir_, dir->d_name);
                i++;
            }
        }
        closedir(d);
    }
    else
    {
        fprintf(stderr, "COULD NOT OPEN DIR: %s: %s\n", dir_, strerror(errno));
    }
    return files;
}

void ParseFilesInsideDir(const char *dir, int *nfiles, int **qnt, char ****parsed)
{
    char **files = FindFilesInsideDir(dir, nfiles);
    char **stringsF = new char*[*nfiles];

    *qnt = new int[*nfiles];
    *parsed = new char**[*nfiles];

    for (int i = 0; i < *nfiles; i++)
    {
        FILE *f = fopen(files[i], "rb");
        stringsF[i] = ReadFile(f);
        fclose(f);
        (*parsed)[i] = Parse(stringsF[i], &(*qnt)[i]);
    }
}


int InitParticles(Particle **p)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/particles", &nfiles, &qnt, &parsed);
    int nP = 0;

    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Positions", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 2)
        {
            nP++;
        }
    }
    *p = new Particle[nP];

    int j_ = 0;
    for (int i = 0; i < nfiles; i++)
    {
        double U0 = GetValueDouble("U0", parsed[i], qnt[i]);
        double betadamp = GetValueDouble("BETADAMP", parsed[i], qnt[i]);
        int ip = GetIndexOfTag("Positions", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 2)
        {
            double x = strtod(parsed[i][j], NULL);
            double y = strtod(parsed[i][j + 1], NULL);
            (*p)[j_] = Particle(betadamp, U0, x, y);
            j_++;
        }
    }
    return nP;
}

int InitLines(LineSegment **l)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/lines", &nfiles, &qnt, &parsed);
    int nL = 0;

    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Positions", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 4)
        {
            nL++;
        }
    }
    *l = new LineSegment[nL];

    int j_ = 0;
    for (int i = 0; i < nfiles; i++)
    {
        double U0 = GetValueDouble("Potential", parsed[i], qnt[i]);
        double R = GetValueDouble("Radius", parsed[i], qnt[i]);
        int ip = GetIndexOfTag("Positions", parsed[i], qnt[i]);
        int start = 0;
        if (strlen(parsed[i][ip + 1]) == 1)
        {
            start = ip + 2;
        }
        else
        {
            start = ip + 1;
        }
        for (int j = start; j < qnt[i]; j += 4)
        {
            double x = strtod(parsed[i][j], NULL);
            double y = strtod(parsed[i][j + 1], NULL);
            double x1 = strtod(parsed[i][j + 2], NULL);
            double y1 = strtod(parsed[i][j + 3], NULL);
            (*l)[j_] = LineSegment(x, y, x1, y1, R, U0);
            j_++;
        }
    }
    return nL;
}


Matrix<Box> CreateBoxes(double range, size_t nPar, double Lx, double Ly)
{
    if (range > Lx || range > Ly)
    {
        fprintf(stderr, "Range Bigger Than Box\n");
        exit(1);
    }

    size_t nx = floor(Lx / range);
    size_t ny = floor(Ly / range);
    double lx = Lx / (double)nx;
    double ly = Ly / (double)ny;

    Matrix<Box>out(nx, ny);
    for (size_t i = 0; i < nx; ++i)
    {
        double x = lx * i;
        for (size_t j = 0; j < ny; ++j)
        {
            double y = ly * j;
            Box aux(lx, ly, x, y, nPar);
            out.SetData(&aux, j, i);
        }
    }
    return out;
}

template <typename E>
void AttBoxesOld(const size_t &nPar, const E* p, Matrix<Box> *out)
{
    for (size_t i = 0; i < out->nCols; ++i)
    {
        for (size_t j = 0; j < out->nRows; ++j)
        {
            (*out)(j, i).AttBox(p, nPar);
        }
    }
}

template <typename E>
void AttBoxes(const size_t &nPar, const E* p, Matrix<Box> *out)
{
    for (size_t i = 0; i < out->nCols; ++i)
    {
        for (size_t j = 0; j < out->nRows; ++j)
        {
            (*out)(j, i).SetIn(0);
        }
    }

    for (size_t i = 0; i < nPar; ++i)
    {
        int iX = p[i].x / (*out)(0, 0).GetLx();
        int iY = p[i].y / (*out)(0, 0).GetLy();
        int ind = (*out)(iY, iX).GetIn();
        (*out)(iY, iX).SetIndex(ind, i);
        ind++;
        (*out)(iY, iX).SetIn(ind);
    }
}
#endif //HEADERGUARD