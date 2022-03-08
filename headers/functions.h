#ifndef __FUNCS
#define __FUNCS
#include <dirent.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <pinning.h>
#include <parser.h>
#include <particle.h>

int FindNumberFilesInsideDir(const char *dir);
char **FindFilesInsideDir(const char *dir, int *nf);
void ParseFilesInsideDir(const char *dir, int *nfiles, int **qnt, char ****parsed);
int InitPinnings(Pinning **p);
int InitParticles(Particle **p);


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
                int n = snprintf(NULL, 0, "%s/%s", dir_, dir->d_name);
                files[i] = new char[n];
                snprintf(files[i], n + 1, "%s/%s", dir_, dir->d_name);
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
        ReadFile(files[i], &stringsF[i]);
        (*parsed)[i] = Parse(stringsF[i], &(*qnt)[i]);
    }
}

int InitPinnings(Pinning **p)
{
    int nfiles;
    int *qnt;
    char ***parsed;
    ParseFilesInsideDir("./input/pinnings", &nfiles, &qnt, &parsed);
    int nPi = 0;

    for (int i = 0; i < nfiles; i++)
    {
        int ip = GetIndexOfTag("Positions", parsed[i], qnt[i]);
        for (int j = ip + 1; j < qnt[i]; j += 2)
        {
            nPi++;
        }
    }
    *p = new Pinning[nPi];

    int j_ = 0;
    for (int i = 0; i < nfiles; i++)
    {
        double U0 = GetValueDouble("U0", parsed[i], qnt[i]);
        double R0 = GetValueDouble("R0", parsed[i], qnt[i]);
        int ip = GetIndexOfTag("Positions", parsed[i], qnt[i]);
        for (int j = ip + 1; j < qnt[i]; j += 2)
        {
            double x = strtod(parsed[i][j], NULL);
            double y = strtod(parsed[i][j + 1], NULL);
            (*p)[j_] = InitPinning(x, y, U0, R0);
            j_++;
        }
    }
    return nPi;
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
        for (int j = ip + 1; j < qnt[i]; j += 2)
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
        for (int j = ip + 1; j < qnt[i]; j += 2)
        {
            double x = strtod(parsed[i][j], NULL);
            double y = strtod(parsed[i][j + 1], NULL);
            (*p)[j_] = InitParticle(betadamp, U0, x, y);
            j_++;
        }
    }
    return nP;
}


#endif //HEADERGUARD