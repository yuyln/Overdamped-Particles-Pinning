#ifndef __FUNCS
#define __FUNCS
#include <dirent.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>

int FindNumberFilesInsideDir(const char *dir);
char **FindFilesInsideDir(const char *dir, int *nf);


#endif

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
                files[i] = dir->d_name;
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