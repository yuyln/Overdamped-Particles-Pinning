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
#include <simulator.h>

int main()
{
    Simulator s;
    //TODO: mena velocity out and save/load system
    s.FixCurrent(2.0);
    Integration(s);
    return 0;
}