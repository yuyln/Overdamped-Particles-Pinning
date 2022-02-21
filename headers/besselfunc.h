#ifndef __BESSELFUNC
#define __BESSELFUNC
#include <cstdlib>
#include <cmath>

double BESSI0(double x)
{
    double y, ax, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7, q8, q9, bx;
    p1 = 1.e0;
    p2 = 3.5156229e0;
    p3 = 3.0899424e0;
    p4 = 1.2067429e0;
    p5 = 0.2659732e0;
    p6 = 0.360768e-1;
    p7 = 0.45813e-2;
    q1 = 0.39894228e0;
    q2 = 0.1328592e-1;
    q3 = 0.225319e-2;
    q4 = -0.157565e-2;
    q5 = 0.916281e-2;
    q6 = -0.2057706e-1;
    q7 = 0.2635537e-1;
    q8 = -0.1647633e-1;
    q9 = 0.392377e-2;
    if (fabs(x) < 3.75)
    {
        y = (x / 3.75) * (x / 3.75);
        return p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))));
    }
    else
    {
        ax = fabs(x);
        y = 3.75 / x;
        bx = exp(ax) / sqrt(ax);
        ax = q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9)))))));
        return ax * bx;
    }
}

double BESSK0(double x)
{
    double y, ax, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7;
    p1 = -0.57721566e0;
    p2 = 0.42278420e0;
    p3 = 0.23069756e0;
    p4 = 0.3488590e-1;
    p5 = 0.262698e-2;
    p6 = 0.10750e-3;
    p7 = 0.74e-5;
    q1 = 1.25331414e0;
    q2 = -0.7832358e-1;
    q3 = 0.2189568e-1;
    q4 = -0.1062446e-1;
    q5 = 0.587872e-2;
    q6 = -0.251540e-2;
    q7 = 0.53208e-3;
    if (x == 0)
    {
        return INFINITY;
    }
    if (x <= 2.0)
    {
        y = x * x / 4.0;
        ax = -log(x / 2.0) * BESSI0(x);
        return ax + (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else
    {
        y = 2.0 / x;
        ax = exp(-x) / sqrt(x);
        return ax * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
}

double BESSI1(double x)
{
    double y, ax, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7, q8, q9, bx;
    p1 = 0.5e0;
    p2 = 0.87890594e0;
    p3 = 0.51498869e0;
    p4 = 0.15084934e0;
    p5 = 0.2658733e-1;
    p6 = 0.301532e-2;
    p7 = 0.32411e-3;
    q1 = 0.39894228e0;
    q2 = -0.3988024e-1;
    q3 = -0.362018e-2;
    q4 = 0.163801e-2;
    q5 = -0.1031555e-1;
    q6 = 0.2282967e-1;
    q7 = -0.2895312e-1;
    q8 = 0.1787654e-1;
    q9 = -0.420059e-2;
    if (fabs(x) < 3.75)
    {
        y = (x / 3.75) * (x / 3.75);
        return x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else
    {
        ax = fabs(x);
        y = 3.75 / ax;
        bx = exp(ax) / sqrt(ax);
        ax = q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9)))))));
        return ax * bx;
    }
}

double BESSK1(double x)
{
    double y, ax, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7;
    p1 = 1.e0;
    p2 = 0.15443144e0;
    p3 = -0.67278579e0;
    p4 = -0.18156897e0;
    p5 = -0.1919402e-1;
    p6 = -0.110404e-2;
    p7 = -0.4686e-4;
    q1 = 1.25331414e0;
    q2 = 0.23498619e0;
    q3 = -0.3655620e-1;
    q4 = 0.1504268e-1;
    q5 = -0.780353e-2;
    q6 = 0.325614e-2;
    q7 = -0.68245e-3;
    if (x == 0)
    {
        return INFINITY;
    }
    if (x <= 2.0)
    {
        y = x * x / 4.0;
        ax = log(x / 2.0) * BESSI1(x);
        return ax + (1.e0 / x) * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else
    {
        y = 2.0 / x;
        ax = exp(-x) / sqrt(x);
        return ax * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
}

double BESSK(int n, double x)
{
    double TOX, BK, BKM, BKP;
    if (n == 0)
    {
        return BESSK0(x);
    }
    if (n == 1)
    {
        return BESSK1(x);
    }
    if (x == 0.0)
    {
        return INFINITY;
    }
    TOX = 2.0 / x;
    BK = BESSK1(x);
    BKM = BESSK0(x);
    for (int j = 1; j <= n - 1; j++)
    {
        BKP = BKM + (double)j * TOX * BK;
        BKM = BK;
        BK = BKP;
    }
    return BK;
}
#endif