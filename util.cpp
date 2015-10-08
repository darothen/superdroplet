
#include <cmath.h>

inline double dmax(double a, double b)
{
    return a >= b ? a : b;
}

inline long lmax(long a, long b)
{
    return a >= b ? a : b;
}

inline double dmin(double a, double b)
{
    return a <= b ? a : b;
}

inline int ifloor(long x)
{
    return int(floor(x));
}

double exp_dist_moments(double x, double n0, double x0, double l)
{
    return pow(x, l)*(n0/x0)*exp(-x/x0);
}
