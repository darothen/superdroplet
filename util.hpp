
#ifndef UTIL_H_
#define UTIL_H_

#include <math.h>

#include "droplet.hpp"

inline double dmax(double a, double b) { return a >= b ? a : b; };
inline long lmax(long a, long b) { return a >= b ? a : b; };
inline double dmin(double a, double b) { return a <= b ? a : b; };
inline int ifloor(long i) { return (int) floor(i); };

double exp_dist_moments(double x, double n0, double x0, double l=0.);
std::tuple<int, int> s_to_min_s(double seconds);
bool smaller(const Droplet &, const Droplet &);
double urand(void);
double total_water(const Droplet *, int);

#endif