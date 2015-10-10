
#ifndef UTIL_H_
#define UTIL_H_

#include "droplet.hpp"

inline double dmax(double, double);
inline long lmax(long, long);
inline double dmin(double, double);
inline int ifloor(long);

double exp_dist_moments(double x, double n0, double x0, double l=0.);
std::tuple<int, int> s_to_min_s(double seconds);
bool smaller(const Droplet &, const Droplet &);
double urand(void);
double total_water(const Droplet *, int);

#endif