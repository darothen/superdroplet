
#include <math.h>
#include <tuple>

#include "droplet.hpp"

inline double dmax(double a, double b) {
    return a >= b ? a : b;
}

inline long lmax(long a, long b) {
    return a >= b ? a : b;
}

inline double dmin(double a, double b) {
    return a <= b ? a : b;
}

inline int ifloor(long x) {
    return int(floor(x));
}

double exp_dist_moments(double x, double n0, double x0, double l) {
    return pow(x, l)*(n0/x0)*exp(-x/x0);
}

std::tuple<int, int> s_to_min_s(double seconds) {
    double seconds_over = fmod(seconds, 60.);
    return std::make_tuple( (seconds - seconds_over)/60., seconds_over );
}


// Compare two droplets to see which is smaller
bool smaller(const Droplet & d1, const Droplet & d2) {
    return (d1.get_radius() < d2.get_radius());
}

double total_water(const Droplet * droplets, int n) {

    double wm0 = 0.0;

    for (int i=0; i < n; i++) {
        Droplet d = droplets[i];
        wm0 += d.get_mass()*d.get_multi()/1e3;
    }
    return wm0;
}

