
#include <boost/random.hpp>

#include "droplet.hpp"

double exp_dist_moments(double x, double n0, double x0, double l) {
    return pow(x, l)*(n0/x0)*exp(-x/x0);
}

std::tuple<int, int> s_to_min_s(double seconds) {
    double seconds_over = fmod(seconds, 60.);
    return std::make_tuple( (seconds - seconds_over)/60., seconds_over );
}


// Compare two droplets to see which is smaller
bool smaller(const Droplet & d1, const Droplet & d2) {
    return (d1._radius < d2._radius);
}

// Generate a random number b/t [0, 1]
double urand(void) {
    // http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
    boost::mt19937 rng;
    static boost::uniform_01<boost::mt19937> zeroone(rng);
    return zeroone();
}

double total_water(const Droplet * droplets, int n) {

    double wm0 = 0.0;

    for (int i=0; i < n; i++) {
        Droplet d = droplets[i];
        wm0 += d._mass*d._multi/1e3;
    }
    return wm0;
}

