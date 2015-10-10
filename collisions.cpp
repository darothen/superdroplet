
#include <algorithm>
#include <iostream>
#include <math.h>

#include "collisions.hpp"
#include "droplet.hpp"
#include "util.hpp"

using namespace std;

void collision_step(Droplet * droplets, int n_part, double t_c, double delta_V) {

    cout << "PREP STEPS" << endl;

    // 1) Make the random permutation of the droplet list
    cout << "   SHUFFLE LIST" << endl;
    random_shuffle(droplets, droplets+n_part);

    // 2) Make the candidate pairs
    cout << "   GEN PAIRS" << endl;


    // 3) Generate the uniform random numbers
    cout << "PROBABILITY LOOP" << endl;
    double scaling = (n_part*(n_part-1)/2.)/floor(n_part/2.);

    cout << "PROB / COLLISION LOOP" << endl;
    int counter = 0;
    int half_n_part = n_part/2;

    for (int i=0; i < half_n_part; i++) {
        Droplet sd_j = droplets[i];
        Droplet sd_k = droplets[i + half_n_part];

        double prob = urand();

    }

}