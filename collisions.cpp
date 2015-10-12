
#include <algorithm>
#include <iostream>
#include <math.h>

#include "collisions.hpp"
#include "util.hpp"

using namespace std;


inline double kernel(const Droplet &sd_j, const Droplet &sd_k) {
    // TODO: implement alternative collision kernels
    const double golovin_b = 1.5e3;
    return golovin_b*(sd_j._rcubed + sd_k._rcubed)*4.*M_PI/3.;
}

void multi_coalesce(Droplet &sd_j, Droplet &sd_k, double gamma) {
//    cout << ">>>" << endl;
//    cout << gamma << " " << sd_j << " " << sd_k << endl;

    double gamma_tilde = dmin(gamma, (double) floor(sd_j._multi / sd_k._multi));

    long excess = sd_j._multi - (long) floor(gamma_tilde*sd_k._multi);
//    cout << "Excess = " << excess << endl;

    long multi_j_p, multi_k_p;
    double rcubed_j_p, rcubed_k_p;
    double solute_j_p, solute_k_p;

    if (excess > 0) {
        multi_j_p = excess;
        multi_k_p = sd_k._multi;

        rcubed_j_p =sd_j._rcubed;
        rcubed_k_p = gamma_tilde*rcubed_j_p + sd_k._rcubed;

        solute_j_p = sd_j._solute;
        solute_k_p = gamma_tilde*solute_j_p + sd_k._solute;

        sd_j = Droplet(multi_j_p, rcubed_j_p, solute_j_p);
        sd_k = Droplet(multi_k_p, rcubed_k_p, solute_k_p);
    } else {
        multi_j_p = (long) floor(sd_k._multi / 2);
        multi_k_p = sd_k._multi - multi_j_p;

        rcubed_j_p = rcubed_k_p =
            gamma_tilde*sd_j._rcubed + sd_k._rcubed;

        solute_j_p = solute_k_p =
            gamma_tilde*sd_j._solute + sd_k._solute;

        sd_j = Droplet(multi_k_p, rcubed_j_p, solute_j_p);
        sd_k = Droplet(multi_j_p, rcubed_k_p, solute_k_p);
    }

//    cout << "---" << endl;
//    cout << sd_j << " " << sd_k << endl;
//    cout << "<<<" << endl;
}

void collision_step(std::vector<Droplet> & droplets, double t_c, double delta_V) {

    cout << "PREP STEPS" << endl;

    // 1) Make the random permutation of the droplet list
    cout << "   SHUFFLE LIST" << endl;
    random_shuffle(droplets.begin(), droplets.end());

    // 2) Make the candidate pairs
    cout << "   GEN PAIRS" << endl;

    // 3) Generate the uniform random numbers
    cout << "PROBABILITY LOOP" << endl;
    long n_part = droplets.size();
    double scaling = (n_part*(n_part-1)/2.)/floor(n_part/2.);

    cout << "PROB / COLLISION LOOP" << endl;
    unsigned int counter = 0;
    long half_n_part = n_part/2;

    unsigned int big_probs = 0;
    double max_prob = 0., min_prob = 1.0;

    for (int i=0; i < half_n_part; i++) {
        Droplet & sd_j = droplets[i];
        Droplet & sd_k = droplets[i + half_n_part];

        double phi = urand();

        double K_ij = kernel(sd_j, sd_k);
        long max_xi = max(sd_j._multi, sd_k._multi);

        double prob = scaling*max_xi*(t_c/delta_V)*K_ij;

        if (prob > max_prob)
            max_prob = prob;
        if (prob < min_prob)
            min_prob = prob;
        if (prob > 1)
            big_probs++;

        // Check for collision and coalescence if necessary
        if ( (prob - floor(prob)) >= phi ) {
            double gamma = floor(prob) + 1;

            if (sd_j < sd_k)
                multi_coalesce(sd_k, sd_j, gamma);
            else
                multi_coalesce(sd_j, sd_k, gamma);

            counter++;
        }

    }

    cout << counter << " collisions simulated" << endl;
    cout << " Max/min probabilities (count): "
         << min_prob << " "
         << max_prob << " "
         << big_probs << endl;

}