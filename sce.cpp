
#include <boost/random.hpp>
#include <algorithm>
#include <iostream>
#include <math.h>

#include "constants.hpp"
#include "droplet.hpp"
#include "util.hpp"
#include "collisions.hpp"

using namespace std;
using namespace constants;

const bool DEBUG = false;

int main() {

    double delta_V = 1e6; // Cell volume, m^3
    double t_c = 1.0; // timestep, seconds
    const resolution n_part = lo;

    // CASE SETTINGS
    int t_end = 3601;
    int plot_dt = 1200;
    double n_0 = pow(2., 23);
    double R_0 = 30.531e-6;
    double X_0 = (4. * M_PI / 3.) * (R_0 * R_0 * R_0);
    double M_0 = X_0 * RHO_WATER;
    double m_tot_ana = 1.0;

    /* ************************************************ */

    // Create a droplet array
    Droplet droplets[n_part];

    // Pre-compute the super-droplet multiplicity
    double total_droplets = delta_V * n_0;
    long xi_i = long(floor(total_droplets / n_part));

    // Using an exponential distribution, sample droplet masses and initialize
    // the droplets.
    boost::random::mt19937 rng;
    boost::random::exponential_distribution<> dist(1. / X_0);
    for (int i = 0; i < n_part; i++) {
        double x = dist(rng);
        double r = pow(x * 3. / M_PI / 4., 1. / 3.);
        droplets[i] = Droplet(xi_i, pow(r, 3.));
//        cout << droplets[i] << "\n";
    }
    sort(droplets, droplets + n_part, smaller); // Sort using a comparator function

    cout << "GRID SETUP" << endl;
    cout << "   radii: " << droplets[0].get_radius() << " - "
    << droplets[n_part - 1].get_radius() << " m" << endl;
    cout << "  volume: " << droplets[0].get_volume() << " - "
    << droplets[n_part - 1].get_volume() << " m^3" << endl;
    cout << "    mass: " << droplets[0].get_mass() << " - "
    << droplets[n_part - 1].get_mass() << " kg" << endl;

    cout << "SD SETUP" << endl;
    cout << "   N_s: " << n_part << endl;
    cout << "  xi_i: " << xi_i << endl;
    long N_per_SD = total_droplets / xi_i / n_part;
    cout << " N per SD_xi: " << N_per_SD << endl;
    cout << "(Initialized " << Droplet::global_droplet_count() << " droplets)"
         << endl;

    // ----------------------------------------------------

    cout << "\nBEGINNING MAIN ROUTINE\n" << endl;

    double wm0 = total_water(droplets, n_part);
    cout << "Initial water mass = " << wm0 << " kg" << endl;

    double t = 0.;
    int ti = 0;

    while (t < t_end) {

        // Get/print timestep info
        auto min_sec = s_to_min_s(t);
        int minutes = std::get<0>(min_sec);
        int seconds = std::get<1>(min_sec);

        cout << "\nSTEP " << ti
             << " (" << minutes << " min";
        if (seconds > 0)
            cout << " " << seconds << " sec";
        cout << ")" << endl;

        // --------------------------------
        collision_step(droplets, n_part, t_c, delta_V);
        // --------------------------------

        ti++;
        t = ti*t_c;

        cout <<  Droplet::global_droplet_count() << " droplets remain" << endl;


        if (DEBUG) break;
    }

}