
#include <boost/format.hpp>
#include <boost/random.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

#include "constants.hpp"
#include "droplet.hpp"
#include "util.hpp"
#include "collisions.hpp"

using namespace std;
using namespace constants;
using std::string;
using boost::format;

const bool DEBUG = false;

int main() {

    double delta_V = 1e6; // Cell volume, m^3
    double t_c = 1.0; // timestep, seconds
    const resolution n_part = med_lo; // lo (2^11), med_lo (2^13), med_hi (2^15), hi (2^17)

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
    std::vector<Droplet> droplets;

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
        droplets.push_back(Droplet(xi_i, pow(r, 3.)));
//        cout << droplets[i] << "\n";
    }
    std::sort(droplets.begin(), droplets.end(), smaller); // Sort using a comparator function

    cout << "GRID SETUP" << endl;
    cout << "   radii: " << droplets[0]._radius << " - "
    << droplets[n_part - 1]._radius << " m" << endl;
    cout << "  volume: " << droplets[0]._volume << " - "
    << droplets[n_part - 1]._volume << " m^3" << endl;
    cout << "    mass: " << droplets[0]._mass << " - "
    << droplets[n_part - 1]._mass << " kg" << endl;

    cout << "SD SETUP" << endl;
    cout << "   N_s: " << n_part << endl;
    cout << "  xi_i: " << xi_i << endl;
    long N_per_SD = total_droplets / xi_i / n_part;
    cout << " N per SD_xi: " << N_per_SD << endl;
    cout << "(Initialized " << Droplet::global_droplet_count() << " droplets)"
         << endl;

    // ----------------------------------------------------

    cout << "\nBEGINNING MAIN ROUTINE\n" << endl;

    double wm0 = total_water(droplets);
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

        if ( (floor(t/plot_dt) == (t/plot_dt)) || (ti == 0) ) {
            cout << endl << "Writing output... ";
            string out_fn = str(format("%5.1fn_output.txt") % t);
            cout << "(" << out_fn << ")" << endl;
            std::ofstream out_file(out_fn);
            for (Droplet &d : droplets) {
                out_file << d._rcubed << "," << d._multi << endl;
            }

            cout << " done." << endl;
        }

        // --------------------------------
        collision_step(droplets, t_c, delta_V);
        // --------------------------------

        ti++;
        t = (int) ti*t_c;

        cout <<  Droplet::global_droplet_count() << " droplets remain" << endl;

        if (DEBUG) break;
    }
}