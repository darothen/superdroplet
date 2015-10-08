
#include <iostream>
#include <math.h>
#include <vector>

#include "droplet.hpp"
#include "util.hpp"

int main()
{

    using std::vector;
    using std::cout;
    using std::endl;

    double delta_V = 1e6; // Cell volume, m^3
    double t_c     = 1.0; // timestep, seconds
    int n_part     = pow(2, 11); // # superdroplets to use

    // CASE SETTINGS
    int t_end = 60*60+1;
    int plot_dt = 30*60;
    double m_tot_ana = 1.0;
    double f = 1.0;
    double R_0 = 10.e-6/f;
    double X_0 = (4.*M_PI/3.)*(pow(R_0, 3.));
    double M_0 = X_0*1000.0;

    /* ************************************************ */

    vector<Droplet> droplets;



}