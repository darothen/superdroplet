
#include <math.h>
#include "droplet.hpp"


Droplet::Droplet()
{
    _multi = 0;
    _rcubed = _solute = _density = 0.;
}

Droplet::Droplet(long multi, double rcubed, double solute,
                 double density) {
    _multi = multi;
    _rcubed = rcubed;
    _solute = solute;
    _density = density;
}

double Droplet::get_mass() const {
    return _density *this->get_volume();
}

long Droplet::get_multi() const {
    return _multi;
}

double Droplet::get_radius() const {
    return pow(_rcubed, 1./3.);
}

double Droplet::get_terminal_v() const {
    // BEARD, 1976
    double rcubed = _rcubed;
    double mass = this->get_mass();

    double r = this->get_radius();
    double d = 2.*r*1e6; // diameter, m -> micron
    double x = mass * 1e3; // convert kg -> g

    double alpha, x_to_beta;

    if (d <= 134.43)
    {
        alpha = 4.5795e5;
        x_to_beta = pow(x, 2./3.);
    }
    else if (134.43 < d <= 1511.64)
    {
        alpha = 4962.0;
        x_to_beta = pow(x, 1./3.);
    }
    else if (1511.64 < d <= 3477.84)
    {
        alpha = 1732.0;
        x_to_beta = pow(x, 1./6.);
    }
    else
    {
        alpha = 917.0;
        x_to_beta = 1.0;
    }

    return (1e-2 * alpha * x_to_beta); // from cm/s -> m/s

}

double Droplet::get_volume() const {
    return _rcubed *4.*M_PI/3.;
}

bool operator< (const Droplet & d1, const Droplet & d2) {
    return (d1.get_multi() <= d2.get_multi());
}