
#include <math.h>
#include "droplet.hpp"

int Droplet::num_droplets = 0;

Droplet::Droplet() {
    _multi = 0;
    _rcubed = _solute = _density = 0.;

    _mass = _volume = 0.;

    num_droplets++;
}

Droplet::Droplet(long multi, double rcubed, double solute,
                 double density) {
    _multi = multi;
    _rcubed = rcubed;
    _solute = solute;
    _density = density;

    _volume = (4.*M_PI/3.)*_rcubed;
    _mass = _volume*_density;

    num_droplets++;
}

Droplet::Droplet(const Droplet &d) {
    _multi = d._multi;
    _rcubed = d._rcubed;
    _solute = d._solute;
    _density = d._density;

    _volume = d._volume;
    _mass = d._mass;

    num_droplets++;
}

Droplet::~Droplet() {
    --num_droplets;
}

double Droplet::get_mass() const {
    return _mass;
}

long Droplet::get_multi() const {
    return _multi;
}

double Droplet::get_radius() const {
    return pow(_rcubed, 1./3.);
}

double Droplet::get_solute() const {
    return _solute;
}

double Droplet::get_volume() const {
    return _volume;
}

double Droplet::get_terminal_v() const {
    // BEARD, 1976
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

Droplet & Droplet::operator= (const Droplet & other) {
    if (this == &other)
        return *this;

    _multi = other._multi;
    _rcubed = other._rcubed;
    _solute = other._solute;
    _density = other._density;

    _volume = other._volume;
    _mass = other._mass;

    return *this;
}

bool operator< (const Droplet & d1, const Droplet & d2) {
    return (d1._multi < d2._multi);
}

// Basic string output
std::ostream &operator<<(std::ostream &os, const Droplet &d) {
    os << "Droplet( "
       << "multi=" << d._multi << ", "
       << "radius=" << d.get_radius()*1e6 << " micron"
       << " )";
    return os;
}
