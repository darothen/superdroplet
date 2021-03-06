
#ifndef DROPLET_H_
#define DROPLET_H_

#include <iostream>

#include "constants.hpp"

using constants::RHO_WATER;

struct Droplet {
private:
    static int num_droplets;

public:
    long _multi;
    double _rcubed, _solute, _density;
    double _radius, _volume, _mass;

    // Constructors
    Droplet();
    Droplet(long multi, double rcubed, double solute=0.0, double density=RHO_WATER);
    Droplet(const Droplet &d);
    ~Droplet();

    // Accessors / mutators
//    double get_mass() const;
//    long   get_multi() const;
//    double get_radius() const;
//    double get_solute() const;
    double get_terminal_v() const;
    void set_droplet(long multi, double rcubed, double solute);
//    double get_volume() const;

    // Static methods
    static int global_droplet_count() { return num_droplets; }

    // Operator overloads and friends
    Droplet & operator= (const Droplet &);
    friend bool operator< (const Droplet & d1, const Droplet & d2);
    friend std::ostream & operator<<(std::ostream & os, const Droplet & d);

};

#endif