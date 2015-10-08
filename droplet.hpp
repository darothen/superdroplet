
#ifndef DROPLET_H_
#define DROPLET_H_

class Droplet
{
private:
    long _multi;
    double _rcubed;
    double _solute;
    double _density;

// TODO: superdroplet_count?

public:
    Droplet();
    Droplet(long multi, double rcubed, double solute, double density=1000.0);

    double get_mass() const;
    long   get_multi() const;
    double get_radius() const;
    double get_terminal_v() const;
    double get_volume() const;

    friend bool operator< (const Droplet & d1, const Droplet & d2);

};

#endif