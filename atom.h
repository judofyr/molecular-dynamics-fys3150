#pragma once
#include <math/vec3.h>

using CompPhys::vec3;

#define ATOMBLOCKSIZE 16

struct AtomVec {
    float x[ATOMBLOCKSIZE];
    float y[ATOMBLOCKSIZE];
    float z[ATOMBLOCKSIZE];

    CompPhys::vec3 vec3(int i) {
        CompPhys::vec3 vec(x[i], y[i], z[i]);
        return vec;
    }
};

struct AtomBlock {
    int counter = 0;
    AtomVec position;
    AtomVec force;
    AtomVec *velocity;
};

class Atom
{
private:
    double m_mass;
public:
    Atom(double mass);
    ~Atom();
    vec3 velocityMaxwellian(double temperature);

    inline double mass() { return m_mass; }
    inline void setMass(double mass) { m_mass = mass; }
};
