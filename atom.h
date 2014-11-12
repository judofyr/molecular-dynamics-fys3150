#pragma once
#include <config.h>
#include <math/vec3.h>
#include <vector>
#include <functional>
#include <algorithm>

using CompPhys::vec3;

struct AtomVec {
    float x[MD_BLOCKSIZE];
    float y[MD_BLOCKSIZE];
    float z[MD_BLOCKSIZE];

    CompPhys::vec3 vec3(int i) {
        CompPhys::vec3 vec(x[i], y[i], z[i]);
        return vec;
    }
};

enum class AtomBlockType { REAL, GHOST };

struct AtomBlock {
    short count = 0;
    AtomBlockType type = AtomBlockType::REAL;
    AtomVec position;
    AtomVec force;
    AtomVec *velocity;
    std::vector<AtomBlock *> neighbours;

    void addNeighbour(AtomBlock *block) {
        for (auto &otherBlock : neighbours) {
            if (block == otherBlock) {
                return;
            }
        }

        neighbours.push_back(block);
    }

    void clearNeighbours() {
        neighbours.clear();
    }

    short neighbourCount() {
        return neighbours.size();
    }

    void each_neighbour(std::function<void(AtomBlock *block)> action) {
        for (auto &block : neighbours) {
            action(block);
        }
    }
};

struct AtomRef {
    AtomBlock *block;
    unsigned char idx;
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
