#pragma once
#include <math/vec3.h>
#include <vector>
#include <functional>
#include <algorithm>

using CompPhys::vec3;

#define ATOMBLOCKSIZE 16
#define ATOMBLOCKNEIGHBOURINLINE 32

struct AtomVec {
    float x[ATOMBLOCKSIZE];
    float y[ATOMBLOCKSIZE];
    float z[ATOMBLOCKSIZE];

    CompPhys::vec3 vec3(int i) {
        CompPhys::vec3 vec(x[i], y[i], z[i]);
        return vec;
    }
};

enum class AtomBlockType { REAL, GHOST };

struct AtomBlock {
    short counter = 0;
    short inlineNeighbourCount = 0;
    AtomBlockType type = AtomBlockType::REAL;
    AtomVec position;
    AtomVec force;
    AtomVec *velocity;
    AtomBlock *neighbours[ATOMBLOCKNEIGHBOURINLINE];
    std::vector<AtomBlock *> moreNeighbours;

    void addNeighbour(AtomBlock *block) {
        for (int i = 0; i < inlineNeighbourCount; i++) {
            if (neighbours[i] == block) {
                // Found it inline
                return;
            }
        }

        if (inlineNeighbourCount < ATOMBLOCKNEIGHBOURINLINE) {
            neighbours[inlineNeighbourCount++] = block;
            return;
        }

        for (auto &otherBlock : moreNeighbours) {
            if (block == otherBlock) {
                return;
            }
        }

        moreNeighbours.push_back(block);
    }

    void clearNeighbours() {
        inlineNeighbourCount = 0;
        moreNeighbours.clear();
    }

    short neighbourCount() {
        return inlineNeighbourCount + moreNeighbours.size();
    }

    void each_neighbour(std::function<void(AtomBlock *block)> action) {
        for (int i = 0; i < inlineNeighbourCount; i++) {
            action(neighbours[i]);
        }

        for (auto &block : moreNeighbours) {
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
