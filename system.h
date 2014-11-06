#pragma once
#include <vector>
#include <atom.h>
#include <math/vec3.h>
#include <functional>

class Potential; class Integrator;
using std::vector;
using CompPhys::vec3;

class System
{
private:
    vec3 m_systemSize;
    Atom *m_atom;
    vector<AtomVec> m_velocities;
    Potential *m_potential;
    Integrator *m_integrator;
    double m_rCutOff;
    double m_currentTime;
    int m_steps;
    void createBlock();

public:
    // TODO: Either make getters, or find a better API that doesn't expose these details
    vector<AtomBlock> m_ghostBlocks;
    vector<AtomBlock> m_atomBlocks;

    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant);
    void applyPeriodicGhostBlocks();
    void applyPeriodicBoundaryConditions();
    size_t atomCount();
    void removeMomentum();
    void calculateForces();
    void step(double dt);
    void for_each(std::function<void(AtomBlock &block, int atomIdx)>);

    // Setters and getters
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    void setCurrentTime(double currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    Atom *atom() const  { return m_atom;  }
    void setAtom(Atom *atom) { m_atom = atom; }

    double rCutOff() const { return m_rCutOff; }
    void setRCutOff(double rCutOff) { m_rCutOff = rCutOff; }

private:
    void addAtom(double x, double y, double z);
};
