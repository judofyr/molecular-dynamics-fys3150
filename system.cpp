#include <system.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <unitconverter.h>

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0)
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    for (auto &atom : m_atoms) {
        for (int i = 0; i < 3; i++) {
            auto pos = &atom->position[i];
            double size = m_systemSize[i];
            if (*pos < -size*0.5) *pos += size;
            if (*pos >= size*0.5) *pos -= size;
        }
    }
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeMomentum() {
    vec3 netMomentum;
    for (auto &atom : m_atoms) {
        netMomentum.addAndMultiply(atom->velocity, atom->mass());
    }

    // We want to remove this much momentum from each atom
    vec3 perMomentum = netMomentum / m_atoms.size();

    for (auto &atom : m_atoms) {
        atom->velocity.addAndMultiply(perMomentum, 1/atom->mass());
    }
}

void System::resetForcesOnAllAtoms() {

}

void System::addAtom(double x, double y, double z) {
    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); // Argon mass
    atom->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
    atom->position.set(x, y, z);
    m_atoms.push_back(atom);
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant) {
    double b = latticeConstant;
    for (int x_i = 0; x_i < numberOfUnitCellsEachDimension; x_i++) {
    for (int y_i = 0; y_i < numberOfUnitCellsEachDimension; y_i++) {
    for (int z_i = 0; z_i < numberOfUnitCellsEachDimension; z_i++) {
        addAtom(b*x_i + 0,     b*y_i + 0,     b*z_i + 0);
        addAtom(b*x_i + (b/2), b*y_i + (b/2), b*z_i + 0);
        addAtom(b*x_i + 0,     b*y_i + (b/2), b*z_i + (b/2));
        addAtom(b*x_i + (b/2), b*y_i + 0,     b*z_i + (b/2));

    }}}
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
