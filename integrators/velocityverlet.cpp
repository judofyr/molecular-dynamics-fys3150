#include <integrators/velocityverlet.h>
#include <system.h>

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::firstKick(System *system, double dt)
{
    m_firstStep = false;
    system->calculateForces();
    halfKick(system, dt);
}

void VelocityVerlet::halfKick(System *system, double dt)
{
    for (auto &atom : system->atoms()) {
        atom->velocity.addAndMultiply(atom->force, dt/2/atom->mass());
    }
}

void VelocityVerlet::move(System *system, double dt)
{
    for (auto &atom : system->atoms()) {
        atom->position.addAndMultiply(atom->velocity, dt);
    }
}

void VelocityVerlet::integrate(System *system, double dt)
{
    if(m_firstStep) {
        firstKick(system, dt);
    } else halfKick(system, dt);
    move(system, dt);
    system->applyPeriodicBoundaryConditions();
    system->calculateForces();
    halfKick(system, dt);
}
