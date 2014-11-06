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
    double dtHalfOverMass = dt * 0.5 / system->atom()->mass();

    system->for_each([dtHalfOverMass](AtomBlock &block, int i) {
        block.velocity->x[i] += block.force.x[i] * dtHalfOverMass;
        block.velocity->y[i] += block.force.y[i] * dtHalfOverMass;
        block.velocity->z[i] += block.force.z[i] * dtHalfOverMass;
    });
}

void VelocityVerlet::move(System *system, double dt)
{
    system->for_each([dt](AtomBlock &block, int i) {
        block.position.x[i] += block.velocity->x[i] * dt;
        block.position.y[i] += block.velocity->y[i] * dt;
        block.position.z[i] += block.velocity->z[i] * dt;
    });
}

void VelocityVerlet::integrate(System *system, double dt)
{
    if(m_firstStep) {
        firstKick(system, dt);
    } else halfKick(system, dt);
    move(system, dt);
    system->calculateForces();
    halfKick(system, dt);
}
