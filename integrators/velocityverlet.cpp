#include <integrators/velocityverlet.h>
#include <system.h>
#include "cpelapsedtimer.h"

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::halfKick(System *system, double dt)
{
    CPElapsedTimer::getInstance().halfKick().start();
    double dtHalfOverMass = dt * 0.5 / system->atom()->mass();

    system->for_each([dtHalfOverMass](AtomBlock &block, int i) {
        block.velocity->x[i] += block.force.x[i] * dtHalfOverMass;
        block.velocity->y[i] += block.force.y[i] * dtHalfOverMass;
        block.velocity->z[i] += block.force.z[i] * dtHalfOverMass;
    });
    CPElapsedTimer::getInstance().halfKick().stop();
}

void VelocityVerlet::move(System *system, double dt)
{
    CPElapsedTimer::getInstance().move().start();

    system->for_each([dt](AtomBlock &block, int i) {
        block.position.x[i] += block.velocity->x[i] * dt;
        block.position.y[i] += block.velocity->y[i] * dt;
        block.position.z[i] += block.velocity->z[i] * dt;
    });

    system->for_each_ghost([dt](AtomBlock &block, int i, AtomRef &ref) {
        block.position.x[i] += ref.block->velocity->x[ref.idx] * dt;
        block.position.y[i] += ref.block->velocity->y[ref.idx] * dt;
        block.position.z[i] += ref.block->velocity->z[ref.idx] * dt;
    });

    CPElapsedTimer::getInstance().move().stop();
}

void VelocityVerlet::integrate(System *system, double dt)
{
    if(m_firstStep) {
        system->calculateForces();
        m_firstStep = false;
    }

    halfKick(system, dt);
    move(system, dt);
    system->calculateForces();
    halfKick(system, dt);
}
