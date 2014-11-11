#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System *system, double dt);
    void move(System *system, double dt);
    bool m_firstStep;
public:
    VelocityVerlet();
    ~VelocityVerlet();
    virtual void integrate(System *system, double dt);
};
