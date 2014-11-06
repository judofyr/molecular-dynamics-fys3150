#include <statisticssampler.h>


double StatisticsSampler::kineticEnergy()
{
    return m_kineticEnergy;
}

double StatisticsSampler::instantaneousTemperature()
{
    return (2.0/3) * m_kineticEnergy / m_numAtoms;
}


double StatisticsSampler::dt()
{
    return m_dt;
}

StatisticsSampler::StatisticsSampler()
{

}

StatisticsSampler::~StatisticsSampler()
{

}

void StatisticsSampler::sample(System *system, double dt)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    m_numAtoms = system->atomCount();
    m_dt = dt;
    // ...
}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    system->for_each([&](AtomBlock &block, int i) {
        float vx = block.velocity->x[i];
        float vy = block.velocity->y[i];
        float vz = block.velocity->z[i];
        float lensq = vx*vx + vy*vy + vz*vz;
        m_kineticEnergy += 0.5 * system->atom()->mass() * lensq;
    });
}
