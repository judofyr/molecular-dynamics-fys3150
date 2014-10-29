#include <statisticssampler.h>


double StatisticsSampler::kineticEnergy()
{
    return m_kineticEnergy;
}

double StatisticsSampler::instantaneousTemperature()
{
    return (2.0/3) * m_kineticEnergy / m_numAtoms;
}

StatisticsSampler::StatisticsSampler()
{

}

StatisticsSampler::~StatisticsSampler()
{

}

void StatisticsSampler::sample(System *system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    m_numAtoms = system->atoms().size();
    // ...
}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    for (auto &atom : system->atoms()) {
        m_kineticEnergy += 0.5 * atom->mass() * atom->velocity.lengthSquared();
    }
}
