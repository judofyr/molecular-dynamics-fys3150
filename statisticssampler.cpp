#include <statisticssampler.h>


double StatisticsSampler::kineticEnergy()
{
    return m_kineticEnergy;
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
    // ...
}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    for (auto &atom : system->atoms()) {
        m_kineticEnergy += 0.5 * atom->mass() * atom->velocity.lengthSquared();
    }
}
