#include "berendsenthermostat.h"
#include "statisticssampler.h"
#include "system.h"
#include <cmath>

BerendsenThermostat::BerendsenThermostat(double temperature, double relaxationFactor) :
    m_temperature(temperature),
    m_relaxationFactor(relaxationFactor)
{
}

void BerendsenThermostat::apply(System *system, StatisticsSampler *sampler)
{
    MDDataType_t scaleFactor = sqrt(1 + m_relaxationFactor*(sampler->temperature()/m_temperature - 1));
    Atoms &atoms = system->atoms();
    for(int i=0; i<atoms.numberOfAtoms; i++) {
        atoms.vx[i] *= scaleFactor;
        atoms.vy[i] *= scaleFactor;
        atoms.vz[i] *= scaleFactor;
    }
}
