#include "berendsenthermostat.h"
#include <cmath>

BerendsenThermostat::BerendsenThermostat(double T_bath, double relaxation) :
    m_T_bath(T_bath), m_relaxation(relaxation)
{
}

void BerendsenThermostat::adjustTemperature(System *system, StatisticsSampler *sampler)
{
    double scale = sampler->dt() / m_relaxation;
    double T = sampler->instantaneousTemperature();
    double gamma = sqrt(1 + scale*(m_T_bath/T - 1));

    system->for_each([&](AtomBlock &block, int i){
        block.velocity->x[i] *= gamma;
        block.velocity->y[i] *= gamma;
        block.velocity->z[i] *= gamma;
    });
}
