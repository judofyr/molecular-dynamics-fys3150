#pragma once
#include "thermostat.h"

class BerendsenThermostat : public Thermostat
{
private:
    double m_T_bath;
    double m_relaxation;

public:
    BerendsenThermostat(double T_bath, double relaxation);
    void adjustTemperature(System *system, StatisticsSampler *sampler);
};
