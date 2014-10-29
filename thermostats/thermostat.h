#pragma once
#include <statisticssampler.h>
#include <system.h>

class Thermostat
{
public:
    Thermostat();
    virtual ~Thermostat() {}
    virtual void adjustTemperature(System *system, StatisticsSampler *sampler) = 0;
};

