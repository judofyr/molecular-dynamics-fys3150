#pragma once
#include <system.h>

class StatisticsSampler
{
private:
    double m_kineticEnergy;
    int m_numAtoms;

public:
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
    double kineticEnergy();
    double instantaneousTemperature();

};
