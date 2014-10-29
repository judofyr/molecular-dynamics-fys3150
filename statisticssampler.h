#pragma once
#include <system.h>

class StatisticsSampler
{
private:
    double m_kineticEnergy;
    int m_numAtoms;
    double m_dt;

public:
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system, double dt);
    void sampleKineticEnergy(System *system);
    double kineticEnergy();
    double instantaneousTemperature();
    double dt();
};
