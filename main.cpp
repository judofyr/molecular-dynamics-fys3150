#include <config.h>
#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/velocityverlet.h>
#include <thermostats/berendsenthermostat.h>
#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>
#include "cpelapsedtimer.h"

using namespace std;

int main()
{
    int numTimeSteps = 1000;
    double dt = UnitConverter::timeFromSI(1e-14); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;

    System system;
    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); // Argon mass
    double b = 5.26;
    double sigma = 3.405;
    int numCells = 8;
    system.setSystemSize(UnitConverter::lengthFromAngstroms(vec3(b*numCells, b*numCells, b*numCells)));
    system.setRCutOff(2.5*sigma); // 2.5 * sigma is the case where the force is tiny
    system.setRShell(2.8*sigma);
    system.setAtom(atom);
    system.createFCCLattice(numCells, UnitConverter::lengthFromAngstroms(b));
    system.setPotential(new LennardJones(sigma, 1.0)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.removeMomentum();
    system.applyPeriodicGhostBlocks();

    auto thermostat = new BerendsenThermostat(UnitConverter::temperatureFromSI(1000), dt*10);

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    for(int timestep=0; timestep<numTimeSteps; timestep++) {
        system.step(dt);
        // movie->saveState(&system, false);

        if(!(timestep % 10)) {
#if MD_SAMPLE
            statisticsSampler->sample(&system, dt);
            double energy = system.potential()->potentialEnergy() + statisticsSampler->kineticEnergy();
            double temp = statisticsSampler->instantaneousTemperature();
            cout << "step=" << timestep << " energy=" << energy << " temp=" << UnitConverter::temperatureToSI(temp) << endl;
#else
            cout << "step=" << timestep << endl;
#endif

#if MD_NEIGHBOURS
            int neighbourCount = 0;
            for (auto &block : system.m_atomBlocks) {
                neighbourCount += block.neighbourCount();
            }
            cout << "neighbours=" << neighbourCount << endl;
#endif
        }

        if (timestep > 200 && timestep < 400) {
            //thermostat->adjustTemperature(&system, statisticsSampler);
        }

    }

    movie->close();


    float calculateForcesFraction = CPElapsedTimer::calculateForces().elapsedTime() / CPElapsedTimer::totalTime();
    float halfKickFraction = CPElapsedTimer::halfKick().elapsedTime() / CPElapsedTimer::totalTime();
    float moveFraction = CPElapsedTimer::move().elapsedTime() / CPElapsedTimer::totalTime();
    float createGhostsFraction = CPElapsedTimer::createGhosts().elapsedTime() / CPElapsedTimer::totalTime();

    float updateNeighborListFraction = CPElapsedTimer::updateNeighborList().elapsedTime() / CPElapsedTimer::totalTime();
    float updateCellListFraction = CPElapsedTimer::updateCellList().elapsedTime() / CPElapsedTimer::totalTime();
    float periodicBoundaryConditionsFraction = CPElapsedTimer::periodicBoundaryConditions().elapsedTime() / CPElapsedTimer::totalTime();
    // float samplingFraction = CPElapsedTimer::sampling().elapsedTime() / CPElapsedTimer::totalTime();
    // float timeEvolutionFraction = CPElapsedTimer::timeEvolution().elapsedTime() / CPElapsedTimer::totalTime();

    cout << endl << "Program finished after " << CPElapsedTimer::totalTime() << " seconds. Time analysis:" << endl;
    cout << fixed
         //<< "      Time evolution    : " << CPElapsedTimer::timeEvolution().elapsedTime() << " s ( " << 100*timeEvolutionFraction << "%)" <<  endl
         << "      Force calculation : " << CPElapsedTimer::calculateForces().elapsedTime() << " s ( " << 100*calculateForcesFraction << "%)" <<  endl
         << "      Moving            : " << CPElapsedTimer::move().elapsedTime() << " s ( " << 100*moveFraction << "%)" <<  endl
         << "      Half kick         : " << CPElapsedTimer::halfKick().elapsedTime() << " s ( " << 100*halfKickFraction << "%)" <<  endl
         << "      Create ghosts     : " << CPElapsedTimer::createGhosts().elapsedTime() << " s ( " << 100*createGhostsFraction << "%)" <<  endl
         << "      Update neighbors  : " << CPElapsedTimer::updateNeighborList().elapsedTime() << " s ( " << 100*updateNeighborListFraction << "%)" <<  endl
         << "      Update cells      : " << CPElapsedTimer::updateCellList().elapsedTime() << " s ( " << 100*updateCellListFraction << "%)" <<  endl
         << "      Periodic boundary : " << CPElapsedTimer::periodicBoundaryConditions().elapsedTime() << " s ( " << 100*periodicBoundaryConditionsFraction << "%)" <<  endl;
         // << "      Sampling          : " << CPElapsedTimer::sampling().elapsedTime() << " s ( " << 100*samplingFraction << "%)" <<  endl;
    cout << endl << numTimeSteps / CPElapsedTimer::totalTime() << " timesteps / second. " << endl;
    cout << system.atomCount()*numTimeSteps / (1000*CPElapsedTimer::totalTime()) << "k atom-timesteps / second. " << endl;
    cout << (system.atomCount() + system.m_ghostBlocks.size()*MD_BLOCKSIZE) *numTimeSteps / (1000*CPElapsedTimer::totalTime()) << "k atom-timesteps / second (including ghosts). " << endl;

    return 0;
}
