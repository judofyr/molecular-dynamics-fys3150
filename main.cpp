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

using namespace std;

int main()
{
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
    int numCells = 4;
    system.setSystemSize(UnitConverter::lengthFromAngstroms(vec3(b*numCells, b*numCells, b*numCells)));
    system.setRCutOff(2.5*sigma); // 2.5 * sigma is the case where the force is tiny
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

    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt);
        movie->saveState(&system, false);

        statisticsSampler->sample(&system, dt);
        double energy = system.potential()->potentialEnergy() + statisticsSampler->kineticEnergy();
        double temp = statisticsSampler->instantaneousTemperature();
        cout << "step=" << timestep << " energy=" << energy << " temp=" << UnitConverter::temperatureToSI(temp) << endl;

        if (timestep > 200 && timestep < 400) {
            //thermostat->adjustTemperature(&system, statisticsSampler);
        }

    }

    movie->close();

    return 0;
}
