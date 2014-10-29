#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/velocityverlet.h>
#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>

using namespace std;

int main()
{
    double dt = UnitConverter::timeFromSI(1e-15); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;

    System system;
    double b = 5.26;
    int numCells = 5;
    system.setSystemSize(UnitConverter::lengthFromAngstroms(vec3(b*numCells, b*numCells, b*numCells)));
    system.createFCCLattice(numCells, UnitConverter::lengthFromAngstroms(b));
    system.setPotential(new LennardJones(3.405, 1.0)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.removeMomentum();

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt);
        statisticsSampler->sample(&system);
        double energy = system.potential()->potentialEnergy() + statisticsSampler->kineticEnergy();
        cout << "step=" << timestep << " energy=" << energy << endl;


        movie->saveState(&system);
    }

    movie->close();

    return 0;
}
