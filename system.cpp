#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "cpelapsedtimer.h"
#include <stdlib.h>
#include <cmath>
#include <cassert>

using namespace std;

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0),
    m_initialized(false),
    m_shouldSample(false)
{
    m_miniAtoms = new MiniAtoms();
    m_atoms = new Atoms();
}

System::~System()
{
    delete m_potential;
    delete m_integrator;
}

void System::initialize(MDDataType_t cutoffRadius) {
    m_rCut = cutoffRadius;
    m_rShell = m_rCut * 1.12;
    m_neighborList.setup(this, m_rShell);
    printStatus();
    validate();
    m_neighborList.update();
    m_initialized = true;
}

void System::validate() {
#ifdef BENCHMARK
    return;
#endif
    CellList &cellList = m_neighborList.cellList();
    if(cellList.numberOfCellsX() < 3 || cellList.numberOfCellsY() < 3 || cellList.numberOfCellsZ() < 3) {
        cout << "Error, system size too small to have at least 3 cells in each dimension, aborting." << endl;
        cout << "Minimum system size: " << UnitConverter::lengthToAngstroms(3*m_rShell) << " Å" << endl;
        exit(1);
    }
}

void System::printStatus() {
    cout << "System initialized with " << m_atoms->numberOfAtoms << " atoms." << endl;
    cout << "System size: " << m_systemSize << endl;
    CellList &cellList = m_neighborList.cellList();
    cout << "Number of cells: " << cellList.numberOfCellsX() << ", " << cellList.numberOfCellsY() << ", " << cellList.numberOfCellsZ() << endl;
    cout << "Cell size: " << UnitConverter::lengthToAngstroms(cellList.lengthX()) << " x " << UnitConverter::lengthToAngstroms(cellList.lengthY()) << " x " << UnitConverter::lengthToAngstroms(cellList.lengthZ()) << " Å^3" << endl;
}

void System::applyPeriodicBoundaryConditions() {
    CPElapsedTimer::periodicBoundaryConditions().start();
#ifdef MD_SIMD
#pragma simd
#endif
    for(unsigned int i=0; i<m_atoms->numberOfAtoms; i++) {
#ifdef SINGLEPRECISION
        // Fixes the bug where the position is slightly negative, but due to float precision, both tests (<0) and >=size) pass
        if(m_atoms->x[i] < 0 && m_atoms->x[i] > -1e-6) m_atoms->x[i] = 0;
        if(m_atoms->y[i] < 0 && m_atoms->y[i] > -1e-6) m_atoms->y[i] = 0;
        if(m_atoms->z[i] < 0 && m_atoms->z[i] > -1e-6) m_atoms->z[i] = 0;
#endif
        if(m_atoms->x[i] < 0) {
            m_atoms->x[i] += m_systemSize[0];
            m_neighborList.xInitial[i] += m_systemSize[0];
        } else if(m_atoms->x[i] >= m_systemSize[0]) {
            m_atoms->x[i] -= m_systemSize[0];
            m_neighborList.xInitial[i] -= m_systemSize[0];
        }
        if(m_atoms->y[i] < 0) {
            m_atoms->y[i] += m_systemSize[1];
            m_neighborList.yInitial[i] += m_systemSize[1];
        } else if(m_atoms->y[i] >= m_systemSize[1]) {
            m_atoms->y[i] -= m_systemSize[1];
            m_neighborList.yInitial[i] -= m_systemSize[1];
        }
        if(m_atoms->z[i] < 0) {
            m_atoms->z[i] += m_systemSize[2];
            m_neighborList.zInitial[i] += m_systemSize[2];
        } else if(m_atoms->z[i] >= m_systemSize[2]) {
            m_atoms->z[i] -= m_systemSize[2];
            m_neighborList.zInitial[i] -= m_systemSize[2];
        }
    }
    CPElapsedTimer::periodicBoundaryConditions().stop();
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
    StatisticsSampler sampler(0);
    vec3 momentum = sampler.sampleMomentum(this);
    
    momentum /= m_atoms->numberOfAtoms;
    for(int i=0; i<m_atoms->numberOfAtoms; i++) {
        m_atoms->vx[i] -= momentum[0]*m_atoms->inverseMass[i];
        m_atoms->vy[i] -= momentum[1]*m_atoms->inverseMass[i];
        m_atoms->vz[i] -= momentum[2]*m_atoms->inverseMass[i];
    }
}

void System::resetForcesOnAllAtoms() {
#ifdef MD_SIMD
#pragma simd
#endif
    for(unsigned int i=0; i<m_atoms->numberOfAtoms; i++) {
        m_atoms->fx[i] = 0;
        m_atoms->fy[i] = 0;
        m_atoms->fz[i] = 0;
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, float latticeConstant, float temperature) {
    MDDataType_t xCell[4] = {0, 0.5, 0.5, 0};
    MDDataType_t yCell[4] = {0, 0.5, 0, 0.5};
    MDDataType_t zCell[4] = {0, 0, 0.5, 0.5};
    int numAtoms = 4*numberOfUnitCellsEachDimension*numberOfUnitCellsEachDimension*numberOfUnitCellsEachDimension;
    assert(MAXNUMATOMS >= numAtoms && "Too many atoms, increase MAXNUMATOMS in atoms.h");

    for(int i=0; i< numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j< numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k< numberOfUnitCellsEachDimension; k++) {
                for(int l=0; l<4; l++) {
                    int atomIndex = m_atoms->numberOfAtoms;
                    m_atoms->inverseMass[atomIndex] = 1.0/UnitConverter::massFromSI(6.63352088e-26);
                    MDDataType_t x = (i+xCell[l])*latticeConstant;
                    MDDataType_t y = (j+yCell[l])*latticeConstant;
                    MDDataType_t z = (k+zCell[l])*latticeConstant;
                    m_atoms->x[atomIndex] = x;
                    m_atoms->y[atomIndex] = y;
                    m_atoms->z[atomIndex] = z;
                    
                    vec3 velocity;
                    MDDataType_t boltzmannConstant = 1.0; // In atomic units, the boltzmann constant equals 1
                    MDDataType_t standardDeviation = sqrt(boltzmannConstant*temperature*m_atoms->inverseMass[atomIndex]);
                    velocity.randomGaussian(0, standardDeviation);
                    m_atoms->vx[atomIndex] = velocity[0];
                    m_atoms->vy[atomIndex] = velocity[1];
                    m_atoms->vz[atomIndex] = velocity[2];
                    m_atoms->index[atomIndex] = atomIndex;
                    
                    m_atoms->numberOfAtoms++;
                }
            }
        }
    }
    
    float sideLength = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(sideLength, sideLength, sideLength));
    cout << "Added " << m_atoms->numberOfAtoms << " atoms in an FCC lattice." << endl;
}

void System::setShouldSample(bool shouldSample)
{
    m_shouldSample = shouldSample;
}

void System::createGhostAtoms()
{
    //    m_atoms->numberOfGhostAtoms = 0; // Reset all ghosts
    
    //    for(unsigned short dimension=0; dimension<3; dimension++) {
    //        unsigned int newGhosts = m_atoms->numberOfGhostAtoms;
    
    //        for(unsigned int i=0; i<m_atoms->numberOfAtoms+newGhosts; i++) {
    //            bool shouldCopy = false;
    //            unsigned int atomIndex = m_atoms->numberOfAtomsIncludingGhosts();
    
    //            if(dimension==0) {
    //                if(m_atoms->x[i] < m_rShell) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i] + m_systemSize[0]; shouldCopy = true;
    //                } else if (m_atoms->x[i] > m_systemSize[0] - m_rShell) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i] - m_systemSize[0]; shouldCopy = true;
    //                }
    
    //                if(shouldCopy) {
    //                    m_atoms->y[atomIndex] = m_atoms->y[i];
    //                    m_atoms->z[atomIndex] = m_atoms->z[i];
    //                }
    //            } else if(dimension==1) {
    //                if(m_atoms->y[i] < m_rShell) {
    //                    m_atoms->y[atomIndex] = m_atoms->y[i] + m_systemSize[1]; shouldCopy = true;
    //                } else if (m_atoms->y[i] > m_systemSize[1] - m_rShell) {
    //                    m_atoms->y[atomIndex] = m_atoms->y[i] - m_systemSize[1]; shouldCopy = true;
    //                }
    
    //                if(shouldCopy) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i];
    //                    m_atoms->z[atomIndex] = m_atoms->z[i];
    //                }
    //            } else if(dimension==2) {
    //                if(m_atoms->z[i] < m_rShell) {
    //                    m_atoms->z[atomIndex] = m_atoms->z[i] + m_systemSize[2]; shouldCopy = true;
    //                } else if (m_atoms->z[i] > m_systemSize[2] - m_rShell) {
    //                    m_atoms->z[atomIndex] = m_atoms->z[i] - m_systemSize[2]; shouldCopy = true;
    //                }
    
    //                if(shouldCopy) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i];
    //                    m_atoms->y[atomIndex] = m_atoms->y[i];
    //                }
    //            }
    
    //            if(shouldCopy) {
    //                m_atoms->vx[atomIndex] = m_atoms->vx[i];
    //                m_atoms->vy[atomIndex] = m_atoms->vy[i];
    //                m_atoms->vz[atomIndex] = m_atoms->vz[i];
    
    //                m_atoms->mass[atomIndex] = m_atoms->mass[i];
    //                m_atoms->index[atomIndex] = atomIndex;
    //                m_atoms->numberOfGhostAtoms++;
    //            }
    //        }
    //    }
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    if(m_shouldSample) m_potential->calculateForcesAndEnergyAndPressure(this);
    else m_potential->calculateForces(this);
}

void System::step(double dt) {
    if(!m_initialized) {
        cout << "System not initialized, aborting. Remember to call initialize()." << endl;
        exit(1);
    }
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
