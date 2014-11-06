#include <system.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <unitconverter.h>
#include <cassert>
#include "cpelapsedtimer.h"

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0)
{
}


void System::createBlock()
{
    AtomVec vel;
    AtomBlock block;
    m_velcoties.push_back(vel);
    block.velocity = &m_velcoties.back();
    m_atomBlocks.push_back(block);
}

System::~System()
{
    delete m_potential;
    delete m_integrator;
}

void System::removeMomentum() {
    // TODO: Why do we even first multiply by mass, and then divide by mass?
    vec3 netMomentum;
    int atomCount = 0;
    for (AtomBlock &block : m_atomBlocks) {
        for (int i = 0; i < ATOMBLOCKSIZE; i++) {
            // TODO: Don't convert to vec3
            vec3 vel = block.velocity->vec3(i);
            netMomentum.addAndMultiply(vel, m_atom->mass());
            atomCount++;
        }
    }

    // We want to remove this much momentum from each atom
    vec3 perMomentum = netMomentum / atomCount;

    for (AtomBlock &block : m_atomBlocks) {
        for (int i = 0; i < ATOMBLOCKSIZE; i++) {
            block.velocity->x[i] += perMomentum.x() / m_atom->mass();
            block.velocity->y[i] += perMomentum.y() / m_atom->mass();
            block.velocity->z[i] += perMomentum.z() / m_atom->mass();
        }
    }
}

void System::resetForcesOnAllAtoms() {
    for (AtomBlock &block : m_atomBlocks) {
        for (int i = 0; i < ATOMBLOCKSIZE; i++) {
            block.force.x[i] = 0;
            block.force.y[i] = 0;
            block.force.z[i] = 0;
        }
    }
}

void System::addAtom(double x, double y, double z) {
    AtomBlock *block = &m_atomBlocks.back();
    if (block->counter == ATOMBLOCKSIZE) {
        // This block is full. Create a new block.
        createBlock();
        block = &m_atomBlocks.back();
    }

    // Current position in block
    int i = block->counter++;
    block->position.x[i] = x;
    block->position.y[i] = y;
    block->position.z[i] = z;

    vec3 vel = m_atom->velocityMaxwellian(UnitConverter::temperatureFromSI(300));

    block->velocity->x[i] = vel.x();
    block->velocity->y[i] = vel.y();
    block->velocity->z[i] = vel.z();

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant) {
    int atomCount = numberOfUnitCellsEachDimension*numberOfUnitCellsEachDimension*numberOfUnitCellsEachDimension*4;
    assert(atomCount % ATOMBLOCKSIZE == 0 && "Atoms should be multiple of blocksize");

    int requiredBlocks = atomCount / ATOMBLOCKSIZE + 1;
    m_velcoties.reserve(requiredBlocks);
    m_atomBlocks.reserve(requiredBlocks);

    // Create the initial block
    createBlock();

    double b = latticeConstant;
    for (int x_i = 0; x_i < numberOfUnitCellsEachDimension; x_i++) {
    for (int y_i = 0; y_i < numberOfUnitCellsEachDimension; y_i++) {
    for (int z_i = 0; z_i < numberOfUnitCellsEachDimension; z_i++) {
        addAtom(b*x_i + 0,     b*y_i + 0,     b*z_i + 0);
        addAtom(b*x_i + (b/2), b*y_i + (b/2), b*z_i + 0);
        addAtom(b*x_i + 0,     b*y_i + (b/2), b*z_i + (b/2));
        addAtom(b*x_i + (b/2), b*y_i + 0,     b*z_i + (b/2));

    }}}
}

void System::applyPeriodicGhostBlocks()
{
    CPElapsedTimer::getInstance().createGhosts().start();
    AtomBlock ghostBlock;
    ghostBlock.counter = 0;

    m_ghostBlocks.clear();

    for (int dim = 0; dim < 3; dim++) {
        auto handleAtom = [&](AtomBlock &block, int i){
            float *posArray = block.position.x + ATOMBLOCKSIZE*dim;
            float pos = posArray[i];
            bool isGhost = false;

            if (pos < m_rCutOff) {
                pos += m_systemSize[dim];
                isGhost = true;
            } else if (pos >= m_systemSize[dim] - m_rCutOff) {
                pos -= m_systemSize[dim];
                isGhost = true;
            }

            if (isGhost) {
                if (ghostBlock.counter == ATOMBLOCKSIZE) {
                    m_ghostBlocks.push_back(ghostBlock);
                    // Reset our version
                    ghostBlock.counter = 0;
                }

                int gi = ghostBlock.counter++;

                // Copy the atom position
                ghostBlock.position.x[gi] = block.position.x[i];
                ghostBlock.position.y[gi] = block.position.y[i];
                ghostBlock.position.z[gi] = block.position.z[i];

                // Override the current dimension
                float *ghostPosArray = ghostBlock.position.x + ATOMBLOCKSIZE*dim;
                ghostPosArray[gi] = pos;
            }
        };

        // Figure out how far into the ghost blocks we are
        size_t ghostCount = m_ghostBlocks.size();
        AtomBlock currentGhostBlock = ghostBlock;

        // Move the current ghost block we're working on
        for (int i = 0; i < currentGhostBlock.counter; i++) {
            handleAtom(currentGhostBlock, i);
        }

        // Move over the other ghost blocks
        for (size_t blockIdx = 0; blockIdx < ghostCount; blockIdx++) {
            for (int i = 0; i < ATOMBLOCKSIZE; i++) {
                handleAtom(m_ghostBlocks[blockIdx], i);
            }
        }

        // Move over all regular atoms
        for_each(handleAtom);
    }

    m_ghostBlocks.push_back(ghostBlock);
    CPElapsedTimer::getInstance().createGhosts().stop();
}

void System::applyPeriodicBoundaryConditions()
{
    CPElapsedTimer::getInstance().periodicBoundaryConditions().start();
    float sizeX = m_systemSize.x();
    float sizeY = m_systemSize.y();
    float sizeZ = m_systemSize.z();

    for_each([&](AtomBlock &block, int i) {
        // Apply periodic boundary conditions
        if (block.position.x[i] < 0) {
            block.position.x[i] += sizeX;
        } else if (block.position.x[i] >= sizeX) {
            block.position.x[i] -= sizeX;
        }

        if (block.position.y[i] < 0) {
            block.position.y[i] += sizeY;
        } else if (block.position.y[i] >= sizeY) {
            block.position.y[i] -= sizeY;
        }

        if (block.position.z[i] < 0) {
            block.position.z[i] += sizeZ;
        } else if (block.position.z[i] >= sizeZ) {
            block.position.z[i] -= sizeZ;
        }
    });

    CPElapsedTimer::getInstance().periodicBoundaryConditions().stop();
}

size_t System::atomCount()
{
    return m_atomBlocks.size() * ATOMBLOCKSIZE;
}

void System::calculateForces() {
    applyPeriodicBoundaryConditions();
    applyPeriodicGhostBlocks();
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}

void System::for_each(std::function<void (AtomBlock &, int)> action)
{
    for (auto &block : m_atomBlocks) {
        for (int i = 0; i < ATOMBLOCKSIZE; i++) {
            action(block, i);
        }
    }
}
