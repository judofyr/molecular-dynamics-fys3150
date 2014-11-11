#include <system.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <unitconverter.h>
#include <cassert>
#include <algorithm>
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
    m_velocities.push_back(vel);
    block.velocity = &m_velocities.back();
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
        for (int i = 0; i < block.count; i++) {
            // TODO: Don't convert to vec3
            vec3 vel = block.velocity->vec3(i);
            netMomentum.addAndMultiply(vel, m_atom->mass());
            atomCount++;
        }
    }

    // We want to remove this much momentum from each atom
    vec3 perMomentum = netMomentum / atomCount;

    for (AtomBlock &block : m_atomBlocks) {
        for (int i = 0; i < block.count; i++) {
            block.velocity->x[i] += perMomentum.x() / m_atom->mass();
            block.velocity->y[i] += perMomentum.y() / m_atom->mass();
            block.velocity->z[i] += perMomentum.z() / m_atom->mass();
        }
    }
}

void System::resetForcesOnAllAtoms() {
    for (AtomBlock &block : m_atomBlocks) {
        for (int i = 0; i < block.count; i++) {
            block.force.x[i] = 0;
            block.force.y[i] = 0;
            block.force.z[i] = 0;
        }
    }
}

void System::addAtom(double x, double y, double z) {
    AtomBlock *block = &m_atomBlocks.back();
    if (block->count == ATOMBLOCKSIZE) {
        // This block is full. Create a new block.
        createBlock();
        block = &m_atomBlocks.back();
    }

    // Current position in block
    int i = block->count++;
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

    int requiredBlocks = atomCount / ATOMBLOCKSIZE + 1;
    m_velocities.reserve(requiredBlocks);
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
    ghostBlock.type = AtomBlockType::GHOST;
    AtomRef ref;
    AtomRef *refPtr = NULL;

    m_ghostedAtoms.clear();
    m_ghostBlocks.clear();

    for (int dim = 0; dim < 3; dim++) {
        auto handleAtom = [&](AtomBlock &block, int i){
            float *posArray = block.position.x + ATOMBLOCKSIZE*dim;
            float pos = posArray[i];
            bool isGhost = false;

            if (pos < m_rShell) {
                pos += m_systemSize[dim];
                isGhost = true;
            } else if (pos >= m_systemSize[dim] - m_rShell) {
                pos -= m_systemSize[dim];
                isGhost = true;
            }

            if (isGhost) {
                // Store reference to the original atom
                if (refPtr) {
                    ref = *refPtr;
                } else {
                    ref.block = &block;
                    ref.idx = i;
                }
                m_ghostedAtoms.push_back(ref);

                if (ghostBlock.count == ATOMBLOCKSIZE) {
                    m_ghostBlocks.push_back(ghostBlock);
                    // Reset our version
                    ghostBlock.count = 0;
                }

                int gi = ghostBlock.count++;


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
        int refI = 0;

        // Move over the other ghost blocks
        for (size_t blockIdx = 0; blockIdx < ghostCount; blockIdx++) {
            for (int i = 0; i < ATOMBLOCKSIZE; i++) {
                refPtr = &m_ghostedAtoms[refI];
                handleAtom(m_ghostBlocks[blockIdx], i);
                refI++;
            }
        }

        // Move the current ghost block we're working on
        for (int i = 0; i < currentGhostBlock.count; i++) {
            refPtr = &m_ghostedAtoms[refI];
            handleAtom(currentGhostBlock, i);
            refI++;
        }

        // Move over all regular atoms
        refPtr = NULL;
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

        assert(block.position.x[i] >= 0);
        assert(block.position.y[i] >= 0);
        assert(block.position.z[i] >= 0);
        assert(block.position.x[i] < sizeX);
        assert(block.position.y[i] < sizeY);
        assert(block.position.z[i] < sizeZ);

    });

    CPElapsedTimer::getInstance().periodicBoundaryConditions().stop();
}

size_t System::atomCount()
{
    return m_atomBlocks.size() * ATOMBLOCKSIZE;
}

void System::calculateForces() {

    if (m_steps % 10 == 0) {
        applyPeriodicBoundaryConditions();
        applyPeriodicGhostBlocks();
        buildCellLists();
        buildNeighbourLists();
    }

    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::buildCellLists()
{
    CPElapsedTimer::getInstance().updateCellList().start();

    m_cellLists.clear();

    // We add 3 extra cells in each dimension: one for each side, and one for round-off.
    int cellSizeX = m_systemSize.x()/m_rShell + 3;
    int cellSizeY = m_systemSize.y()/m_rShell + 3;
    int cellSizeZ = m_systemSize.z()/m_rShell + 3;

    m_cellLists.resize(cellSizeX*cellSizeY*cellSizeZ);

    auto handleAtom = [&](AtomBlock &block, int i) {
        int cellX = block.position.x[i]/m_rShell+1;
        int cellY = block.position.y[i]/m_rShell+1;
        int cellZ = block.position.z[i]/m_rShell+1;
        int index = cellX + cellY*cellSizeX + cellZ*cellSizeX*cellSizeY;
        auto &blocks = m_cellLists[index];
        AtomRef ref;
        ref.block = &block;
        ref.idx = i;
        blocks.push_back(ref);
    };

    for_each(handleAtom);

    for (auto &block : m_ghostBlocks) {
        for (int i = 0; i < block.count; i++) {
            handleAtom(block, i);
        }
    }

    CPElapsedTimer::getInstance().updateCellList().stop();
}

void System::buildNeighbourLists()
{
    CPElapsedTimer::updateNeighborList().start();

    int cellSizeX = m_systemSize.x()/m_rShell + 3;
    int cellSizeY = m_systemSize.y()/m_rShell + 3;
    int cellSizeZ = m_systemSize.z()/m_rShell + 3;

    float cutOff2 = m_rShell*m_rShell;

    for (auto &block : m_atomBlocks) {
        block.clearNeighbours();
    }

    for (auto &block : m_ghostBlocks) {
        block.clearNeighbours();
    }

    for (int cellX = 1; cellX < cellSizeX - 1; cellX++)
    for (int cellY = 1; cellY < cellSizeY - 1; cellY++)
    for (int cellZ = 1; cellZ < cellSizeZ - 1; cellZ++) {
        int cellIndex = cellX + cellY*cellSizeX + cellZ*cellSizeX*cellSizeY;
        auto &refs = m_cellLists[cellIndex];

        auto ref1 = &refs[0];
        auto ref1end = ref1 + refs.size();

        for (; ref1 != ref1end; ref1++) {
            for (int cdx = -1; cdx <= 1; cdx++)
            for (int cdy = -1; cdy <= 1; cdy++)
            for (int cdz = -1; cdz <= 1; cdz++) {
                AtomRef *ref2, *ref2end;

                int otherCellIndex = (cellX+cdx) + (cellY+cdy)*cellSizeX + (cellZ+cdz)*cellSizeX*cellSizeY;
                auto &otherRefs = m_cellLists[otherCellIndex];
                ref2 = &otherRefs[0];
                ref2end = ref2 + otherRefs.size();

                for (; ref2 != ref2end; ref2++) {
                    if (ref1->block >= ref2->block)
                        continue;

                    float x1 = ref1->block->position.x[ref1->idx];
                    float y1 = ref1->block->position.y[ref1->idx];
                    float z1 = ref1->block->position.z[ref1->idx];

                    float x2 = ref2->block->position.x[ref2->idx];
                    float y2 = ref2->block->position.y[ref2->idx];
                    float z2 = ref2->block->position.z[ref2->idx];

                    float dx = x1 - x2;
                    float dy = y1 - y2;
                    float dz = z1 - z2;

                    float dr2 = dx*dx + dy*dy + dz*dz;
                    if (dr2 < cutOff2) {
                        ref1->block->addNeighbour(ref2->block);
                    }
                }
            }
        }
    }

    CPElapsedTimer::updateNeighborList().stop();
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}

void System::for_each(std::function<void (AtomBlock &, int)> action)
{
    for (auto &block : m_atomBlocks) {
        for (int i = 0; i < block.count; i++) {
            action(block, i);
        }
    }
}

void System::for_each_ghost(std::function<void (AtomBlock &, int, AtomRef &)> action)
{
    auto origRef = m_ghostedAtoms.begin();
    for (auto &block : m_ghostBlocks) {
        for (int i = 0; i < block.count; i++) {
            action(block, i, *origRef);
            origRef++;
        }
    }
}
