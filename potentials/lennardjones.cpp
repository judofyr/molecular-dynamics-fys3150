#include <potentials/lennardjones.h>
#include <cmath>
#include <cassert>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0;
    float rCut2 = system->rCutOff() * system->rCutOff();

    auto calculateForce = [&](AtomBlock &block, int i, AtomBlock &otherBlock, int j, bool isGhost) {
        float x1 = block.position.x[i];
        float y1 = block.position.y[i];
        float z1 = block.position.z[i];

        float x2 = otherBlock.position.x[j];
        float y2 = otherBlock.position.y[j];
        float z2 = otherBlock.position.z[j];

        float dx = x1 - x2;
        float dy = y1 - y2;
        float dz = z1 - z2;

        float dr2 = dx*dx + dy*dy + dz*dz;

        assert(dr2 != 0 && "dr should never be exactly zero");
        if (dr2 > rCut2) return;

        float oneOverDr2 = 1.0/dr2;
        float oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

        float sigmaSixth = pow(m_sigma, 6.0);
        float force = -24*m_epsilon*sigmaSixth*oneOverDr6*(2*sigmaSixth*oneOverDr6 - 1)*oneOverDr2;

        float oneOverDrCut2 = 1.0/rCut2;
        float oneOverDrCut6 = oneOverDrCut2*oneOverDrCut2*oneOverDrCut2;
        float potentialEnergyCutoff = 4*m_epsilon*sigmaSixth*oneOverDrCut6*(sigmaSixth*oneOverDrCut6 - 1);
        float potentialEnergy = 4*m_epsilon*sigmaSixth*oneOverDr6*(sigmaSixth*oneOverDr6 - 1) - potentialEnergyCutoff;
        if (isGhost) potentialEnergy *= 0.5;
        m_potentialEnergy += potentialEnergy;

        block.force.x[i] += dx * -force;
        block.force.y[i] += dy * -force;
        block.force.z[i] += dz * -force;
        otherBlock.force.x[j] += dx * force;
        otherBlock.force.y[j] += dy * force;
        otherBlock.force.z[j] += dz * force;
    };

    auto &blocks = system->m_atomBlocks;
    auto endBlock = blocks.end();
    auto &ghostBlocks = system->m_ghostBlocks;
    auto endGhostBlock = ghostBlocks.end();

    // Go over each atom
    for (auto block = blocks.begin(); block != endBlock; block++) {
        for (int i = 0; i < ATOMBLOCKSIZE; i++) {
            // Compare against the other atoms
            auto otherBlock = block;
            int j = i+1;

            while (true) {
                if (j == ATOMBLOCKSIZE) {
                    otherBlock++;
                    if (otherBlock == endBlock)
                        break;
                    j = 0;
                }
                calculateForce(*block, i, *otherBlock, j, false);
                j++;
            }

            // Compare against all ghost atoms
            for (otherBlock = ghostBlocks.begin(); otherBlock != endGhostBlock; otherBlock++) {
                for (int j = 0; j < otherBlock->counter; j++) {
                    calculateForce(*block, i, *otherBlock, j, true);
                }
            }
        }
    }
}
