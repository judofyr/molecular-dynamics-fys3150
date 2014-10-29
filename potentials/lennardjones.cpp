#include <potentials/lennardjones.h>
#include <cmath>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0;
    auto &atoms = system->atoms();

    for (size_t i = 0; i < atoms.size(); i++) {
        for (size_t j = i+1; j < atoms.size(); j++) {
            auto atom1 = atoms[i];
            auto atom2 = atoms[j];
            vec3 r = system->distanceVectorBetweenAtoms(atom1, atom2);
            double dr = r.length();

            double U = 4*m_epsilon*(pow(m_sigma/dr, 12) - pow(m_sigma/dr, 6));
            double F = -4*m_epsilon*(-12*pow(m_sigma/dr, 12) + 6*pow(m_sigma/dr, 6))/dr;
            m_potentialEnergy += U;

            r.normalize();
            atom1->force.addAndMultiply(r, F);
            atom2->force.addAndMultiply(r, -F);
        }
    }
}
