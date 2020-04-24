#pragma once
#include "hamiltonian.h"
#include <vector>

class interactingHO : public Hamiltonian {
public:
    interactingHO(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};
