#pragma once
#include "hamiltonian.h"
#include <vector>

class H_RBM_Simple : public Hamiltonian {
public:
    interactingHO(System* system, double omega,double a);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
    double m_a=0;
};
