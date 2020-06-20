#pragma once
#include "hamiltonian.h"
#include "../RBM.h"
#include <vector>

class H_RBM_Simple : public Hamiltonian {
public:
    H_RBM_Simple(System* system,double omega,double sigma);
    double computeLocalEnergy(std::vector<Particle*> particles);
    double compute_v_j(int j) ;
    double compute_nabla_ln_psi(int k) ;
    double compute_nabla2_ln_psi(int k);
    

private:
    std::vector<double>             m_W=std::vector<double>();
    std::vector<double>             m_a=std::vector<double>();
    std::vector<double>             m_b=std::vector<double>();
    std::vector<double>             m_X;
    double m_M=0;
    double m_N=0;

    double m_sigma=0;
    double m_omega=0;
};
