#pragma once
#include "wavefunction.h"
#include "../RBM.h"
#include "../Hamiltonians/H_RBM_Simple.h"

class RBM_SimpleWF : public WaveFunction {
public:
    RBM_SimpleWF(System* system,  int numberOfDimensions,double omega,double sigma);
    double evaluate(std::vector<class Particle*> particles);
    double compute_v_j(int j) ;
    std::vector<double> updateQForce(std::vector<class Particle*> particles,int particle);

  private:
    //  class Distances*                m_distances = nullptr;
      std::vector<double>             m_W=std::vector<double>();
      std::vector<double>             m_a=std::vector<double>();
      std::vector<double>             m_b=std::vector<double>();
      std::vector<double>             m_X=std::vector<double>();

      double m_M=0;
      double m_N=0;

      double m_sigma=0;
      double m_omega=0;
};
