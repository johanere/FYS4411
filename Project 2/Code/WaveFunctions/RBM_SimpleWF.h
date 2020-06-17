#pragma once
#include "wavefunction.h"

class RBM_SimpleWF : public WaveFunction {
public:
    RBM_SimpleWF(class System* system,  int numberOfDimensions, std::vector<double> w, std::vector<double> a, std::vector<double> b,double sigma));
    double evaluate(std::vector<class Particle*> particles);

  private:
      class Distances*                m_distances = nullptr;
      std::vector<double>             m_w=std::vector<double>();
      std::vector<double>             m_a=std::vector<double>();
      std::vector<double>             m_b=std::vector<double>();
      double                          m_sigma = 1;
};
