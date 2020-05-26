#pragma once
#include "wavefunction.h"

class interactingWF : public WaveFunction {
public:
    interactingWF(class System* system,  int numberOfDimensions, double alpha, double beta, double a);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    std::vector<double> updateQForce(std::vector<class Particle*> particles,int particle);

  private:
      class Distances* m_distances = nullptr;
      double m_a =0;
};
