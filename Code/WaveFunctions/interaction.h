#pragma once
#include "wavefunction.h"

class interactingWF : public WaveFunction {
public:
    interactingWF(class System* system, int numberOfDimensions, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    std::vector<double> updateQForce(std::vector<class Particle*> particles,int particle);
};
