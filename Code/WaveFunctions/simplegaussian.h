#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, int numberOfDimensions, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    std::vector<double> updateQForce(std::vector<class Particle*> particles,int particle);
};
