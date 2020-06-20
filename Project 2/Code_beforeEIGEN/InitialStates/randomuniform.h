#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles,double a, int interaction );
    void setupInitialState(double a, int interaction);
    std::vector<double> GeneratePosition();
};
