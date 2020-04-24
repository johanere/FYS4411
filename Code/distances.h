#pragma once
#include <vector>
#include "../particle.h"
class Distances {
public:
    Distances(class System* system,int numberOfDimensions, int numberOfParticles);
    void InitiateR(std::vector<Particle*> particles);
    void calculateR_jk(std::vector<Particle*> particles, int particle);
    double getR_jk(int j, int k);
private:
    int     m_numberOfDimensions = 0;
    int     m_numberOfParticles = 0;
    std::vector<double> m_particleDistances = std::vector<double>();
    class System* m_system = nullptr;
};
