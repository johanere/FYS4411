#include "wavefunction.h"
#include <iostream>
using namespace std;

WaveFunction::WaveFunction(System* system, int numberOfDimensions) {
    m_system = system;
    m_qForce = std::vector<double> (numberOfDimensions,0);
}
