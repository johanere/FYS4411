#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <cmath>

#include <chrono>
#include "sampler.h"
#include "writefile.h"
#include <vector>

using namespace std;


int main() {
    bool write=false;
    int seed = 2020;
    int method =1; // 0 - brute force metropolis, 1 - metropolis hastings


    int numberOfDimensions  = 3;
    int numberOfParticles   = 5;
    int numberOfSteps       = (int) 1e6;  //(int) 2*1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;//0.5;  //omega/2.0; //0.5;          // Variational parameter.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.3;          // Amount of the total steps used
    // for equilibration.
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, numberOfDimensions,alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    system->runMetropolisSteps          (numberOfSteps,method);

    system->getWaveFunction()->evaluate(system->getParticles());

    chrono::steady_clock::time_point end = chrono::steady_clock::now();


    cout<<"Closed form "<<omega*numberOfDimensions*numberOfParticles/2<<endl;
    cout<<"Time used "<<chrono::duration_cast<chrono::milliseconds>(end - begin).count()<< " ms"<<endl;
    std::vector<double> energySamples = system->getSampler()->getEnergySamples();

    if (write==true){
      write_LocalEnergy(energySamples,"../../Results/Raw/output_test");
    }


    return 0;
}


// remember to change seed in randomuniform
