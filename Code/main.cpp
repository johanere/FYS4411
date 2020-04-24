#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"

#include "WaveFunctions/interaction.h"

#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"

#include "Hamiltonians/interactingHO.h"

#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <cmath>

#include "distances.h"

#include <chrono>
#include "sampler.h"
#include "writefile.h"
#include <vector>

using namespace std;


int main() {
    bool write=false;
    int seed = 12341;
    int method =0; // 0 - brute force metropolis, 1 - metropolis hastings
    int wf=0; // 0 - non-interacting, 1- interacting

    int numberOfDimensions  = 3;
    int numberOfParticles   = 4;
    int numberOfSteps       = (int) pow(2,6);    //(int) 2*1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.3;          // Variational parameter.
    double stepLength       = 0.3;          // Metropolis step length.
    double equilibration    = 0.3;          // Amount of the total steps used
    // for equilibration.
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    System* system = new System(seed);


    if (wf==0){
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, numberOfDimensions,alpha));
    }
    else if (wf==1){
    system->setHamiltonian              (new interactingHO(system, omega));
    system->setWaveFunction             (new interactingWF(system, numberOfDimensions,alpha));
    }
    else {
      cout<<"Provide valid wf selection"<<endl;
      abort();
    }

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setDistances                (new Distances(system, numberOfDimensions,numberOfParticles));

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    system->runMetropolisSteps          (numberOfSteps,method);
    system->getWaveFunction()->evaluate(system->getParticles());






    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout<<"Closed form "<<omega*numberOfDimensions*numberOfParticles/2<<endl;
    cout<<"Time used "<<chrono::duration_cast<chrono::milliseconds>(end - begin).count()<< " ms"<<endl;
    std::vector<double> energySamples = system->getSampler()->getEnergySamples();

    if (write==true){
      write_LocalEnergy(energySamples,"../../Results/Raw/Energies");
    }


    return 0;
}


// remember to change seed in randomuniform
