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
#include "sampler.h"
#include "writefile.h"
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

int main(int argc, char* argv[]) {



  for (int run=1; run<=runs;run++)
    {
    System* system = new System();

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, a,  interaction ));
    system->setDistances                (interaction,a);

    // Pick the right wave function
    if (interaction==0){
    system->setHamiltonian              (new HarmonicOscillator(system, 1));
    system->setWaveFunction             (new SimpleGaussian(system, numberOfDimensions,alpha,beta));
    }
    else if (interaction==1){
    system->setHamiltonian              (new interactingHO(system, 1,a));
    system->setWaveFunction             (new interactingWF(system, numberOfDimensions,alpha,beta,a));
    }

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    // Run metropolis
    system->runMetropolisSteps          (numberOfSteps,method,GD_iters);

    system->getWaveFunction()->evaluate(system->getParticles());

    // update weights here

    } // end of run loop

  return 0;
}


// remember to change seed in randomuniform
