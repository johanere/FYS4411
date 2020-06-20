#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interaction.h"

#include "WaveFunctions/RBM_SimpleWF.h"

#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/interactingHO.h"

#include "Hamiltonians/H_RBM_Simple.h"

#include "RBM.h"

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

double a = 0.0000000001; // hard shell exclusion

double omega = 1.0;
double sigma = 1.0;
int numberOfParticles = 2;
int numberOfDimensions = 2;
int interaction = 0;
int n=2;
int GD_iters = 1;
int runs = 1;
int m= numberOfParticles*numberOfDimensions ;
double equilibration = 0.3;
double stepLength = 0.5;
int numberOfSteps = (int) pow(2,3);
int method=0;

int seed = 1;
double alpha = 0.5;
double beta = 1;
double learningrate=0.01;


RBM* rbm = new RBM(GD_iters,m,n,learningrate,sigma,seed);


  for (int run=1; run<=runs;run++)
    {
    System* system = new System(seed);

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, a,  interaction ));
    system->setDistances                (interaction,a);


    system->setRBM                      (rbm,run);

    // Pick the right wave function
    if (interaction==0){
    system->setHamiltonian              (new H_RBM_Simple(system,omega,sigma));
    system->setWaveFunction             (new RBM_SimpleWF(system, numberOfDimensions,omega,sigma) );

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
    system->getSampler()->printOutputToTerminal();
    //system->getSampler()->Return_gradients();

    } // end of run loop
    return 0;
}



// remember to change seed in randomuniform
