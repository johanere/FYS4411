#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include "system.h"
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "RBM.h"
#include "sampler.h"
#include "random.h"
#include "writefile.h"


using namespace std;


int main(int argc, char* argv[]) {


double omega = 1.0;
double sigma = 1.0;


int numberOfParticles = 2;
int numberOfDimensions = 3;
int m= numberOfParticles*numberOfDimensions ;

int n=2;

int interaction = 0;

int GD_iters = 10;

double equilibration = 0.0;
double stepLength = 0.5;
int numberOfSteps = (int) pow(2,15);
int method=1;

srand((unsigned int) time(0));

double learningrate=0.5;


RBM* rbm = new RBM(GD_iters,m,n,learningrate,sigma,omega);


  for (int run=1; run<=GD_iters;run++)
    {
    System* system = new System();

    system->setRBM                      (rbm,run);

    //system->setDistances                (interaction,a);
    Hamiltonian* ham = new Hamiltonian(system);
    Wavefunction* wf = new Wavefunction(system);
    system->setHamiltonian              (ham);

    system->setWavefunction             (wf);

    system->setEquilibrationFraction    (equilibration);

    system->setStepLength               (stepLength);

    system->runMetropolisSteps          (numberOfSteps,method);

    system->getWavefunction()->evaluate(system->get_X());

    system->getSampler()->printOutputToTerminal();



    rbm->Update_gradients(system->getSampler()->Return_gradients_a(),
    system->getSampler()->Return_gradients_b(), system->getSampler()->
    Return_gradients_W());

    if (run==GD_iters){
    system->getSampler()->blocking();
    }
    delete ham;
    delete wf;
    delete system;
    } // end of run loop

    return 0;
}



// remember to change seed in randomuniform
