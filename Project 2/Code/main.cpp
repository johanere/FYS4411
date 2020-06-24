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

int n=2;
int P = 1;
int D = 1;
int m= P*D ;
double omega = 1.0;
double sigma = 1.0;
int interaction = 0;

int method=0;
double stepLength = 1.0;
int numberOfSteps = (int) pow(2,21);
double equilibration = 0.5;
int GD_iters = 15;
double learningrate=0.3;

string folder_name="../../Results/sim_5/";

if (argc <= 2){
  cout << "Bad Usage: " << argv[0] <<" Too few CL arguments" << endl;
  exit(1);}
if (argc > 2){
  method = atoi(argv[1]);
  n=       atoi(argv[2]);
}

folder_name.append(to_string(method));
folder_name.append("/");

if (method==1) stepLength = 2.0;

srand((unsigned int) time(0));


int checkrunstotal=10;
learningrate=0.24;


std::vector <double> errors_last;
errors_last.reserve(checkrunstotal);

for (int checkrun=1; checkrun<=checkrunstotal;checkrun++)
  {
    learningrate+=0.02;
    cout<<learningrate<<endl;

    string filename= to_string(learningrate);

    std::vector <double> energy_gd;
    std::vector <double> errors;
    energy_gd.reserve(GD_iters);

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

    system->setParameters( equilibration, stepLength, D, P, interaction);

    system->runMetropolisSteps          (numberOfSteps,method);

    system->getWavefunction()->evaluate(system->get_X());

  //  system->getSampler()->printOutputToTerminal();

    rbm->Update_gradients(system->getSampler()->Return_gradients_a(),
    system->getSampler()->Return_gradients_b(), system->getSampler()->
    Return_gradients_W());

    energy_gd.push_back(system->getSampler()->getEnergy());

    if (run==GD_iters){
    system->getSampler()->blocking();
    errors_last.push_back(system->getSampler()->getSE());
    }

    delete ham;
    delete wf;
    delete system;
    } // end of run loop
    write_LocalEnergy(energy_gd,  filename , folder_name) ;

  }
  write_LocalEnergy(errors_last,  "errors" , folder_name) ;
    return 0;
}



// remember to change seed in randomuniform
