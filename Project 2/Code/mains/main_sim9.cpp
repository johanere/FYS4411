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

string folder_name="../../Results/sim_9/";

string tag;
if (argc < 1){
  cout << "Bad Usage: " << argv[0] <<" Too few CL arguments" << endl;
  exit(1);}
if (argc > 1){
  method = atoi(argv[1]);
  tag = argv[2];
}
n=2;

cout<<"Method "<<method<< " eta "<< learningrate<<endl;
cout<<"To folder "<<folder_name<<endl;
//srand((unsigned int) time(0));

learningrate=0.1;

int checkrunstotal=1;
int average_runs=9;


sigma=0.914;
for (int a_run=1; a_run<=average_runs;a_run++)
  {
    cout<<"run "<< a_run<<endl;


  for (int checkrun=0; checkrun<checkrunstotal;checkrun++)
    {
      learningrate+=0.1;
      cout<<" "<<"lr "<<learningrate<<endl;

      Eigen::VectorXd energy_vec       = Eigen::VectorXd::Zero(GD_iters);

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

        rbm->Update_gradients(system->getSampler()->Return_gradients_a(),
        system->getSampler()->Return_gradients_b(), system->getSampler()->
        Return_gradients_W());


        energy_vec(run-1)= system->getSampler()->getEnergy();

        delete ham;
        delete wf;
        delete system;
        } // end of run loop
        write_vector(energy_vec,to_string(learningrate),folder_name);
        delete rbm;
    } //end of checkrun

} // end of a_run


  return 0;
}



// remember to change seed in randomuniform
