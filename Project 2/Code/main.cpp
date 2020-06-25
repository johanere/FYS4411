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

int n=4;
int P = 2;
int D = 2;
int m= P*D ;
double omega = 1.0;
double sigma = 1.0;
int interaction = 1;

int method=0;
double stepLength = 1.0;
int numberOfSteps = (int) pow(2,21);
double equilibration = 0.5;
int GD_iters = 20;
double learningrate=0.3;
double learningrate_initial=0.3;

string folder_name="../../Results/sim_12/";

string tag;
if (argc < 1){
  cout << "Bad Usage: " << argv[0] <<" Too few CL arguments" << endl;
  exit(1);}
if (argc > 1){
  method = atoi(argv[1]);
}
tag=to_string(method);

cout<<"Method "<<method<< " n "<< n<<endl;
cout<<"To folder "<<folder_name<<endl;
//srand((unsigned int) time(0));
if (method == 0) {
  learningrate = 0.38;
  stepLength = 1.0;
 }
else if (method == 1) {
  learningrate = 0.34;
  stepLength = 2.0;
}
else if (method == 2) {
  learningrate = 0.26;
  sigma=0.914;
  }
else { cout<<"Provide valid method selection"<<endl; abort(); }

cout<<" "<<"eta "<<learningrate<<" sigma "<<sigma<<endl;


int checkrunstotal=1;
int average_runs=10;

Eigen::VectorXd avr_error         = Eigen::VectorXd::Zero(GD_iters);
Eigen::VectorXd avr_en_deviation  = Eigen::VectorXd::Zero(GD_iters);
Eigen::VectorXd energy_vector     = Eigen::VectorXd::Zero(GD_iters);

for (int a_run=1; a_run<=average_runs;a_run++)
  {
    Eigen::VectorXd error_for_fitting = Eigen::VectorXd::Zero(GD_iters);
  for (int checkrun=0; checkrun<checkrunstotal;checkrun++)
    {
      RBM* rbm = new RBM(GD_iters,m,n,learningrate,sigma,omega);
      for (int run=1; run<=GD_iters;run++)
        {
        cout<<"Average run"<<a_run<<" run:"<<run<<endl;
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


        system->getSampler()->blocking();
        error_for_fitting(run-1)=system->getSampler()->getEnergy();
        avr_en_deviation(run-1)+= (system->getSampler()->getEnergy()-3.0)*(system->getSampler()->getEnergy()-3.0);
        avr_error(run-1)+=system->getSampler()->getSE();
        cout<<"e="<<system->getSampler()->getEnergy()<<endl;
        cout<<"de2="<<(system->getSampler()->getEnergy()-3.0)*(system->getSampler()->getEnergy()-3.0)<<endl;
        energy_vector(run-1)+=system->getSampler()->getEnergy();

        if (a_run==average_runs){
          if (run==GD_iters){
          cout<<"printing"<<endl;
          write_LocalEnergy(system->getSampler()->getEnergySamples(),tag+"_energy_one_iter",folder_name);}
        }
        if (run==GD_iters){
            cout<<"printing for fitting"<<endl;
        write_vector(error_for_fitting,tag+to_string(a_run)+"e_for_reg",folder_name+tag+"/");
        }

        delete ham;
        delete wf;
        delete system;
        } // end of run loop

        delete rbm;
    } //end of checkrun

} // end of a_run

  avr_en_deviation=avr_en_deviation/average_runs;
  avr_error=avr_error/average_runs;
  energy_vector=energy_vector/average_runs;

  write_vector(avr_en_deviation,tag+"_deltaE",folder_name);
  write_vector(avr_error,tag+"_error",folder_name);
  write_vector(energy_vector,tag+"_energy",folder_name);
  return 0;
}



// remember to change seed in randomuniform
