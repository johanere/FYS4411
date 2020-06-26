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

// system settings
int n;
int P = 2;
int D = 2;
int m= P*D ;
double omega = 1.0;
double sigma = 1.0;
int interaction = 1;

//VMC settings
int method=0;
double stepLength = 1.0;
int numberOfSteps = (int) pow(2,21);
double equilibration = 0.5;

//Optimization settings
int GD_iters = 200;
double learningrate=0.3;
double learningrate_initial=0.3;

//output
string folder_name="../../Results/sim_16/";
string tag;

if (argc < 1){
  cout << "Bad Usage: " << argv[0] <<" Too few CL arguments" << endl;
  exit(1);}
if (argc > 1){
  method = atoi(argv[1]);
  n = atoi(argv[2]);
}
tag=to_string(method);
tag.append("_");
tag.append(to_string(n));
tag.append("_");

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


int average_runs=5; // Used to average over noise
int checkrunstotal=1; // Used to run variations of parameters

for (int a_run=1; a_run<=average_runs;a_run++)
  {
    Eigen::VectorXd avr_error         = Eigen::VectorXd::Zero(GD_iters);
    Eigen::VectorXd avr_en_deviation  = Eigen::VectorXd::Zero(GD_iters);


  for (int checkrun=0; checkrun<checkrunstotal;checkrun++)
    {
      // learningrate +=0.2
      RBM* rbm = new RBM(GD_iters,m,n,learningrate,sigma,omega);
      for (int run=1; run<=GD_iters;run++)
        {
        cout<<"Averaging run"<<a_run<<" run:"<<run<<endl;
        System* system = new System();

        system->setRBM                      (rbm,run);

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

        avr_en_deviation(run-1)= (system->getSampler()->getEnergy()-3.0)*(system->getSampler()->getEnergy()-3.0);
        avr_error(run-1)=system->getSampler()->getSE();

        delete ham;
        delete wf;
        delete system;
        } // end of run loop

        delete rbm;
    } //end of checkrun

  write_vector(avr_en_deviation,tag+to_string(a_run)+"_deltaE",folder_name);
  write_vector(avr_error,tag+to_string(a_run)+"_error",folder_name);

  }
  return 0;
}
