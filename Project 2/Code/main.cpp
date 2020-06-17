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
  double beta=1.0;
  double delta_alpha=0;
  double stepLength=0.5,alpha;
  int numberOfParticles,method,interaction,eliptical,GD_iters;

  // read from command line
  if (argc <= 4){
    cout << "Bad Usage: " << argv[0] <<" Too few CL arguments" << endl;
    exit(1);}
  if (argc > 4){
    numberOfParticles = atoi(argv[1]);
    interaction=atoi(argv[2]);
    eliptical=atoi(argv[3]);
    alpha = atof(argv[4]);}

  // manage output to files and terminal
  string folder_name = "../../Results/Densitites/";
  bool    write_E=false;    //print energies to file
  bool    write_R=true;   //print distances to file
  bool    timing=1;        //set to 1 to print time to terminal, else set to 0
  bool    print_singlerun=1; // set to 1 to print each run to terminal, else set to 0
  bool    write_last=0;

  //system
  double  a                = 0.0043;
  int     numberOfDimensions  = 3;

  //method
          method = 1;                                 // 0 - non-interacting, 1- interacting
  int     numberOfSteps       = (int) pow(2,18);
  double  equilibration    = 0.3;                     // Amount of the total steps used for equilibration.
          stepLength = 0.01;                         // Metropolis step length.

  //manage multiple simulations in one execution
  double  gamma=0.5;        //initial learning rate
          GD_iters=0;  //number of iterations of gradient decent
  int     runs=125; //number of runs
  double  run_delta_alpha=-0.00; // delta alpha on each run when GD is not used


  beta = 1.0;
  if (method==0){ stepLength = 0.5;} //adjust step length if using importance sampling 0.3 / 0.01 / 0.1 teste
  if (eliptical==1){beta = 2.82843;}   // adjust beta if using eliptical potential
  if (GD_iters>0){runs =GD_iters;}  //set runs to be gradient decent runs


  //make folder for results
  string main_folder=folder_name;
  folder_name.append("_");
  folder_name.append(to_string(numberOfParticles));
  folder_name.append("_");
  folder_name.append(to_string(interaction));

  std::string sPath = main_folder;
  mode_t nMode = 0733; // UNIX style permissions
  int nError = 0;
  #if defined(_WIN32)
    nError = _mkdir(sPath.c_str()); // can be used on Windows
  #else
    nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
  #endif
  if (nError != 0) {
    cout<<"Error in creating new directory. Directory may already exist, if not check nError when writing file."<<endl;
  }

  double total_accepted_steps=0;
  double total_sigma=0;
  double total_time=0;

  double E_min=100000;
  double alpha_E_min=0;
  double run_E_min=0;

    for (int run=1; run<=runs;run++)
    {
    alpha=alpha - gamma* delta_alpha;
    cout<<"\n"<<endl;
    cout<<"Running VMC run "<<run<<" for N="<<numberOfParticles<<", d="<<numberOfDimensions<<", alpha= "<<alpha<<" and beta = "<<beta<<endl;
    cout<<"Using method "<<method<<", dt="<<stepLength<<" with interaction = "<<interaction<<" cycles" <<endl;

    string filename = folder_name; filename.append("/");
    filename.append(to_string(run));
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    System* system = new System();


    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, a,  interaction ));
    system->setDistances                (interaction,a);

    if (interaction==0){
    system->setHamiltonian              (new HarmonicOscillator(system, 1));
    system->setWaveFunction             (new SimpleGaussian(system, numberOfDimensions,alpha,beta));
    }
    else if (interaction==1){
    system->setHamiltonian              (new interactingHO(system, 1,a));
    system->setWaveFunction             (new interactingWF(system, numberOfDimensions,alpha,beta,a));
    }
    else {
      cout<<"Provide valid wf selection"<<endl;
      abort();
    }



    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    system->runMetropolisSteps          (numberOfSteps,method,GD_iters);
    if (print_singlerun>0){      system->getSampler()->printOutputToTerminal();}
    system->getWaveFunction()->evaluate(system->getParticles());






    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    if (timing==true){
    cout<<"Time used "<<chrono::duration_cast<chrono::milliseconds>(end - begin).count()<< " ms"<<endl;
    total_time+=chrono::duration_cast<chrono::milliseconds>(end - begin).count();}

    std::vector<double> energySamples = system->getSampler()->getEnergySamples();

    if (runs>0){
       total_sigma+=system->getSampler()->getError();
       total_accepted_steps+=system->getSampler()->getAcceptanceRate();
      }


    if (write_E==true){

      if (run==1){
        std::string sPath = folder_name;
        mode_t nMode = 0733; // UNIX style permissions
        int nError = 0;
        #if defined(_WIN32)
          nError = _mkdir(sPath.c_str()); // can be used on Windows
        #else
          nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
        #endif
        if (nError != 0) {
          cout<<"Error in creating new directory. Directory may already exist, if not check nError when writing file."<<endl;
        }
      }
      if (write_last==true){
         if (run==runs){write_LocalEnergy(energySamples,filename);}
         else {;}
        }
      else {
      write_LocalEnergy(energySamples,filename); }

    }
    if (write_R==true){
              if (run==1){
                std::string sPath = folder_name;
                mode_t nMode = 0733; // UNIX style permissions
                int nError = 0;
                #if defined(_WIN32)
                  nError = _mkdir(sPath.c_str()); // can be used on Windows
                #else
                  nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
                #endif
                if (nError != 0) {
                  cout<<"Error in creating new directory. Directory may already exist, if not check nError when writing file."<<endl;
                }
              }


      std::vector<double> distances = system->getRadialDistances();
      write_distances(distances,filename);
    }
    if (GD_iters>0){delta_alpha=system->getSampler()->getGradAlpha();
    gamma=gamma*0.9;
    cout<<"Run " <<run<<" New Gamma="<<gamma<<endl;
    if (E_min>system->getSampler()->getEnergy()){
      E_min=system->getSampler()->getEnergy();
      alpha_E_min=alpha;
      run_E_min=run;
    }
    }

    else {delta_alpha=run_delta_alpha;}



  } // end of run loop
  if (runs>0){

  cout<<"------ Averages -------"<<endl;
  cout<<"Accepted steps        = "<<total_accepted_steps/((double) runs)<<endl;
  cout << std::scientific;
  cout<<"Sigma (not corrected) = "<<total_sigma/((double) runs)<<endl;
  cout << std::fixed;
  if (timing==true){
  cout<<"Time                  = "<<total_time/((double) runs)<<endl; }
  write_average_of_runs(runs,total_accepted_steps,  total_sigma,  total_time,folder_name);
  }
  if (GD_iters>0){cout<<"Min E="<<E_min<<" at alpha="<<alpha_E_min<<" run nr="<<run_E_min<<endl;}
  return 0;
}


// remember to change seed in randomuniform
