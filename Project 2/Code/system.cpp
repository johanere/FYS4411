#include <cassert>
#include <math.h>
#include <iostream>

#include "system.h"
#include "sampler.h"
#include "Wavefunction.h"
#include "random.h"
#include "RBM.h"

#include "headers.h"

#include <math.h>
#include <iostream>
using std::cout;
using std::endl;


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

System::~System(void) {
   delete m_sampler;
}

bool System::brute_force_Step() {
    assert(m_stepLength>0);
    //pick random particle
    move_index=m_random->nextInt(0,m_M-1);
    // make copy of old pos
    oldState=m_X;
    //move particle
    m_X[move_index] += 2 * (m_random->nextDouble()-0.5) * m_stepLength;
    //calculate new distances for moved particle
    if (m_interaction>0){
    ;//  calculateR_jk(m_particles, particle); //CHECK THIS LATER
    }
    //calculate new WF
    newWF = m_Wavefunction->evaluate(m_X);

  /*  if (m_interaction>0){ //CHECK THIS LATER
        if (m_break_eval>0){
        m_break_eval=0;
        m_particles[particle]->setPosition(oldPosition);
        calculateR_jk(m_particles, particle);
        return false;}
    } */
    // check if move in accepted, if not accepted, reset particle
    if(m_random->nextDouble() > (newWF*newWF) / (oldWF*oldWF)) {
      m_X=oldState;
      //reset distances
      if (m_interaction>0){
      ;//  calculateR_jk(m_particles, particle);
      }
      return false;
    }
    else  {
      oldWF = newWF;
      return true;
    }
} //end of brute force sampling


bool System::importance_sampling_Step() {
    assert(m_stepLength>0);
    //pick random particle
    move_index=m_random->nextInt(0,m_M-1);
    // make copy of old pos
    oldState=m_X;
    double oldQForce = m_Wavefunction->getQForce(m_X,move_index);
    m_X[move_index] +=   0.5*oldQForce*m_stepLength+m_random->nextGaussian(0,1)*m_stepLengthsqrt;
/*
    if (m_interaction>0){
      calculateR_jk(m_particles, particle);
    }*/
    //calculate new WF
    newWF = m_Wavefunction->evaluate(m_X);

    //auto reject new proposal if any particles come with a distance less than a
    /*if (m_interaction>0){
        if (m_break_eval>0){
        m_break_eval=0;
        m_particles[particle]->setPosition(oldPosition);
        calculateR_jk(m_particles, particle);
        return false;}
    }*/


    double newQForce = m_Wavefunction->getQForce(m_X,move_index);

    //Greens
    double gfRatio=0;
    double term1=0;
    double term2=0;


    term1 = oldState[move_index]- m_X[move_index] - 0.5 *m_stepLength*newQForce;
    term2 = m_X[move_index] - oldState[move_index] - 0.5 *m_stepLength*oldQForce;
    gfRatio +=  - (term1*term1) + (term2*term2);
    gfRatio= gfRatio/(2.0*m_stepLength);
    gfRatio = exp(gfRatio);


    if(m_random->nextDouble() > gfRatio*(newWF*newWF)  / (oldWF*oldWF)) {
      //if not accepted, reset particle
      m_X=oldState;
      return false;
      if (m_interaction>0){
        //calculateR_jk(m_particles, particle);
        ;
      }
    }
    else  {
      oldWF = newWF;
      return true;
    }

} //end of improtance sampling step


// functions needed for gibbs sampling
double System::logistic(double x){
  return (1.0/ (1.0+exp(x) ) ) ;
}

double System::v_j(int j,Eigen::VectorXd X){
  double sum = 0;
    for (int i=0; i<getRBM()->get_M();i++)
    {
      sum+= X(i)* getRBM()->get_W()(i,j) ;
    }
  //sum =  X.dot(m_W.col(j)) ;
  return getRBM()->get_b()(j)+sum / (getRBM()->get_sigma() * getRBM()->get_sigma()) ;
}

double System::h_sampling(int j){
  return logistic(-v_j(j,m_X));
}

void System::X_sampling(){
  double mean;
  double std = m_rbm->get_sigma() * m_rbm->get_sigma();
  for (int i=0; i<m_M ;i++)
  {
    mean =0;
    for (int j=0; j<m_N ;j++)
    {
      mean+=  getRBM()->get_W()(i,j)*m_H[j];
    }
    mean +=getRBM()->get_a()(i);
    if (mean>1.0){
    //cout<<mean<<endl;
  ;}
    m_X[i]=m_random->nextGaussian(mean, std);
  }
}

bool System::gibbs_Step() {
  for (int j=0; j<m_N ;j++)
  {
    double y=h_sampling(j);
    if(y > m_random->nextDouble())
    {      m_H[j]=1.0; }
    else { m_H[j]=0.0; }
  }
  X_sampling();
  return true;
} //end of gibbs sampling

void System::runMetropolisSteps(int numberOfMetropolisSteps, int method) {

    m_X                         = Eigen::VectorXd::Random(m_M);
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    if(method  == 2.0) m_gibbsfactor=0.5;
    else  m_gibbsfactor=1;



    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps,
      numberOfMetropolisSteps-m_equilibrationFraction*numberOfMetropolisSteps,true);

    //InitiateR(); //calculate distances

    oldWF = m_Wavefunction->evaluate(m_X);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
      //check method chhoice

        if (method == 0) acceptedStep = brute_force_Step();
        else if (method == 1) acceptedStep = importance_sampling_Step();
        else if (method == 2) acceptedStep = gibbsStep11();
        else { cout<<"Provide valid method selection"<<endl; abort(); }

        // sample

        if (i>=m_equilibrationFraction*numberOfMetropolisSteps)
        {
          m_sampler->sample(acceptedStep);

        }
      }
    m_sampler->computeAverages();
}




void System::setParameters(double equilibrationFraction,double stepLength, double D, double P, int interaction) {
    assert(equilibrationFraction >= 0);
    assert(stepLength >= 0);
    m_stepLength = stepLength;
    m_stepLengthsqrt=sqrt(m_stepLength);
    m_equilibrationFraction = equilibrationFraction;
    m_D=D;
    m_P=P;
    m_interaction=interaction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWavefunction(Wavefunction* Wavefunction) {
    m_Wavefunction = Wavefunction;
}


void System::setRBM(RBM* rbm, int current_run) {
    m_rbm = rbm;
    m_M=m_rbm -> get_M();
    m_N=m_rbm -> get_N();
    m_rbm -> WeightsAndBiases(current_run);
    m_X= 2*(Eigen::VectorXd::Random(m_M)-0.5*Eigen::VectorXd::Ones(m_M));
    m_H= Eigen::VectorXd::Zero(m_N);
}

bool System::gibbsStep11(){
    /*
    Adjust positions of particles according to the Gibbs sampling rule
    */

    // Calculate conditional probabilities of the hidden nodes
    for(int j = 0; j < getRBM()->get_N(); j++){
         double z = getRBM()->get_b()(j) + m_X.dot(getRBM()->get_W().col(j))/(m_rbm->get_sigma()*m_rbm->get_sigma() );
         m_H(j) = m_random->nextDouble() < sigmoid(z);
    }

    // Update particle positions
    for(int i = 0; i < getRBM()->get_M(); i++){

        double meanPos = getRBM()->get_a()(i) + getRBM()->get_W().row(i)*m_H;

        double posDistribution = m_random->nextGaussian(meanPos, m_rbm->get_sigma());
        m_X(i) = posDistribution;
    }
      return true;
    }

double System::sigmoid(double z){
    return 1.0/(1 + exp(-z));
}
