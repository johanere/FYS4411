#include <iostream>
#include <cmath>
#include <vector>

#include "sampler.h"
#include "system.h"
#include "Wavefunction.h"
#include "RBM.h"
#include "Hamiltonian.h"

#include "blocker.h"
#include <fstream>
#include <string>
using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
    m_acceptedSteps = 0;

    m_psi_a = Eigen::VectorXd::Zero(m_system->getRBM()->get_M());
    m_psi_a_EL = Eigen::VectorXd::Zero(m_system->getRBM()->get_M());

    m_psi_b =   Eigen::VectorXd::Zero(m_system->getRBM()->get_N());
    m_psi_b_EL = Eigen::VectorXd::Zero(m_system->getRBM()->get_N());

    m_psi_W=    Eigen::MatrixXd::Zero(m_system->getRBM()->get_M(),m_system->getRBM()->get_N() );
    m_psi_W_EL= Eigen::MatrixXd::Zero(m_system->getRBM()->get_M(),m_system->getRBM()->get_N());

} //end of   Sampler

Sampler::~Sampler(void) {
}

void Sampler::setNumberOfMetropolisSteps(int steps, int steps_after_eq, bool store_samples) {
    m_numberOfMetropolisSteps = steps;
    if (store_samples==true){m_EnergySamples.reserve(steps_after_eq);}
    else {cout<<"not storing energy samples"<<endl;}

} //end of  setNumberOfMetropolisSteps


void Sampler::sample(bool acceptedStep) {

    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergy2 = 0;}
    m_DeltaE=m_system->getHamiltonian()->computeLocalEnergy(m_system->get_X() );

    m_cumulativeEnergy  += m_DeltaE;
    m_cumulativeEnergy2 += m_DeltaE*m_DeltaE;
    m_stepNumber++;
    m_EnergySamples.push_back(m_DeltaE);

    Update_expectations(m_DeltaE);

    if (acceptedStep==true){m_acceptedSteps++;}
    else {;}

} //end of sample

void Sampler::printOutputToTerminal() {
    cout << "E:     " << m_energy << endl;
  //  cout << std::scientific;
  //  cout << "S:     " << m_error << endl;
  //  cout << std::fixed;
  //  cout << "A:     " <<  ((double) m_acceptedSteps)/m_stepNumber << endl;
    cout << " \n" <<endl;
} //end of printOutputToTerminal

void Sampler::computeAverages() {
    m_energy    =   m_cumulativeEnergy/m_stepNumber;
    m_energy2   =   m_cumulativeEnergy2/m_stepNumber;
    m_variance  =   m_energy2 - (m_energy*m_energy);
    m_error     =   sqrt(m_variance/m_stepNumber);
    m_A         =   ((double) m_acceptedSteps)  / ((double) m_stepNumber);
} //end of computeAverages

void Sampler::blocking() {
   Blocker block(m_EnergySamples);
   m_mse_mean=block.mse_mean;
   m_stdErr=block.stdErr;
   m_mse_stdErr=block.mse_stdErr;
   printf("Expected value = %g (with mean sq. err. = %g) \n", block.mean, block.mse_mean);
   printf("Standard error = %g (with mean sq. err. = %g) \n", block.stdErr, block.mse_stdErr);
}

void Sampler::Update_expectations(double m_DeltaE){

  Eigen::VectorXd X= m_system->get_X();
  Eigen::VectorXd a= m_system->getRBM()->get_a();
  Eigen::VectorXd b= m_system->getRBM()->get_b();
  Eigen::MatrixXd W= m_system->getRBM()->get_W();
  double sigma =  m_system->getRBM()->get_sigma();

  m_psi_a    +=  m_system->get_gibbsfactor()* (X - a)/(sigma*sigma);
  m_psi_a_EL +=  m_system->get_gibbsfactor()* (X - a)/(sigma*sigma)   * m_DeltaE;

  for (int j=0; j<m_system->getRBM()->get_N();j++)
  {

    m_psi_b(j)    += m_system->get_gibbsfactor()* logistic(-v_j(j,X));
    m_psi_b_EL(j) += m_system->get_gibbsfactor()* logistic(-v_j(j,X)) * m_DeltaE;
  }

  for (int i=0; i<m_system->getRBM()->get_M();i++)
  {
    for (int j=0; j<m_system->getRBM()->get_N();j++)
    {
      m_psi_W(i,j)    += m_system->get_gibbsfactor()* X(i) * logistic(-v_j(j,X))  /(sigma*sigma) ;
      m_psi_W_EL(i,j) += m_system->get_gibbsfactor()* X(i) * logistic(-v_j(j,X))  /(sigma*sigma) * m_DeltaE;
    }
  }


}

double Sampler::logistic(double x){
  return (1.0/ (1.0+exp(x) ) ) ;
}

double Sampler::v_j(int j,Eigen::VectorXd X){
  double sum = 0;
    for (int i=0; i<m_system->getRBM()->get_M();i++)
    {
      sum+= X(i)* m_system->getRBM()->get_W()(i,j) ;
    }
  //sum =  X.dot(m_W.col(j)) ;
  return m_system->getRBM()->get_b()(j)+sum / (m_system->getRBM()->get_sigma() * m_system->getRBM()->get_sigma()) ;
  }

Eigen::VectorXd Sampler::Return_gradients_a(){

  return  2*(m_psi_a_EL/m_stepNumber-m_energy*m_psi_a/m_stepNumber);
}
Eigen::VectorXd Sampler::Return_gradients_b(){

  return  2*(m_psi_b_EL/m_stepNumber-m_energy*m_psi_b/m_stepNumber);
}

Eigen::MatrixXd Sampler::Return_gradients_W(){

  return  2*(m_psi_W_EL/m_stepNumber-m_energy*m_psi_W/m_stepNumber);
}
