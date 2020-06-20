#include "H_RBM_Simple.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

#include <cmath>

using std::cout;
using std::endl;

H_RBM_Simple::H_RBM_Simple(System* system,double omega,double sigma) :
        Hamiltonian(system) {

    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        m_X.push_back(0.0);
      }
    }

    m_M = m_system ->getRBM()  -> get_m() ;
    m_N = m_system ->getRBM()  -> get_n() ;


    assert(m_M == (int) m_X.size());

    m_a =  m_system->getRBM() -> get_a() ;
    m_b =  m_system->getRBM() -> get_b() ;
    m_W =  m_system->getRBM() -> get_W() ;

    m_sigma = sigma;
    m_omega = omega;

    cout<<m_W[0]<<endl;
} // end of constructor

double H_RBM_Simple::computeLocalEnergy(std::vector<Particle*> particles) {
    double nabla,nabla2,E_K=0,E_P=0;

    //update X vector
    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        m_X[i+dim]= particles[i]->getPosition()[dim];
      }
    }




    for (int k=0; k<m_M;k++)  // loop over visible nodes
    {
      nabla=compute_nabla_ln_psi(k);
      nabla2=compute_nabla2_ln_psi(k);
      E_K+=nabla*nabla+nabla2 ;
      E_P+=m_X[k]*m_X[k] ;
    }



    double EL= -0.5 * E_K  + 0.5 * m_omega * m_omega * E_P;

    return EL;
} // end of local energy


double H_RBM_Simple::compute_v_j(int j) {
  double sum=0;
  for (int i=0; i<m_M;i++)
  {
    sum+=m_X[i]*m_W[i*m_M+j]/(m_sigma*m_sigma);
      //cout<<"asdasd"<<m_W[i*m_M+j]<<endl;;
  }

  return m_b[j]+sum ;
} // end of v_j



double H_RBM_Simple::compute_nabla_ln_psi(int k) {
  double sum=0;
  for (int j=0; j<m_N;j++) //w_kj / sigma (1+exp(-bj - ))
  {
    sum+= m_W[k*m_M+j]/(m_sigma*m_sigma) * (1.0/( 1.0 + exp(- compute_v_j(j) ) ));
  }
  return -(m_X[k]-m_a[k])/(m_sigma*m_sigma)+sum;
} //end of nabla ln psi



double H_RBM_Simple::compute_nabla2_ln_psi(int k) {
  double sum=0;
  double exponential;
  for (int j=0; j<m_M;j++) //w_kj / sigma (1+exp(-bj - ))
  {
    exponential =  exp(- compute_v_j(j) ) ;
    sum+= m_W[k*m_M+j]*m_W[k*m_M+j]/(m_sigma*m_sigma*m_sigma*m_sigma) *
    exponential/( (1+exponential)*(1+exponential) )  ;
  }
  cout<<"EK i H "<<exponential/( (1+exponential)*(1+exponential) )<<endl;
  return -1/(m_sigma*m_sigma)+sum ;
} // end of nabla2 ln psi
