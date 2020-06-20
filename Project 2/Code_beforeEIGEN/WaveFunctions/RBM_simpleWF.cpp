#include "RBM_SimpleWF.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

//remove
#include <iostream>
using namespace std;
//

RBM_SimpleWF::RBM_SimpleWF(System* system,  int numberOfDimensions,double omega,double sigma) :
    WaveFunction(system,numberOfDimensions) {

      for (int i=0; i<m_system->getNumberOfParticles();i++)
      {
        for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
        {
          m_X.push_back(0.0);
        }
      }

      m_M = ( m_system->getRBM() )-> get_m() ;
      m_N = ( m_system->getRBM() )-> get_n() ;


      assert(m_M == (int) m_X.size());

      m_a = ( m_system->getRBM() )-> get_a() ;
      m_b = ( m_system->getRBM() )-> get_b() ;
      m_W = ( m_system->getRBM() )-> get_W() ;
      m_sigma = sigma ;
      m_omega = omega ;
} // end of constructor


double RBM_SimpleWF::evaluate(std::vector<class Particle*> particles) {
    double exponent1=0;
    double f=1.0;

    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        m_X[i+dim]= particles[i]->getPosition()[dim];
      }
    }

    //first exp factor
    for (int i=0; i<m_M;i++)
    {
        exponent1+= ( m_X[i]-m_a[i] )*( m_X[i]-m_a[i]  ) / (2*m_sigma*m_sigma) ;
    }
    //second exp factor
    for (int j=0; j<m_N;j++)
    {
      f = f * ( 1.0 + exp( compute_v_j(j) ) ) ;
    }

    // the whole shabang
    return exp(-exponent1)*f;

} // end of evaluate

double RBM_SimpleWF::compute_v_j(int j) {
  double sum=0;
  for (int i=0; i<m_M;i++)
  {
    sum+=m_X[i]*m_W[i*m_M+j]/(m_sigma*m_sigma);
  }
  return m_b[j]+sum ;
} // end of v_j


std::vector<double> RBM_SimpleWF::updateQForce(std::vector<class Particle*> particles,int particle) {

    for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        m_qForce.at(dim) =0;
      }



     return m_qForce;
}
