#include <cmath>
#include <cassert>
#include "Wavefunction.h"
#include "system.h"

#include <iostream>
using namespace std;


Wavefunction::Wavefunction(System* system){
m_system  =system;

m_a = m_system->getRBM()-> get_a();
m_b = m_system->getRBM()-> get_b();
m_W = m_system->getRBM()-> get_W();

m_sigma = m_system->getRBM()-> get_sigma();
m_sigma4=m_sigma*m_sigma*m_sigma*m_sigma;

m_M=m_system->getRBM()->get_M();
m_N=m_system->getRBM()->get_N();


} // end of constructor

Wavefunction::~Wavefunction(void) {
}

double Wavefunction::evaluate(Eigen::VectorXd X) {
  double exponent1=0;
  double f=1.0;
  for (int i=0; i<m_M;i++)
  {
    exponent1+=(X(i)-m_a(i))*(X(i)-m_a(i)) / (2*m_sigma*m_sigma);
  }
  for (int j=0; j<m_N;j++)
  {
    f *= (1.0+exp( v_j(j,X)   )) ;
  }

  return exp(-exponent1)*f;
} // end of evaluate


double Wavefunction::logistic(double x){
  return (1.0/ (1.0+exp(x) ) ) ;
}

double Wavefunction::v_j(int j,Eigen::VectorXd X){
  double sum = 0;
  for (int i=0; i<m_M;i++)
  {
    sum+= X(i)* m_W(i,j) ;
  }
  //sum =  X.dot(m_W.col(j)) ;
  return m_b(j)+sum / (m_sigma*m_sigma) ;
  }


double Wavefunction::getQForce(Eigen::VectorXd X,int k) {
  double sum=0;
  for (int j=0; j<m_N;j++)
  {
    sum+= m_W(k,j)*logistic( - v_j(j,X) );
  }

  return -2*(X(k)-m_a(k))/(m_sigma*m_sigma) + 2*sum/(m_sigma*m_sigma);
}
