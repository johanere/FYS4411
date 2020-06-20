#include <cassert>
#include <iostream>

#include "Hamiltonian.h"
#include "system.h"

#include "Eigen/Dense"

#include <cmath>

using std::cout;
using std::endl;

Hamiltonian::Hamiltonian(System* system) {
          m_system=system;

          m_a = m_system->getRBM()-> get_a();
          m_b = m_system->getRBM()-> get_b();
          m_W = m_system->getRBM()-> get_W();

          m_sigma = m_system->getRBM()-> get_sigma();
          m_sigma4=m_sigma*m_sigma*m_sigma*m_sigma;

          m_M=m_system->getRBM()->get_M();
          m_N=m_system->getRBM()->get_N();
          m_omega= m_system->getRBM()->get_omega()  ;


} // end of constructor

Hamiltonian::~Hamiltonian(void) {
}


double Hamiltonian::computeLocalEnergy(Eigen::VectorXd X) {
    double e_K=0,e_P=0;
    double term1;

    for (int i=0; i<m_M;i++)
    {
      term1=nabla(i,X);
      e_K+= term1*term1+nabla2(i,X);
    ;
    }

    for (int i=0; i<m_M;i++)
    {
      e_P+=X(i)*X(i)*m_omega*m_omega;
    }
    return -0.5*e_K + 0.5*e_P;
} // end of local energy



double Hamiltonian::logistic(double x){
  return (1.0/ (1.0+exp(x) ) ) ;
}

double Hamiltonian::v_j(int j,Eigen::VectorXd X){
  double sum = 0;
  for (int i=0; i<m_M;i++)
  {
    sum+= X(i)* m_W(i,j) ;
  }
  //sum =  X.dot(m_W.col(j)) ;
  return m_b(j)+sum / (m_sigma*m_sigma) ;
}

double Hamiltonian::nabla(int k,Eigen::VectorXd X){
  double sum=0;
  for (int j=0; j<m_N;j++)
  {
    sum+= m_W(k,j)*logistic( - v_j(j,X) );
  }

  return -(X(k)-m_a(k))/(m_sigma*m_sigma) + sum/(m_sigma*m_sigma);
}

double Hamiltonian::nabla2(int k,Eigen::VectorXd X){
  double sum=0;
  for (int j=0; j<m_N;j++)
  {
    sum+= m_W(k,j)*m_W(k,j)*logistic( - v_j(j,X) ) * logistic( - v_j(j,X) ) * exp(-v_j(j,X));
  }
  return -1.0/(m_sigma*m_sigma) +sum/m_sigma4;
}
