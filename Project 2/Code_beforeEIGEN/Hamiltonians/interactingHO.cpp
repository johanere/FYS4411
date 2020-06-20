#include "interactingHO.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

#include <cmath>

using std::cout;
using std::endl;

interactingHO::interactingHO(System* system, double omega,double a) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_a=a;
}

double interactingHO::computeLocalEnergy(std::vector<Particle*> particles) {

    std::vector<double> paramters=(m_system->getWaveFunction())->getParameters();

    int n_particles=m_system->getNumberOfParticles();

    assert(m_system->getNumberOfDimensions() ==3);

    double alpha = paramters[0];
    double beta  = paramters[1];
    double sum = 0;

    double term1=0; // del_k^2 phi_k / phi_k
    double term2=0; // del_k phi_k / phi_k sum (rk-rj)/rkj*u'kj
    double term3=0; // sum (rk-rj)/rkj*u'kj
    double term4=0;// sum (u''kj + 2/rkj * u'kj)

    double rx=0;
    double ry=0;
    double rz=0;

    double delta_rx=0;
    double delta_ry=0;
    double delta_rz=0;

    double r_kj=0;
    double du_kj= 0; //modified to include 1/rkl

    double rx_j=0;
    double ry_j=0;
    double rz_j=0;

    for (int k=0; k<n_particles;k++)
    {

       term1=0; // del_k^2 phi_k / phi_k
       term2=0; // del_k phi_k / phi_k sum (rk-rj)/rkj*u'kj
       term3=0; // sum (rk-rj)/rkj*u'kj
       term4=0;// sum (u''kj + 2/rkj * u'kj)

       rx=particles[k]->getPosition().at(0);
       ry=particles[k]->getPosition().at(1);
       rz=particles[k]->getPosition().at(2);

       delta_rx=0;
       delta_ry=0;
       delta_rz=0;

      term1 = -2*alpha*(2 + beta ) + 4*alpha*alpha*( (rx*rx)+(ry*ry)+(beta*beta*rz*rz) ); //del_k^2 phi_k / phi_k
      for (int j=0; j<n_particles;j++)
      {
        if (j !=  k )
        {
          r_kj=m_system->getR_jk(k,j);
          du_kj= m_a/(r_kj *r_kj*(r_kj-m_a) ); //modified to include 1/rkl

          rx_j=particles[j]->getPosition().at(0);
          ry_j=particles[j]->getPosition().at(1);
          rz_j=particles[j]->getPosition().at(2);

          delta_rx+=(rx-rx_j)*du_kj;
          delta_ry+=(ry-ry_j)*du_kj;
          delta_rz+=(rz-rz_j)*du_kj;

          term4 +=  (  (m_a*m_a) - (2*m_a*r_kj)  ) / (  (r_kj*r_kj) * (r_kj-m_a) * (r_kj-m_a)  ) + 2*du_kj;

        }
        else {;}

      }
      term2 = 2*(-2*alpha*( (rx*delta_rx) + (ry*delta_ry) + (beta*rz*delta_rz) ) );
      term3 = (delta_rx*delta_rx)+(delta_ry*delta_ry)+(delta_rz*delta_rz);
      sum += -(term1+term2+term3+term4)+(rx*rx)+(ry*ry)+(beta*beta*rz*rz);
    }


    double EL=0.5*sum;

    if (EL<0){cout<<EL<<endl;}

    return EL;
}
