#include "interactingHO.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

#include <cmath>

using std::cout;
using std::endl;

H_RBM_Simple::interactingHO(System* system, double omega,double a) :
        Hamiltonian(system) {
    assert(omega > 0);
    // set number of nodes
    // initiate W, b, a
}

double H_RBM_Simple::computeLocalEnergy(std::vector<Particle*> particles) {
    sigma =1
    m_a =1

    std::vector<double> paramters=(m_system->getWaveFunction())->getParameters();

    int n_particles=m_system->getNumberOfParticles();


    for (int i=0; i<m_system->getM();i++)  // loop over visible nodes
    {
      term1=0
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)  // a_k-X_k/sigma
      {
        term1+=( (particles[i]->getPosition()[dim]) -m_a[i+dim] ) / sigma
      }

      for (int j=0; j<n;j++) //w_kj / sigma (1+exp(-bj - ))
      {
        for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
        {
          term1+=( (particles[i]->getPosition()[dim]) -m_a[i+dim] ) / sigma
        }

      }
    }



    double EL=0.5*sum;

    if (EL<0){cout<<EL<<endl;}

    return EL;
}
