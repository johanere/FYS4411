#include "interaction.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

//remove
#include <iostream>
using namespace std;
//

RBM_SimpleWF::RBM_SimpleWF(System* system,  int numberOfDimensions, std::vector<double> w, std::vector<double> a, std::vector<double> b, double sigma) :
    WaveFunction(system,numberOfDimensions) {
    // m_parameters.reserve(1);
    m_w=w;
    m_b=b;
    m_a=a;
    m_sigma=sigma
}

double RBM_SimpleWF::evaluate(std::vector<class Particle*> particles) {

    double exponent1=0,exponent2=0
    double f=1.0;
    int m,n;

    //first exp factor
    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        exponent1+= ( (particles[i]->getPosition()[dim]) -m_a[i+dim] )*( (particles[i]->getPosition()[dim]) -m_a[i+dim] ) / (2*sigma*sigma)
      }
    }

    //second exp factor
    m=m_system->getM();
    n=m_system->getN();
    for (int j=0; j<n;j++)
    {
      exponent2=0
      for (int i=0; i<m_system->getNumberOfParticles();i++)
      {
        for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
        {
          exponent2+= ( (particles[i]->getPosition()[dim]) * w[i*m+j] ) / m_sigma
        }
      }
      f=f*( 1.0+exp(b[j]+exponent2) )
    }

    // the whole shabang
    return exp(-exponent1)*f;

}
