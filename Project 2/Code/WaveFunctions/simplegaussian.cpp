#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

//remove
#include <iostream>
using namespace std;
//

SimpleGaussian::SimpleGaussian(System* system, int numberOfDimensions,double alpha, double beta) :
        WaveFunction(system,numberOfDimensions) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */


    double r_squared,g=1.0,f=1.0;

    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      r_squared=0;
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        if (dim==2){r_squared+=m_parameters[1]*(particles[i]->getPosition()[dim])*(particles[i]->getPosition()[dim]);}
        else{r_squared+= (particles[i]->getPosition()[dim])*(particles[i]->getPosition()[dim]);}
        // if dimension == 2 : beta
      }

    g=g*exp(-m_parameters[0]*r_squared);  //exp ( alpha * r^2)
    }
    return g*f;
}


std::vector<double> SimpleGaussian::updateQForce(std::vector<class Particle*> particles,int particle) {
    for (int dim=0; dim < m_system->getNumberOfDimensions();dim++)
    {
    m_qForce.at(dim)=-4*m_parameters[0] * particles[particle]->getPosition()[dim]  ;
    }

    // F_i = -4 alpha r_i
    return m_qForce;
}
