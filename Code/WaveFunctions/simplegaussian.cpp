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

SimpleGaussian::SimpleGaussian(System* system, int numberOfDimensions,double alpha) :
        WaveFunction(system,numberOfDimensions) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
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
        r_squared+= (particles[i]->getPosition()[dim])*(particles[i]->getPosition()[dim]);
        // if dimension == 2 : beta
      }
    g=g*exp(-m_parameters[0]*r_squared);  //exp ( alpha * r^2)
    }
    return g*f;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    cout<<" DD DD  DD  DD  DD  DD "<<endl;
    return 0;
}

std::vector<double> SimpleGaussian::updateQForce(std::vector<class Particle*> particles,int particle) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    for (int dim=0; dim < m_system->getNumberOfDimensions();dim++)
    {
    m_qForce.at(dim)=-4*m_parameters[0] * particles[particle]->getPosition()[dim]  ;
    }
    // F_i = -4 alpha r_i
    return m_qForce;
}
