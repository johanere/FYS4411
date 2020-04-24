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

interactingWF::interactingWF(System* system, int numberOfDimensions,double alpha) :
        WaveFunction(system,numberOfDimensions) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double interactingWF::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    assert(m_system->getNumberOfDimensions() == 3);

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

    //decleare variables
    double a=1;
    double rj_0,rj_1,rj_2,rk_0,rk_1,rk_2,r_jk,u;
    for (int j=0; j<m_system->getNumberOfParticles();j++)
    {
      rj_0=particles[j]->getPosition().at(0);
      rj_1=particles[j]->getPosition().at(1);
      rj_2=particles[j]->getPosition().at(2);
      for (int k=j+1; k<m_system->getNumberOfParticles();k++)
      {
        rk_0=particles[k]->getPosition().at(0);
        rk_1=particles[k]->getPosition().at(1);
        rk_2=particles[k]->getPosition().at(2);
        r_jk=sqrt((rj_0-rk_0)*(rj_0-rk_0)+(rj_1-rk_1)*(rj_1-rk_1)+(rj_2-rk_2)*(rj_2-rk_2));
        if (r_jk==0) {
          return 0;
        }
        else{
          f=f*(1-a/(r_jk));
        }
      }

    }
    return g*f;
}


double interactingWF::computeDoubleDerivative(std::vector<class Particle*> particles) {
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

std::vector<double> interactingWF::updateQForce(std::vector<class Particle*> particles,int particle) {
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
