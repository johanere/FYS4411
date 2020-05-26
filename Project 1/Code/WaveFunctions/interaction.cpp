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

interactingWF::interactingWF(System* system, int numberOfDimensions,double alpha,double beta,double a) :
        WaveFunction(system,numberOfDimensions) {
    assert(alpha >= 0);
    assert(numberOfDimensions == 3);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_a=a;
}

double interactingWF::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    //single particle
    assert(m_system->getNumberOfDimensions() == 3);

    double r_squared,g=1.0,f=1.0;
    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      r_squared=0;
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        if (dim==2){r_squared+= m_parameters[1]*(particles[i]->getPosition()[dim])*(particles[i]->getPosition()[dim]);}
        else{r_squared+= (particles[i]->getPosition()[dim])*(particles[i]->getPosition()[dim]);}
        // if dimension == 2 : beta
      }
      g=g*exp(-m_parameters[0]*r_squared);  //exp ( alpha * r^2)
    }

    //interaction
    double r_jk;

    for (int j=0; j<m_system->getNumberOfParticles();j++)
    {
      for (int k=j+1; k<m_system->getNumberOfParticles();k++)
      {
        r_jk=m_system->getR_jk(j,k);
        if(r_jk>m_a){
        f=f*( 1 - m_a/m_system->getR_jk(j,k) );}
        else {return 0; }
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
      double factor=0,r_kj=0,phi_dk=0;


      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        m_qForce.at(dim) =0;
      }


      for (int j=0; j<m_system->getNumberOfParticles();j++)
      {

        if (j != int( particle) )
        {
          r_kj= m_system->getR_jk(particle,j);


          factor = m_a/(r_kj *r_kj*(r_kj-m_a) );
          for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
          {
            m_qForce.at(dim) += ( (particles[particle]->getPosition()[dim] - particles[j]->getPosition()[dim]) ) * factor;
          }
        }
        else {;}
      }

      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {

      if (dim==2){phi_dk=-2*m_parameters[0] * m_parameters[1] * particles[particle]->getPosition()[2];}
      else{phi_dk= -2*m_parameters[0] * particles[particle]->getPosition()[dim];}
      m_qForce.at(dim)=2*(phi_dk+m_qForce.at(dim));
      }

     return m_qForce;
}
