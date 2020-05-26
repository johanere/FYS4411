#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"



using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...

    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    return kineticEnergy + potentialEnergy;
    */
    std::vector<double> paramters=(m_system->getWaveFunction())->getParameters();

    int n_particles=m_system->getNumberOfParticles();
    double alpha = paramters[0];
    double beta = paramters[1];
    double r_squared=0;

    for (int i=0; i<n_particles;i++)
    {
      for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++)
      {
        if (dim==2) {r_squared+=(paramters[1]*paramters[1]* particles[i]->getPosition().at(dim) )*( particles[i]->getPosition().at(dim) );}
        else {r_squared+=( particles[i]->getPosition().at(dim) )*( particles[i]->getPosition().at(dim) );}
        // if dimension == 2 : beta
      }
    }

    double EL=alpha*(m_system->getNumberOfDimensions()-1+beta)*n_particles
    +(0.5*m_omega*m_omega-2*alpha*alpha)*r_squared;
    return EL;
}
