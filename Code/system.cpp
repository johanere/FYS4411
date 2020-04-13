#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"



//remove
#include <iostream>
using std::cout;
using std::endl;
//

System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    assert(m_stepLength>0);


    //pick random particle
    particle=m_random->nextInt(0,m_numberOfParticles-1);

    // make copy of old pos
    oldPosition=m_particles[particle]->getPosition();


    for (int dimension=0; dimension < m_numberOfDimensions; dimension++)
    {
      //move particle
      m_particles[particle]->adjustPosition((m_random->nextDouble()-0.5)*m_stepLength,dimension);
    }
    //calculate new WF
    newWF = m_waveFunction->evaluate(m_particles);

    // check if move in accepted
    if(m_random->nextDouble() > (newWF*newWF) / (oldWF*oldWF)) {
      //if not accepted, reset particle
      m_particles[particle]->setPosition(oldPosition);
      return false;
    }
    else  {
      oldWF = newWF;
      return true;
    }

}


bool System::metropolishastingsStep() {
    assert(m_stepLength>0);
    //pick random particle
    particle=m_random->nextInt(0,m_numberOfParticles-1);
    // make copy of old pos
    oldPosition=m_particles[particle]->getPosition();


    std::vector<double> oldQForce = m_waveFunction->updateQForce(m_particles,particle);

    for (int dim=0; dim < m_numberOfDimensions; dim++)
    {
      //move particle
      //  adjust D*qForce[dim]*delta_t*eps*sqrt(delta_t), dim
      m_particles[particle]->adjustPosition(0.5*oldQForce.at(dim)*m_stepLength
      +m_random->nextGaussian(0,1)*m_stepLengthsqrt,dim);
    }
    //calculate new WF
    newWF = m_waveFunction->evaluate(m_particles);

    std::vector<double> newQForce = m_waveFunction->updateQForce(m_particles,particle);
    
    //Greens
    double gfRatio;
    double term1;
    double term2;
    for (int dim=0; dim < m_numberOfDimensions; dim++)
    {
        term1 = oldPosition.at(dim) - m_particles[particle]->getPosition().at(dim) - 0.5 *m_stepLength*newQForce.at(dim);
        term2 = m_particles[particle]->getPosition().at(dim) - oldPosition.at(dim) - 0.5 *m_stepLength*oldQForce.at(dim);
        gfRatio += (term2*term2) - (term1*term1);
    }

    gfRatio=gfRatio/2*m_stepLength;
    gfRatio = exp(gfRatio);
    // check if move in accepted
    if(m_random->nextDouble() > gfRatio*(newWF*newWF) / (oldWF*oldWF)) {
      //if not accepted, reset particle
      m_particles[particle]->setPosition(oldPosition);
      return false;
    }
    else  {
      oldWF = newWF;
      return true;
    }

}

void System::runMetropolisSteps(int numberOfMetropolisSteps, int method) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    oldWF = m_waveFunction->evaluate(m_particles);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
      if (method == 0) acceptedStep = metropolisStep();
      else if (method == 1) acceptedStep = metropolishastingsStep();
      else {
        cout<<"Provide valid method selection"<<endl;
        abort();
      }

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
         if (i>=m_equilibrationFraction*numberOfMetropolisSteps)
          {
          m_sampler->sample(acceptedStep);
          }
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
    m_stepLengthsqrt=sqrt(m_stepLength);
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}
