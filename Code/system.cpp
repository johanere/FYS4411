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
bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
     assert(m_stepLength>0);
     //get old WF
    float oldWF = m_waveFunction->evaluate(m_particles);

    //pick random particle
    int particle=Random::nextInt(m_numberOfParticles);
    // make copy of old pos
    std::vector<double> oldPosition=m_particles[particle]->getPosition();

    for (int dimension=0; dimension < m_numberOfDimensions; dimension++)
    {
      //move particle
      m_particles[particle]->adjustPosition((Random::nextDouble()-0.5)*m_stepLength,dimension);
    }
    //calculate new WF
    float newWF = m_waveFunction->evaluate(m_particles);


    // check if move in accepted
    double q= (newWF*newWF) / (oldWF*oldWF);
    if(Random::nextDouble() > q) {
      //if not accepted, reset particle
      m_particles[particle]->setPosition(oldPosition);
      return false;
    }
    else  {
      return true;
    }

}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

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
