#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system,int GD_iters) {
    m_system = system;
    m_stepNumber = 0;
    m_acceptedSteps = 0;
    m_GDiters=GD_iters;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
    m_EnergySamples.reserve(steps);
}

void Sampler::sample(bool acceptedStep) {
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergy2 = 0;
    }

    m_DeltaE=m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());


    m_cumulativeEnergy  += m_DeltaE;
    m_cumulativeEnergy2 += m_DeltaE*m_DeltaE;
    m_stepNumber++;
    m_EnergySamples.push_back(m_DeltaE);
    if (acceptedStep==true){m_acceptedSteps++;}
    else {;}

    if (m_GDiters>1)
    {
      double r2;
      for (int i=0;i<m_system->getNumberOfParticles();i++)
      {
        r2=0;
        for (int dim=0;dim<m_system->getNumberOfDimensions();dim++)
        {
          if (dim==2){r2-=(m_system->getParticles().at(i))->getPosition().at(dim)
            * (m_system->getParticles().at(i))->getPosition().at(dim);} // add *times parameter[1]
          else{r2-=(m_system->getParticles().at(i))->getPosition().at(dim)
            * (m_system->getParticles().at(i))->getPosition().at(dim);}


        }
      }
        m_ElR+=m_DeltaE*(r2);
        m_sumR2+=r2;
    }
}

void Sampler::printOutputToTerminal() {
    cout << "alpha: " << m_system->getWaveFunction()->getParameters().at(0) << endl;
    cout << "E:     " << m_energy << endl;
    cout << "E/N:   "<<m_energy/m_system->getNumberOfParticles()<<endl;
    cout << std::scientific;
    cout << "S:     " << m_error << endl;
    cout << std::fixed;
    cout << "A:     " <<  ((double) m_acceptedSteps)/m_stepNumber << endl;
}

void Sampler::computeAverages() {
    m_energy = m_cumulativeEnergy/m_stepNumber;
    m_energy2 =  m_cumulativeEnergy2/m_stepNumber;
    m_variance=m_energy2 - (m_energy*m_energy);
    m_error=sqrt(m_variance/m_stepNumber);
    cout<<"\n"<<endl; //FJERNE

    if (m_GDiters>1)
    {
      m_gradAlpha=2*( (m_ElR/m_stepNumber)-(m_energy)*(m_sumR2/m_stepNumber) );
    }
    m_A= ((double) m_acceptedSteps)/((double) m_stepNumber);
}
