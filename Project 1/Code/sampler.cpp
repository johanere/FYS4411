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
    // Make sure the sampling variable(s) are initialized at the first step.
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
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();
    /*
    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 2^" << std::log2(ms) << endl;
    cout << " Number of equilibration steps  : 2^" << std::log2(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " std : " << m_error << endl;
    cout << "accepted steps " << ((double) m_acceptedSteps)/m_stepNumber<<endl;
    cout << endl; */
    cout << "alpha: " << pa.at(0) << endl;
    cout << "E:     " << m_energy << endl;
    cout << "E/N:   "<<m_energy/m_system->getNumberOfParticles()<<endl;
    cout << std::scientific;
    cout << "S:     " << m_error << endl;
    cout << std::fixed;
    cout << "A:     " <<  ((double) m_acceptedSteps)/m_stepNumber << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
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
