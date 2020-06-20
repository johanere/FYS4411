#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include "RBM.h"
#include "Hamiltonians/H_RBM_Simple.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system,int GD_iters) {
    m_system = system;
    m_stepNumber = 0;
    m_acceptedSteps = 0;
    m_GDiters=GD_iters;
    m_sigma = m_system->getRBM() ->get_sigma();

    for (int i=0; i<m_system->getNumberOfParticles();i++)
    {
      for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
      {
        m_X.push_back(0.0);
      }
    }
    m_a = m_system->getRBM() -> get_a();
    m_b = m_system->getRBM() -> get_b();
    m_W = m_system->getRBM() -> get_W();


    // fill in zeros in all gradients IOT accumulate later
    for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
    {
      m_psi_a.push_back( 0.0);
      m_psi_a_EL.push_back( 0.0);
    }

    for (int j=0; j<m_system->getRBM() ->get_n() ;j++)
    {
      m_psi_b.push_back( 0.0 );
      m_psi_b_EL.push_back(0.0 );
    }

    for (int k=0; k<m_system->getRBM() ->get_m() ;k++)
    {
      for (int l=0; l<m_system->getRBM() ->get_n() ;l++)
      {
        m_psi_W.push_back( 0.0 )     ;
        m_psi_W_EL.push_back( 0.0  ) ;
      }
    }

}

Sampler::~Sampler(void) {
   cout << "Sampler is being deleted" << endl;
}

void Sampler::setNumberOfMetropolisSteps(int steps, int steps_after_eq) {
    m_numberOfMetropolisSteps = steps;
    cout<<"Fjernet m_EnergySamples"<<endl;
    //m_EnergySamples.reserve(steps_after_eq);


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
  //m_EnergySamples.push_back(m_DeltaE);
    cout<<m_DeltaE<<endl;
    Gradient_in_step(m_DeltaE);

    if (acceptedStep==true){m_acceptedSteps++;}
    else {;}

    if (m_GDiters>1)
    {

    }
}

void Sampler::printOutputToTerminal() {
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

double Sampler::compute_v_j(int j) {
  double sum=0;


  for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
  {
    sum+=m_X[i]*m_W[ (i* m_system->getRBM() ->get_m() )  +j]/(m_sigma*m_sigma);
  }
  return m_b[j]+sum ;
} // end of v_j



void Sampler::Gradient_in_step(double m_DeltaE) {
  //update X
  for (int i=0; i<m_system->getNumberOfParticles();i++)
  {
    for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
    {
      m_X[i+dim]= m_system->getParticles()[i]->getPosition()[dim];
    }
  }



  for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
  {
    m_psi_a[i]      +=     (m_W[i]-m_a[i]) / (m_sigma*m_sigma)  ;
    m_psi_a_EL[i]   +=     (m_W[i]-m_a[i]) / (m_sigma*m_sigma) * m_DeltaE ;
  }

  for (int j=0; j<m_system->getRBM() ->get_n() ;j++)
  {
    m_psi_b[j]      +=    1.0 / (1.0 + exp(-compute_v_j(j)) )  ;
    m_psi_b_EL[j]   +=     1.0 / (1.0 + exp(-compute_v_j(j)) ) * m_DeltaE ;
  }
  for (int k=0; k<m_system->getRBM() ->get_m() ;k++)
  {
    for (int l=0; l<m_system->getRBM() ->get_n() ;l++)
    {
      m_psi_W[k*    (m_system->getRBM() ->get_m())+l ]     +=  m_X[k] / ((1.0 + exp(-compute_v_j(l)) ) * m_sigma*m_sigma )  ;
      m_psi_W_EL[k* (m_system->getRBM() ->get_m())+l ]  +=  m_X[k] / ((1.0 + exp(-compute_v_j(l)) ) * m_sigma*m_sigma ) * m_DeltaE ;
    }
  }
} // End of Gradient_in_step

void Sampler::Return_gradients(){
  cout<<"START I RET GRAD"<<endl;
  //Take averages
  for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
  {
    m_psi_a[i]      = m_psi_a[i]     /m_stepNumber   ;
    m_psi_a_EL[i]   = m_psi_a_EL[i]  /m_stepNumber    ;
    m_a[i]=( 2*(m_psi_a_EL[i]-m_psi_a[i]*m_energy));
  }

  for (int j=0; j<m_system->getRBM() ->get_n() ;j++)
  {
    m_psi_b[j]      =    m_psi_b[j]     /m_stepNumber ;
    m_psi_b_EL[j]   =    m_psi_b_EL[j]  /m_stepNumber ;
    m_b[j]          = ( 2*(m_psi_b_EL[j]-m_psi_b[j]*m_energy));
  }
  for (int k=0; k<m_system->getRBM() ->get_m() ;k++)
  {
    for (int l=0; l<m_system->getRBM() ->get_n() ;l++)
    {
      m_psi_W[k*    (m_system->getRBM() ->get_m())+l ]     = m_psi_W[k*    (m_system->getRBM() ->get_m())+l ]  /m_stepNumber ;
      m_psi_W_EL[k* (m_system->getRBM() ->get_m())+l ]     = m_psi_W_EL[k* (m_system->getRBM() ->get_m())+l ]  /m_stepNumber ;
      m_W[k* (m_system->getRBM() ->get_m())+l ]= ( 2*(m_psi_W_EL[k*    (m_system->getRBM() ->get_m())+l]-m_psi_W[k*    (m_system->getRBM() ->get_m())+l]*m_energy));
    }
  }

  cout<<"FERDIG I RET GRAD"<<endl;
}


/*
void Sampler::Return_gradient_a(double m_DeltaE) {
std::vector <double>  grad_a;
for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
{
  grad_a
}

} // End of Gradient_in_step

std::vector <double> get_exp_val_grad(std::vector <double>  psi_beta_EL,
   std::vector <double>  psi_beta_EL) {

std::vector <double>  grad_beta;
for (int l=0; l< (int) psi_beta_EL.size() ;l++)
{
  grad_beta.push_back()
}


} */
