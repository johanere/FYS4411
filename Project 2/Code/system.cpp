#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"

#include "RBM.h"

//#include "WaveFunctions/RBM_SimpleWF.h"
//#include "Hamiltonians/H_RBM_Simple.h"

#include <math.h>


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

bool System::brute_force_Step() {
    assert(m_stepLength>0);
    //pick random particle
    particle=m_random->nextInt(0,m_numberOfParticles-1);
    // make copy of old pos
    oldPosition=m_particles[particle]->getPosition();
    //move particle
    for (int dimension=0; dimension < m_numberOfDimensions; dimension++){
      m_particles[particle]->adjustPosition(2*(m_random->nextDouble()-0.5)*m_stepLength,dimension);
    }

    //calculate new distances for moved particle
        if (m_interaction>0){
      calculateR_jk(m_particles, particle);
    }
    //calculate new WF
    newWF = m_waveFunction->evaluate(m_particles);

    if (m_interaction>0){
        if (m_break_eval>0){
        m_break_eval=0;
        m_particles[particle]->setPosition(oldPosition);
        calculateR_jk(m_particles, particle);
        return false;}
    }
    // check if move in accepted, if not accepted, reset particle
    if(m_random->nextDouble() > (newWF*newWF) / (oldWF*oldWF)) {
      m_particles[particle]->setPosition(oldPosition);
      //reset distances
      if (m_interaction>0){
        calculateR_jk(m_particles, particle);
      }
      return false;
    }
    else  {
      oldWF = newWF;
      return true;
    }
} //end of brute force sampling


bool System::importance_sampling_Step() {
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
    if (m_interaction>0){
      calculateR_jk(m_particles, particle);
    }
    //calculate new WF
    newWF = m_waveFunction->evaluate(m_particles);

    //auto reject new proposal if any particles come with a distance less than a
    if (m_interaction>0){
        if (m_break_eval>0){
        m_break_eval=0;
        m_particles[particle]->setPosition(oldPosition);
        calculateR_jk(m_particles, particle);
        return false;}
    }


    std::vector<double> newQForce = m_waveFunction->updateQForce(m_particles,particle);

    //Greens
    double gfRatio=0;
    double term1=0;
    double term2=0;

    for (int dim=0; dim < m_numberOfDimensions; dim++)
    {
        term1 = oldPosition.at(dim) - m_particles[particle]->getPosition().at(dim) - 0.5 *m_stepLength*newQForce.at(dim);
        term2 = m_particles[particle]->getPosition().at(dim) - oldPosition.at(dim) - 0.5 *m_stepLength*oldQForce.at(dim);
        gfRatio +=  - (term1*term1) + (term2*term2);
        gfRatio= gfRatio/(2.0*m_stepLength);
        gfRatio = exp(gfRatio);
    }

    if(m_random->nextDouble() > gfRatio*(newWF*newWF)  / (oldWF*oldWF)) {
      //if not accepted, reset particle
      m_particles[particle]->setPosition(oldPosition);
      return false;
      if (m_interaction>0){
        calculateR_jk(m_particles, particle);
      }
    }
    else  {
      oldWF = newWF;
      return true;
    }

} //end of improtance sampling step

void System::runMetropolisSteps(int numberOfMetropolisSteps, int method,int GD_iters) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this,GD_iters);

    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    InitiateR(); //calculate distances

    oldWF = m_waveFunction->evaluate(m_particles);


    for (int i=0; i < numberOfMetropolisSteps; i++) {
      if (method == 0) acceptedStep = brute_force_Step();
      else if (method == 1) acceptedStep = importance_sampling_Step();
      else {
        cout<<"Provide valid method selection"<<endl;
        abort();
      }
         if (i>=m_equilibrationFraction*numberOfMetropolisSteps)
          {
          m_sampler->sample(acceptedStep);
          }
    }

    m_sampler->computeAverages();
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

void System::setDistances(int interaction,double a) {
  m_interaction=interaction;
  m_a=a;
}

void System::setRBM(RBM* rbm) {
    m_rbm = rbm;
    m_M=m_rbm -> get_m() ;
    m_N=m_rbm -> get_n();
    m_rbm -> InitiateWeightsAndBiases();
}

void System::InitiateR( ){
  if (m_interaction>0)
  {
  // Calculates the distance between all particles and stores the values in m_particleDistances
    double distance;
    std::vector<double> r_i,r_j;

    for (int i=0; i<m_numberOfParticles;i++)
    {
      r_i=m_particles[i]->getPosition();
      for (int j=i; j<m_numberOfParticles;j++)
      {
        distance=0;
        if (i<j)
        {
          r_j=m_particles[j]->getPosition();
          for (int dim=0; dim < m_numberOfDimensions; dim++)
          {
            distance=distance+( r_i.at(dim)-r_j.at(dim) )*( r_i.at(dim)-r_j.at(dim) );
          }
        }
        m_particleDistances.push_back(sqrt(distance));
     }
   }
  }
  else {;}
}

void System::calculateR_jk(std::vector<Particle*> particles, int particle_i){
  // Calculates the distance between particle i and all other particles
  if (m_interaction==0) return;
  double distance;
  int index;
  std::vector<double> r_i,r_j;

  r_i=particles[particle_i]->getPosition();

  for (int j=0; j<m_numberOfParticles;j++)
  {
    distance=0;
    if (particle_i!=j)
    {
      r_j=particles[j]->getPosition();
      for (int dim=0; dim < m_numberOfDimensions; dim++)
      {
        distance+=(r_i[dim]-r_j[dim])*(r_i[dim]-r_j[dim]);
      }
      if (distance<=m_a){m_break_eval=1;}
    }

    if (particle_i <= j)
    {
      index=particle_i * (m_numberOfParticles) - (particle_i - 1) * particle_i / 2 + j - particle_i;
    }
    else
    {
      index=j * (m_numberOfParticles) - (j - 1) * j / 2 + particle_i - j;
    }
    m_particleDistances.at(index)=sqrt(distance);
  }
}

double System::getR_jk(int i, int j){
  // returns the distance between particle i and particle j
  assert(m_particleDistances.size() == (unsigned int) m_numberOfParticles*(m_numberOfParticles+1)/2.0);
  if (m_interaction>0)
  {
  int index;
  if (i <= j)
  {
    index= i * (m_numberOfParticles) - (i - 1) * i / 2 + j - i;
  }
  else
  {
    index= j * (m_numberOfParticles) - (j - 1) * j / 2 + i - j;
  }
  if ( m_particleDistances.at(index))
  return  m_particleDistances.at(index);
}
else {return 0;}
}

std::vector<double> System::getRadialDistances(){
  std::vector<double> radial_distances = std::vector<double>();
  for (int i=0; i<m_numberOfParticles;i++)
  {
    double r=0;
    for (int dim=0; dim <m_numberOfDimensions; dim++)
    {
          r+=(m_particles[i]->getPosition().at(dim) )*( m_particles[i]->getPosition().at(dim) );
    }

    radial_distances.push_back(sqrt(r));
  }
    return radial_distances;
}
