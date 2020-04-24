#include "distances.h"
#include <cassert>
#include <math.h>
#include "../particle.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

#include <iostream>
using namespace std;

//Vector containing all distances between particles
Distances::Distances(class System* system,
  int numberOfDimensions,
  int numberOfParticles) {
assert(numberOfDimensions > 0 && numberOfParticles > 0);
m_system = system;
m_numberOfDimensions = numberOfDimensions;
m_numberOfParticles = numberOfParticles;
}

void Distances::InitiateR(std::vector<Particle*> particles){
// Calculates the distance between all particles and stores the values in m_particleDistances
  double distance;
  std::vector<double> r_i,r_j;

  for (int i=0; i<m_numberOfParticles-1;i++)
  {
    r_i=particles[i]->getPosition();
    for (int j=i; j<m_numberOfParticles;j++)
    {
      distance=0;
      if (i<j){
        r_j=particles[j]->getPosition();
        for (int dim=0; dim < m_numberOfDimensions; dim++)
        {
          distance=distance+( r_i.at(dim)-r_j.at(dim) )*( r_i.at(dim)-r_j.at(dim) );
        }
      }
      m_particleDistances.push_back(sqrt(distance));
   }
  }

}

void Distances::calculateR_jk(std::vector<Particle*> particles, int particle_i){
  // Calculates the distance between particle i and all other particles


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

double Distances::getR_jk(int i, int j){
  // returns the distance between particle i and particle j
  int index;
  if (i <= j)
  {
    index= i * (m_numberOfParticles) - (i - 1) * i / 2 + j - i;
  }
  else
  {
    index= j * (m_numberOfParticles) - (j - 1) * j / 2 + i - j;
  }
  return m_particleDistances[index];
}
