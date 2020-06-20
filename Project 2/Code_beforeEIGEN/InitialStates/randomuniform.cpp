#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"



#include <random>

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles, double a, int interaction )  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */

    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState(a,interaction);
}

void RandomUniform::setupInitialState(double a, int interaction ) {
    //set up random number generator and distribution
    /*std::mt19937  generator (42);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    */

    std::vector<double> r_j = std::vector<double>();
    double distance=1;


    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = GeneratePosition();
        if (interaction>0 )
        {
          for (int j=0; j<i;j++)
          {
            distance=0;
              r_j=m_particles[j]->getPosition();
              for (int dim=0; dim < m_numberOfDimensions; dim++)
              {
                distance=distance+( position.at(dim)-r_j.at(dim) )*( position.at(dim)-r_j.at(dim) );
              }
            distance=sqrt(distance);

          if (distance<a)
          {cout<<"particle "<<i<<" is within a="<<a<<" of another particle ( "<<distance<<")"<<endl;
          exit(0);}
        }
      }





        m_particles.push_back(new Particle());
        m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i)->setPosition(position);
    }
}

std::vector<double>  RandomUniform::GeneratePosition(){
    double pos=0;
  std::vector<double> position = std::vector<double>();
  for (int j=0; j < m_numberOfDimensions; j++) {
      pos=(m_system->getRandomEngine())->nextDouble();
      position.push_back(pos);
  }
  return position;
}
