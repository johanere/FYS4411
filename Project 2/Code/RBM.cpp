#include <iostream>
#include <cmath>
#include <vector>
#include "RBM.h"
#include "system.h"
#include "Math/random.h"
using std::cout;
using std::endl;


RBM::RBM(System* system,int GD_iters,int m, int n) {
  m_GDiters=GD_iters;
  m_m=m;
  m_n=n;
  m_system=system;
}

void RBM::InitiateWeightsAndBiases() {
  double r=0;

  for (int i=0; i<m_m;i++)
  {
    r=m_system->getRandomEngine()->nextDouble();
    m_a.push_back( r ) ;
  }

  for (int j=0; j<m_n;j++)
  {
    r=m_system->getRandomEngine()->nextDouble();
    m_b.push_back(r) ;
  }

  for (int i=0; i<m_m;i++)
  {
    for (int j=0; j<m_n;j++)
    {
      r=m_system->getRandomEngine()->nextDouble();
      m_W.push_back( r );
    }
  }

} // end of InitiateWeightsAndBiases
