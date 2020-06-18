#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include "RBM.h"
#include "system.h"
#include "Math/random.h"
using std::cout;
using std::endl;


RBM::RBM(int GD_iters,int m, int n, double learningrate) {
  m_random = new Random();
  m_GDiters=GD_iters;
  m_M=m;
  m_N=n;
  m_learningrate=learningrate;
}

RBM::RBM(int GD_iters,int m, int n, double learningrate,int seed){
    m_random = new Random(seed);
    m_GDiters=GD_iters;
    m_M=m;
    m_N=n;
    m_learningrate=learningrate;
}

void RBM::WeightsAndBiases(int current_run) {
  if (current_run==1){
    InitiateWeightsAndBiases();
  }
  else{
    ;
  }
} //end of WeightsAndBiases

void RBM::Update_gradients(std::vector<double> grad_a,
std::vector<double> grad_b,std::vector<double> grad_W) {
assert(m_learningrate =! 0);
assert(m_a.size() == grad_a.size() );
assert(m_b.size() == grad_b.size() );
assert(m_W.size() == grad_W.size() );

for (int i=0; i <m_M ;i++)
{
  m_a[i]-=m_learningrate*grad_a[i];
}
for (int i=0; i <m_N ;i++)
{
  m_b[i]-=m_learningrate*grad_b[i];
}
for (int i=0; i< m_M*m_N ;i++)
{
  m_W[i]-=m_learningrate*grad_W[i];
}

} //end of WeightsAndBiases

void RBM::InitiateWeightsAndBiases() {
  double r=0;

  for (int i=0; i<m_M;i++)
  {
    r=m_random->nextDouble();
    m_a.push_back( r ) ;
  }

  for (int j=0; j<m_N;j++)
  {
    r=m_random->nextDouble();
    m_b.push_back(r) ;
  }

  for (int i=0; i<m_M;i++)
  {
    for (int j=0; j<m_N;j++)
    {
      r=m_random->nextDouble();
      m_W.push_back( r );
    }
  }

} // end of InitiateWeightsAndBiases
