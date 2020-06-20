#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

#include "RBM.h"
#include "system.h"
#include "random.h"

using std::cout;
using std::endl;


RBM::RBM(int GD_iters,int m, int n, double learningrate, double sigma,double omega) {
  m_random = new Random();
  m_GDiters=GD_iters;
  m_M=m;
  m_N=n;
  m_learningrate=learningrate;
  m_sigma=sigma;
  m_omega=omega;


}

RBM::RBM(int GD_iters,int m, int n, double learningrate,double sigma, double omega,int seed){
    m_random = new Random(seed);
    m_GDiters=GD_iters;
    m_M=m;
    m_N=n;
    m_learningrate=learningrate;
    m_sigma=sigma;
    m_omega=omega;
}

void RBM::WeightsAndBiases(int current_run) {
  if (current_run==1){
    InitiateWeightsAndBiases();
  }
  else{
    m_learningrate*=0.9;
  }
} //end of WeightsAndBiases

void RBM::InitiateWeightsAndBiases() {
     m_a = Eigen::VectorXd::Random(m_M).cwiseAbs()*0.5;
     m_b = Eigen::VectorXd::Random(m_N).cwiseAbs()*0.5;
     m_W = Eigen::MatrixXd::Random(m_M,m_N).cwiseAbs()*0.5;

}// end of InitiateWeightsAndBiases

void RBM::Update_gradients(Eigen::VectorXd grad_a,
Eigen::VectorXd grad_b, Eigen::MatrixXd grad_W) {

  assert(m_learningrate > 0);
  assert(m_a.size() == grad_a.size() );
  assert(m_b.size() == grad_b.size() );
  assert(m_W.size() == grad_W.size() );

  m_a-=m_learningrate*grad_a;
  m_b-=m_learningrate*grad_b;
  m_W-=m_learningrate*grad_W;

} //end of WeightsAndBiases
