#pragma once
#include "../RBM.h"
#include "Eigen/Dense"

class Hamiltonian{
public:
    Hamiltonian(class System* system);
    ~Hamiltonian();
    double computeLocalEnergy(Eigen::VectorXd X);

    double logistic(double x);
    double v_j(int j,Eigen::VectorXd X);

    double nabla(int k,Eigen::VectorXd X);
    double nabla2(int k,Eigen::VectorXd X);

private:
  class System* m_system = nullptr;
  Eigen::VectorXd m_a;
  Eigen::VectorXd m_b;
  Eigen::MatrixXd m_W;
  double m_omega;
  double m_sigma;
  double m_sigma4;
  int m_M;
  int m_N;
};
