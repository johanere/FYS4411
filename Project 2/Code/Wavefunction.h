#pragma once
#include "Eigen/Dense"

class Wavefunction{
public:
    Wavefunction(class System* system);
    ~Wavefunction();
    double evaluate(Eigen::VectorXd X);
    double getQForce(Eigen::VectorXd X,int k);
    double logistic(double x);
    double v_j(int j,Eigen::VectorXd X);

  private:
      class System* m_system = nullptr;

      Eigen::VectorXd m_a;
      Eigen::VectorXd m_b;
      Eigen::MatrixXd m_W;
      double m_sigma;
      double m_sigma4;
      int m_M;
      int m_N;
};
