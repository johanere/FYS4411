#pragma once
#include <vector>
#include "Eigen/Dense"

class RBM {
public:
    RBM(int GD_iters, int m, int n, double learningrate, double sigma,  double omega);
    RBM(int GD_iters, int m, int n, double learningrate, double sigma,  double omega, int seed);
    void InitiateWeightsAndBiases();
    void WeightsAndBiases(int current_run);
    void Update_gradients(Eigen::VectorXd grad_a, Eigen::VectorXd grad_b, Eigen::MatrixXd grad_W);

    Eigen::VectorXd                get_a()      { return m_a; }
    Eigen::VectorXd                get_b()      { return m_b; }
    Eigen::MatrixXd                get_W()      { return m_W; }

    int get_M()                    { return m_M; }
    int get_N()                    { return m_N;  }
    double get_sigma()             { return m_sigma;  }
    double get_omega()             { return m_sigma;  }
private:
    int                         m_GDiters = 0;

    int                         m_M=0;
    int                         m_N=0;

    double                      m_learningrate=0;

    Eigen::VectorXd m_a;
    Eigen::VectorXd m_b;
    Eigen::MatrixXd m_W;

    class Random*               m_random = nullptr;

    double m_sigma=0;
    double m_omega=0;
};
