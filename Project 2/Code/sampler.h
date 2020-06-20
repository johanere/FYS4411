#pragma once
#include <vector>
#include "Eigen/Dense"

class Sampler {
public:
    Sampler(class System* system);
    ~Sampler();  

    void setNumberOfMetropolisSteps(int steps,int steps_after_eq,bool store_samples);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    void Update_expectations(double m_DeltaE);
    double getEnergy()          { return m_energy; }
    double getError()          { return m_error; }
    double getAcceptanceRate()  {return m_A;}
    std::vector <double> getEnergySamples() { return m_EnergySamples; }

    void blocking();

    Eigen::VectorXd                Return_gradients_a();
    Eigen::VectorXd                Return_gradients_b();
    Eigen::MatrixXd                Return_gradients_W();

    double logistic(double x);
    double v_j(int j,Eigen::VectorXd X);

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_energy2 = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergy2 = 0;
    double  m_DeltaE=0;
    double  m_variance=0;
    double  m_error=0;
    double  m_A=0;
    int     m_acceptedSteps=0;
    class System* m_system = nullptr;
    std::vector <double> m_EnergySamples;


    double mse_mean;
    double stdErr;
    double mse_stdErr;


    Eigen::VectorXd m_psi_a;
    Eigen::VectorXd m_psi_a_EL;

    Eigen::VectorXd m_psi_b;
    Eigen::VectorXd m_psi_b_EL;

    Eigen::MatrixXd m_psi_W;
    Eigen::MatrixXd m_psi_W_EL;

    double m_sigma=0;
};
