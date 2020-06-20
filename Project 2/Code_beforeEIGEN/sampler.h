#pragma once
#include <vector>

class Sampler {
public:
    Sampler(class System* system,int GD_iters);
    ~Sampler();  // This is the destructor: declaration
    void setNumberOfMetropolisSteps(int steps,int steps_after_eq);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double getError()          { return m_error; }
    double getGradAlpha()          { return m_gradAlpha; }
    double getAcceptanceRate()  {return m_A;}
    //std::vector <double> getEnergySamples() { return m_EnergySamples; }

    double compute_v_j(int j);

    std::vector<double>             get_grad_a()      { return m_a; }
    std::vector<double>             get_grad_b()      { return m_b; }
    std::vector<double>             get_grad_W()      { return m_W; }



    void Gradient_in_step(double m_DeltaE);


    void Return_gradients();
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
    double  m_ElR=0;
    double  m_sumR2=0;
    double  m_gradAlpha=0;
    double  m_A=0;
    int     m_acceptedSteps=0;
    int     m_GDiters=0;
    class System* m_system = nullptr;
    //std::vector <double> m_EnergySamples;


    double m_sigma=0;

    std::vector <double>        m_psi_a_EL;
    std::vector <double>        m_psi_b_EL;
    std::vector <double>        m_psi_W_EL;

    std::vector <double>        m_psi_a;
    std::vector <double>        m_psi_b;
    std::vector <double>        m_psi_W;

    std::vector <double>        m_a;
    std::vector <double>        m_b;
    std::vector <double>        m_W;

    std::vector <double>        m_X;
};
