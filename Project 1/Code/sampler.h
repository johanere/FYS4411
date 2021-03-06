#pragma once
#include <vector>

class Sampler {
public:
    Sampler(class System* system,int GD_iters);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double getError()          { return m_error; }
    double getGradAlpha()          { return m_gradAlpha; }
    double getAcceptanceRate()  {return m_A;}
    std::vector <double> getEnergySamples() { return m_EnergySamples; }

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
    double   m_A=0;
    int     m_acceptedSteps=0;
    int     m_GDiters=0;
    class System* m_system = nullptr;
    std::vector <double> m_EnergySamples;


};
