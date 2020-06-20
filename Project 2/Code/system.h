#pragma once
#include <vector>
#include "Eigen/Dense"
#include "RBM.h"
#include "headers.h"
class System {
public:
    System();
    System(int seed);
    ~System();  
    bool brute_force_Step             ();
    bool importance_sampling_Step     ();
    void runMetropolisSteps           (int numberOfMetropolisSteps,int method);

    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWavefunction            (class Wavefunction* Wavefunction);
    void setRBM                     (class RBM* rbm,int current_run);

    class Wavefunction*             getWavefunction()   { return m_Wavefunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class Random*                   getRandomEngine()   { return m_random; }
    class RBM*                      getRBM()            { return m_rbm; }

    Eigen::VectorXd                 get_X()      { return m_X; }


    int    getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()      { return m_equilibrationFraction; }




private:
    int                             m_M = 0; // # visible nodes
    int                             m_N = 0; // # hidden nodes

    int                             m_numberOfMetropolisSteps = 0;

    int                             move_index=0;
    int                             m_interaction=0;

    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0;
    double                          m_stepLengthsqrt = 0;

    int                             m_break_eval=0;

    float                           oldWF=0;
    float                           newWF=0;

    Eigen::VectorXd                 oldState;
    Eigen::VectorXd                 m_X;

    bool                            acceptedStep;
    std::vector<double>             particle_distance=std::vector<double>();

    class Wavefunction*             m_Wavefunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class Sampler*                  m_sampler = nullptr;
    class Random*                   m_random = nullptr;
    class RBM*                      m_rbm = nullptr;


};
