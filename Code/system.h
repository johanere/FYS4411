#pragma once
#include <vector>

class System {
public:
    System();
    System(int seed);
    bool metropolisStep             ();

    bool metropolishastingsStep     ();

    void runMetropolisSteps         (int numberOfMetropolisSteps,int method);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setDistances               (class Distances* distances);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class Random*                   getRandomEngine()   { return m_random; }
    class Distances*                getDistances()      { return m_distances; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }


private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;

    double                          m_stepLengthsqrt = 0.316;

    float                           oldWF=0;
    int                             particle=0;
    std::vector<double>             oldPosition=std::vector<double>();
    float                           newWF=0;
    bool                            acceptedStep;
    std::vector<double>             particle_distance=std::vector<double>();

    class Distances*                m_distances = nullptr;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    class Random*                   m_random = nullptr;

};
