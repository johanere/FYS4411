#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system,int numberOfDimensions);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> updateQForce(std::vector<class Particle*> particles,int particle)=  0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    std::vector<double> m_qForce=std::vector<double>(); 

};
