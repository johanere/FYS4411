#pragma once
#include <vector>

class RBM {
public:
    RBM(int GD_iters, int m, int n, double learningrate, double sigma);
    RBM(int GD_iters, int m, int n, double learningrate, double sigma, int seed);
    void InitiateWeightsAndBiases();
    void WeightsAndBiases(int current_run);
    void Update_gradients(std::vector<double> grad_a, std::vector<double> grad_b,std::vector<double> grad_makeW);

    std::vector<double>             get_a()      { return m_a; }
    std::vector<double>             get_b()      { return m_b; }
    std::vector<double>             get_W()      { return m_W; }

    int get_m()             { return m_M; }
    int get_n()             { return m_N; }
    double get_sigma()             { return m_sigma; }
private:
    int                         m_GDiters = 0;
    int                         m_M=0;
    int                         m_N=0;
    double                      m_learningrate=0;

    double                      m_sigma=0;

    std::vector <double>        m_a;
    std::vector <double>        m_b;
    std::vector <double>        m_W;



    class Random*               m_random = nullptr;

};
