#pragma once
#include <vector>

class RBM {
public:
    RBM(class System* system,int GD_iters, int m, int n);
    void InitiateWeightsAndBiases();

    std::vector<double>             get_a()      { return m_a; }
    std::vector<double>             get_b()      { return m_b; }
    std::vector<double>             get_W()      { return m_W; }

    int get_m()             { return m_m; }
    int get_n()             { return m_n; }

private:
    int                         m_GDiters = 0;
    int                         m_m=0;
    int                         m_n=0;

    class System*               m_system = nullptr;

    std::vector <double>        m_a;
    std::vector <double>        m_b;
    std::vector <double>        m_W;


};
