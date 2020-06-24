#include <iostream>
#include "Eigen/Dense"

#include <iostream>
using std::cout;
using std::endl;

double Hola(double x){
  return 5*x;
}

int main()
{
  Eigen::VectorXd m_a;
  m_a= Eigen::VectorXd::Random(5);
  cout<<m_a<<endl;
  cout<<Hola(m_a)<<endl;
}
