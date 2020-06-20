double Sampler::compute_v_j(int j) {
  double sum=0;


  for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
  {
    sum+=m_X[i]*m_W[ (i* m_system->getRBM() ->get_m() )  +j]/(m_sigma*m_sigma);
  }
  return m_b[j]+sum ;
} // end of v_j


void Sampler::Gradient_in_step(double m_DeltaE) {
  //update X
  for (int i=0; i<m_system->getNumberOfParticles();i++)
  {
    for (int dim=0; dim <m_system->getNumberOfDimensions();dim++)
    {
      m_X[i+dim]= particles[i]->getPosition()[dim];
    }
  }

  for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
  {
    m_psi_a.push_back( (m_W[i]-m_a[i])/(m_sigma*m_sigma)  );
    m_psi_a_EL.push_back( (m_W[i]-m_a[i])/(m_sigma*m_sigma) * m_DeltaE );
  }

  for (int j=0; j<m_system->getRBM() ->get_n() ;i++)
  {
    m_psi_b.push_back( 1.0 / (1.0 + exp(-v(j)) )  );
    m_psi_b_EL.push_back( m_DeltaE / (1.0 + exp(-v(j)) ) );
  }

  for (int k=0; k<m_system->getRBM() ->get_m() ;k++)
  {
    for (int l=0; l<m_system->getRBM() ->get_n() ;l++)
    {
      m_psi_W.push_back( m_X[k] / ((1.0 + exp(-v(l)) ) * m_sigma*m_sigma  )  );
      m_psi_W_EL.push_back( m_DeltaE* m_X[k] / ((1.0 + exp(-v(l)) ) * m_sigma*m_sigma  ) );
    }
  }
} // End of Gradient_in_step




double compute_v_j(int j);
void Gradient_in_step(double m_DeltaE);

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

---------------

void Sampler::Return_gradients(){

  std::vector <double>  grad_a;
  std::vector <double>  grad_b;
  std::vector <double>  grad_W;
  //Take averages
  for (int i=0; i<m_system->getRBM() ->get_m() ;i++)
  {
    m_psi_a[i]      = m_psi_a[i]     /m_stepNumber   ;
    m_psi_a_EL[i]   = m_psi_a_EL[i]  /m_stepNumber    ;
    grad_a.push_back( 2*(m_psi_a_EL[i]-m_psi_a[i]*m_energy));
  }

  for (int j=0; j<m_system->getRBM() ->get_n() ;j++)
  {
    m_psi_b[j]      =    m_psi_b[j]     /m_stepNumber ;
    m_psi_b_EL[j]   =    m_psi_b_EL[j]  /m_stepNumber ;
    grad_b.push_back( 2*(m_psi_b_EL[j]-m_psi_b[j]*m_energy));
  }
  for (int k=0; k<m_system->getRBM() ->get_m() ;k++)
  {
    for (int l=0; l<m_system->getRBM() ->get_n() ;l++)
    {
      m_psi_W[k*    (m_system->getRBM() ->get_m())+l ]     = m_psi_W[k*    (m_system->getRBM() ->get_m())+l ]  /m_stepNumber ;
      m_psi_W_EL[k* (m_system->getRBM() ->get_m())+l ]     = m_psi_W_EL[k* (m_system->getRBM() ->get_m())+l ]  /m_stepNumber ;
      grad_W.push_back( 2*(m_psi_W_EL[k*    (m_system->getRBM() ->get_m())+l]-m_psi_W[k*    (m_system->getRBM() ->get_m())+l]*m_energy));
    }
  }
  cout<<"In sampler return grad"<<endl;
  m_system->getRBM()->Update_gradients(grad_a,grad_b,grad_W);
  cout<<"End of sampler return grad"<<endl;
}
