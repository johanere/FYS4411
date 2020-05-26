#include "writefile.h"

void write_LocalEnergy(std::vector<double>& energySamples, std::string filename)
{
  /*
  Writes the sampled local energy to a file
  */

filename.append("_E.dat");
ofstream outfile;
/*
outfile.open(filename, ios::out | ios::binary | ios::trunc);
outfile.write(reinterpret_cast<const char*> (energySamples.data()),energySamples.size()*sizeof(double));
*/
outfile.open(filename);
outfile << setiosflags(ios::showpoint | ios::uppercase);

for (unsigned int i=0; i<energySamples.size();i++){
  outfile << setw(15) << setprecision(8) <<energySamples[i]<<endl;
}
outfile.close();
}

void write_distances(std::vector<double>& distances, std::string filename)
{
  /*
  Writes the sampled local energy to a file
  */
filename.append("NI10_r.dat");
ofstream outfile;
/*
outfile.open(filename, ios::out | ios::binary | ios::trunc);
outfile.write(reinterpret_cast<const char*> (energySamples.data()),energySamples.size()*sizeof(double));
*/
outfile.open(filename);
outfile << setiosflags(ios::showpoint | ios::uppercase);

for (unsigned int i=0; i<distances.size();i++){
  outfile << setw(15) << setprecision(8) <<distances[i]<<endl;
}
outfile.close();
}


void write_average_of_runs(double runs, double total_accepted_steps, double total_sigma, double total_time,std::string folder_name)
{
string filename = folder_name;
filename.append("/");
filename.append("averages.txt");
ofstream outfile;

outfile.open(filename);
outfile << setiosflags(ios::showpoint | ios::uppercase);

outfile << setw(15) << setprecision(8) <<"Accepted steps        = "<<total_accepted_steps/((double) runs)<<endl;
outfile << setw(15) << setprecision(8) <<"Time                  = "<<total_time/((double) runs)<<endl;
outfile << setw(15) << setprecision(8) << std::scientific;
outfile << setw(15) << setprecision(8) <<"Sigma (not corrected) = "<<total_sigma/((double) runs)<<endl;

outfile.close();
}
