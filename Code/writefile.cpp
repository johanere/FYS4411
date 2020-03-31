#include "writefile.h"

void write_LocalEnergy(std::vector<double>& energySamples, std::string filename)
{
  /*
  Writes the sampled local energy to a file
  */
filename.append(".txt");
ofstream outfile;
/*
outfile.open(filename, ios::out | ios::binary | ios::trunc);
outfile.write(reinterpret_cast<const char*> (energySamples.data()),energySamples.size()*sizeof(double));
*/
outfile.open(filename);
outfile << setiosflags(ios::showpoint | ios::uppercase);

outfile << setw(15) << setprecision(8) <<"Energies" <<endl;
for (int i=0; i<energySamples.size();i++){
  outfile << setw(15) << setprecision(8) <<energySamples[i]<<endl;
}
outfile.close();
}
