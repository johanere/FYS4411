#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Eigen/Dense"

using namespace std;
void write_LocalEnergy(vector<double> energySamples,string filename,string foldername);
void write_distances(vector<double>&, string);
void write_average_of_runs(double runs,double total_accepted_steps, double total_sigma, double total_time,string);
void write_vector(Eigen::VectorXd vector, std::string filename, std::string foldername);
