#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;
void write_LocalEnergy(vector<double>&,string);
void write_distances(vector<double>&, string);
void write_average_of_runs(double runs,double total_accepted_steps, double total_sigma, double total_time,string);
