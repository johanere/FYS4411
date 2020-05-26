## FYS 4411 - Project 1
The project report ("FYSSTK4411project1_JN.pdf") may be found in this folder.
#### Folders
All tex files used to compile the project report may be found under "Tex".  All programs used to obtain the results described in the report is under "Code", while the results themselves may be found under "Results". "DataFiles" contains two articles, and the original data sets. 

### Abstract

This project aimed to examine using Vartiational Monte Carlo (VMC) in making an estimate, 
<img src="https://render.githubusercontent.com/render/math?math=\bar E">, of the ground state energy, 
<img src="https://render.githubusercontent.com/render/math?math=E_0">, of a trapped, hard sphere Bose gas. The error of the energy estimate was evaluated by using the sample error with "blocking", 
<img src="https://render.githubusercontent.com/render/math?math=\hat \sigma">. This method accounts for sample correlation, which was concluded to be necessary as the raw sample error was largely deceptive. Two approaches to the sampling scheme involved in the Metropolis algorithm were used, "Importance Sampling" and the "Brute Force" approach. "Importance sampling" was found to equilbrate the system faster than the "Brute Force" approach, albeit only for suitable values of 
<img src="https://render.githubusercontent.com/render/math?math=\Delta t">. Several other parameters, such as the number of particles, N and the value 
<img src="https://render.githubusercontent.com/render/math?math=\alpha"> was observed to affect the suitability of different
<img src="https://render.githubusercontent.com/render/math?math=\Delta t"> values. Increases in particle number, N, was generally observed to increase 
<img src="https://render.githubusercontent.com/render/math?math=\hat \sigma">; after 
<img src="https://render.githubusercontent.com/render/math?math=2^{18}"> Monte Carlo cycles, 
<img src="https://render.githubusercontent.com/render/math?math=\hat \sigma \sim 10^3"> (N=10) increased to 
<img src="https://render.githubusercontent.com/render/math?math=\sim 10^2"> (N=50) and 
<img src="https://render.githubusercontent.com/render/math?math=\sim 10^{-1}"> (N=100), both with and without the Jastrow factor, which accounts for particle-particle correlation.  Incorporating the Jastrow factor was concluded to have an increasing effect on the average estimated energy per particle 
<img src="https://render.githubusercontent.com/render/math?math=\bar E /N"> as N increased, from 
<img src="https://render.githubusercontent.com/render/math?math=\bar E/N\approx2.45\approx"> constant for all values of N, to 4.55 (N=50) and 2.66 (N=100). This impact dependence was further verified by qualitative comparison of one-body densities. 
