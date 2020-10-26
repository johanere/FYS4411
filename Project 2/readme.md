## FYS4411 - Project 2
The project report ("Project_paper.pdf") can be found in this folder.
#### Folders
All tex files used to compile the project report may be found under "Tex".  All programs used to obtain the results described in the report is under "Code". Python code used to analyse output from the C++ VMC engine may be found under "data_analysis" in "Results".

### Abstract
In this project, which is an extension paper to of a previous project paper on Variational Monte Carlo (VMC) methods,
 training a restricted Boltzmann machine (RBM), with N hidden nodes, and M visible nodes, using VMC methods is the main focus.
 Two Metropolis Hastings methods, the brute force approach (BF) and Importance Sampling (IS) are tested and compared to Gibbs Sampling (GS).
 The methods are evaluated using the estimated energy 
<img src="https://render.githubusercontent.com/render/math?math=\bar E">
for both a system of one electron in a one-dimensional harmonic oscillator (HO) 
and an interacting system of two electrons in a two-dimensional HO. BF (using learning rate 
<img src="https://render.githubusercontent.com/render/math?math=\eta=0.38">
) and IS (
<img src="https://render.githubusercontent.com/render/math?math=\eta=0.34">
) is 
shown to results in similar optimization trajectories for the RBM; a short initial phase (
<img src="https://render.githubusercontent.com/render/math?math=\sim 10">
 iterations) of rapid optimization, 
followed by (
<img src="https://render.githubusercontent.com/render/math?math=\sim 50">
to 
<img src="https://render.githubusercontent.com/render/math?math=\sim 100">
 training iterations for N=2 and N=4 respectively) of a slow and close-to linear optimization
trajectory (regression coefficient 
<img src="https://render.githubusercontent.com/render/math?math=\sim 10^{-3}">
on iterations 10 to 20). Lastly, a phase of exponentially decaying training rates, 
which converged both for N=2 and N=4, was observed (up to 200 iterations). RBM optimization using Gibbs Sampling (GS) was conclusively
 shown to follow different optimization trajectories than BF and IS. GS (
<img src="https://render.githubusercontent.com/render/math?math=\eta=0.26">)
only exhibited two different training behaviors; rapid 
initial training (<10  iterations), and a phase of close to no optimization with constant 
<img src="https://render.githubusercontent.com/render/math?math=\bar E">
(regression coefficient 
<img src="https://render.githubusercontent.com/render/math?math=\sim 10^{-4}">
persistently up to 200 training iterations). The value of the standard deviation of the position distribution from GS, 
<img src="https://render.githubusercontent.com/render/math?math=\sigma">
, was shown to directly control 
<img src="https://render.githubusercontent.com/render/math?math=\bar E ">
 in the "no training-phase" of GS.



