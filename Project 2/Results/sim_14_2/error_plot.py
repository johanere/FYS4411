import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import matplotlib as mpl
import glob, os

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
plt.rcParams['font.family'] = 'serif'
#rc('text', usetex=True)
#rc('font', size=14)
#rc('legend', fontsize=13)
#rc('text.latex', preamble=r'\usepackage{lmodern}')

label=["Brute Force","Importance Sampling","Gibbs Sampling"]

plt.figure()





#open for reading
file_number=0
for file in glob.glob("*_error.dat"):
    print(file)
    f= open(file)
    if f.mode == 'r':
        lines = f.readlines()
    n=len(lines)

    #placeholders
    cumulative_energy=0
    cumulative_energy2=0
    mcs=np.zeros(n)
    S=np.zeros(n)
    E=np.zeros(n)
    acceptedStep=np.zeros(n)
    E_avr=np.zeros(n)

    #extract data
    step=0
    for i in range(len(lines)):
        E[i]=float(lines[i])
        step+=1
        mcs[i]=step

    #plot
    start_plot=0#int(round(step*0.0005))

    plt.plot(mcs[start_plot:],E[start_plot:],label="%s"%(label[file_number]))
    file_number+=1


plt.legend(frameon=False, loc="best")
plt.grid("on")
plt.xlabel(r'Iterations ')#, fontsize=10)
plt.ylabel(r'Mean $\hat{\sigma}$')#,fontsize=10)
plt.savefig("n2_error_interacting_long.pdf")
plt.show()
