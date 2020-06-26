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

label=["BF $N=2$", "IS $N=2$","GS $N=2$","BF $N=4$","IS $N=4$","GS $N=4$"]

plt.figure()



avr=np.zeros(6)

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
    start_avr=180
    avr[file_number]=np.mean(E[start_avr:])

    if file_number>2:
        plt.plot(mcs[start_plot:],E[start_plot:],":",label="%s"%(label[file_number]))
    else:
        plt.plot(mcs[start_plot:],E[start_plot:],label="%s"%(label[file_number]),alpha=0.7)
    file_number+=1

for i in range(6):
    print("%.2e"%avr[i])
plt.legend(frameon=False, loc="best")
plt.grid("on")
plt.xlabel(r'Iterations ')#, fontsize=10)
plt.ylabel(r'$\hat{\sigma}$')#,fontsize=10)
plt.savefig("n2_error_interacting_long.pdf")
plt.show()
