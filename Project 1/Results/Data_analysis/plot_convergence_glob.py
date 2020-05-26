import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import glob, os

label=["$\\Delta t= 0.005$ ","$\\Delta t= 0.05$ ","$\\Delta t= 0.5$ ","$\\Delta t= 1.0$ "]


plt.figure(num=None, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
plt.ticklabel_format(axis='both', style='sci')
#open for reading
file_number=0
for file in glob.glob("*.dat"):
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

    #extract data
    step=0
    for i in range(len(lines)):
        E[i]=float(lines[i])
        step+=1
        mcs[i]=step
        cumulative_energy+=E[i]
        cumulative_energy2+=E[i]*E[i]
        E_avr=cumulative_energy/(mcs[i])
        S[i]=cumulative_energy2/(mcs[i])-E_avr*E_avr


    #plot
    start_plot=int(round(step*0.5))

    plt.plot(np.log2(mcs[start_plot:]),np.sqrt(S/mcs)[start_plot:],label="%s"%(label[file_number]))
    file_number+=1



plt.xlabel('MC cycles after equilibrium in powers of 2', fontsize=14)
plt.ylabel('$\\sigma$',fontsize=14)
plt.grid('on')
plt.legend(loc=0)
plt.savefig("convergence_dt.pdf")
plt.show()
