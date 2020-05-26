import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import matplotlib as mpl
import glob, os


from matplotlib import rc
rc('text', usetex=True)
rc('font', size=14)
rc('legend', fontsize=13)
rc('text.latex', preamble=r'\usepackage{lmodern}')

label=["Brute Force","$\\Delta t=0.01$","$\\Delta t=0.005$"]
#["$\\Delta t= 0.005$ ","$\\Delta t= 0.05$ ","$\\Delta t= 0.5$ ","$\\Delta t= 1.0$ "]


plt.figure()

#---plt.figure(num=None, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
#----plt.ticklabel_format(axis='both', style='sci')
#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']



#open for reading
file_number=0
for file in glob.glob("*.csv"):
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
        cumulative_energy+=E[i]
        cumulative_energy2+=E[i]*E[i]
        E_avr[i]=cumulative_energy/step



    #plot
    start_plot=0#int(round(step*0.0005))
    E_bar=np.sum(E_avr)/n
    E_bar=np.zeros(n)+E_bar
    E_diff=np.abs(E_avr-E_bar)
    plt.plot(np.log2(mcs[start_plot:]),E_diff[start_plot:],label="%s"%(label[file_number]))
    file_number+=1



plt.legend(frameon=False, loc=0)
plt.grid("on")
plt.xlabel(r'MC cycles')#, fontsize=10)
plt.ylabel(r'$|\bar{E}-\langle \bar{E}\rangle |$')#,fontsize=10)
plt.savefig("convergence_dt.pdf")
plt.show()
