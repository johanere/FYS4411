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

label=["Brute Force","Importance Sampling","Gibbs Sampling $\sigma=1.0$","Gibbs Sampling $\sigma=0.9$"]

#["$\\Delta t= 0.005$ ","$\\Delta t= 0.05$ ","$\\Delta t= 0.5$ ","$\\Delta t= 1.0$ "]
styles=[":k","--k","-.k"]

plt.figure()

#---plt.figure(num=None, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
#----plt.ticklabel_format(axis='both', style='sci')
#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']



#open for reading


f= open("average_error")
if f.mode == 'r':
    lines = f.readlines()
n=len(lines)

#placeholders

mcs=np.asarray([0.03125,0.0625,0.125,0.25,0.5,1.0,2.0])

print(mcs)

E=np.zeros(n)

#extract data
step=0
for i in range(len(lines)):
    E[i]=float(lines[i])

    #plot


plt.plot(mcs,E)

plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.grid("on")
plt.xlabel(r' step length l ')#, fontsize=10)
plt.ylabel(r'$\hat{\sigma}$')#,fontsize=10)
plt.savefig("error.pdf")
plt.show()
