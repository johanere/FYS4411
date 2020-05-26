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
for file in glob.glob("*_r.dat"):
    print(file)
    f= open(file)
    if f.mode == 'r':
        lines = f.readlines()
    n=len(lines)
    r=np.zeros(n)
    #extract data
    for i in range(len(lines)):
        r[i]=float(lines[i])


    M=10
    r_max=2.5#np.max(r)
    r_min= 0
    if(file_number==0):
        bins=list(np.linspace(r_min,r_max,M))
        hist_total=np.zeros(M-1)
    hist,bins = np.histogram(r,bins)
    hist_total=hist_total+hist

    file_number+=1

r_densitites=np.copy(hist_total)/(float(n)*float(file_number))
distances=np.linspace(0,r_max,len(r_densitites))

x=0
for i in range(len(hist_total)):
    r0=bins[i]
    r1=bins[i+1]
    vol=((r1+r0)/2.0)**2#4.0/3.0*np.pi*(r1**3-r0**3)
    r_densitites[i]=r_densitites[i]/vol
plt.plot(distances,r_densitites,label="%s"%(label[file_number]))


print(hist)
print(bins)


plt.legend(frameon=False, loc=2)
plt.grid("on")
plt.xlabel(r'$|r|$')#, fontsize=10)
plt.ylabel(r'$\rho(r)$')#,fontsize=10)
#plt.savefig("convergence_dt.pdf")
plt.show()
