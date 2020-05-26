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

label=["Interacting","Non-interacting"]

M=20
r_max=4.0#np.max(r)
r_min=0.0
plt.figure()
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




    if(file_number==0):
        bins=list(np.linspace(r_min,r_max,M))
        hist_total=np.zeros(M-1)
    hist,bins = np.histogram(r,bins)
    hist_total=hist_total+hist

    file_number+=1
    print(hist_total)

r_densitites=np.copy(hist_total)/(float(n)*float(file_number))
distances=np.linspace(0,r_max,len(r_densitites))

x=0
for i in range(len(hist_total)):
    r0=bins[i]
    r1=bins[i+1]
    vol=((r1+r0)/2.0)**2#4.0/3.0*np.pi*(r1**3-r0**3)
    r_densitites[i]=r_densitites[i]/vol
plt.plot(distances,r_densitites,label="%s"%(label[0]))



file_number=0
for file in glob.glob("*NI_r.dat"):
    print(file)
    f= open(file)
    if f.mode == 'r':
        lines = f.readlines()
    n=len(lines)
    r=np.zeros(n)
    #extract data
    for i in range(len(lines)):
        r[i]=float(lines[i])





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
plt.plot(distances,r_densitites,label="%s"%(label[1]))


plt.legend(frameon=False, loc=0)
plt.grid("on")
plt.xlabel(r'$|r|$')#, fontsize=10)
plt.ylabel(r'$\rho(r)$')#,fontsize=10)
#plt.savefig("convergence_dt.pdf")
plt.show()
