import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
outdata=["Energies1.dat","Energies10.dat","Energies100.dat","Energies5000.dat"]
lab=["N=1", "N=10","N=100", "N=500",]
L=20

figure(num=None, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
plt.ticklabel_format(axis='both', style='sci')

#open for reading
for j in range(len(outdata)):
    f= open(outdata[j])
    if f.mode == 'r':
        lines = f.readlines()
    n=len(lines)

    #placeholders
    mcs=np.zeros(n)
    T=np.zeros(n)
    E=np.zeros(n)
    C=np.zeros(n)
    Mtot=np.zeros(n)
    X=np.zeros(n)
    M=np.zeros(n)
    configs=np.zeros(n)

    #extract data
    for i in range(len(lines)):
        E[i]=float(lines[i].split())
        mcs[i]=i


    avrage_E=E/np.sum(E)
    #plot
    plt.plot(mcs,avrage_E,label="%s"%(lab[j]))


N=float(hold[0])
pow=int(np.log10(N))


plt.xlabel('MC cycles', fontsize=14)
plt.ylabel('$<E>/<E>_{EQ}$',fontsize=14)
plt.grid('on')
plt.legend(loc=0)
plt.title('$<|E|>$ convergence to equilibrium for $L=%s$'%(L), fontsize=14)
plt.show()
