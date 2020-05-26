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

plt.figure()

y_low=np.asarray([0.5,0.499665,0.498037,0.496836,0.495535,0.497029,0.497260,0.497318,0.496750,0.497055])
y_mid=np.asarray([0.5,0.497909,0.497408,0.496520,0.496099,0.496428,0.496589,0.496339,0.496673,0.496342])
y_high=np.asarray([0.5,0.497837,0.497072,0.496357,0.496347,0.496144,0.496117,0.496113,0.496005,0.495985,])
x=np.linspace(0,10,10)
z=np.zeros(10)
z=z+0.5

plt.scatter(x,y_low)
plt.plot(x,y_low,label="$2^{15}$ MC cycles")
plt.scatter(x,y_mid)
plt.plot(x,y_mid,label="$2^{18}$ MC cycles")
plt.scatter(x,y_high)
plt.plot(x,y_high,label="$2^{21}$ MC cycles")


plt.grid("on")
plt.xlabel(r'Runs')#, fontsize=10)
plt.ylim((0.495,0.50))
plt.ylabel(r'$\alpha$')#,fontsize=10)
plt.legend(loc=0)
plt.savefig("GD_interacting.pdf")
plt.show()
