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

y=np.asarray([0.3, 0.94, 0.789,0.6966,0.64322,0.627312,0.532436,0.513021, 0.500806, 0.499836])
x=np.linspace(0,10,10)
z=np.zeros(10)
z=z+0.5

plt.scatter(x,y)
plt.plot(x,z)




plt.grid("on")
plt.xlabel(r'Runs')#, fontsize=10)
plt.ylabel(r'$\alpha$')#,fontsize=10)
plt.savefig("GD.pdf")
plt.show()
