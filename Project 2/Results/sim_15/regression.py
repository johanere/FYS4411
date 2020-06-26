import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import matplotlib as mpl
import glob, os

import sklearn.linear_model
from matplotlib import rc
#rc('text', usetex=True)
#rc('font', size=14)
#rc('legend', fontsize=13)
#rc('text.latex', preamble=r'\usepackage{lmodern}')


from numpy import genfromtxt

l_lim=11
h_lim=200


my_data = genfromtxt("4_2_deltaE.dat", delimiter='\n')
X= my_data[10:]
y = np.linspace(l_lim,h_lim,190)


reg = sklearn.linear_model.LinearRegression().fit(y[:,None], X[:,None]
    )
A=reg.coef_
b=reg.intercept_

at_E=float((b)/A)
print("GS $N=4$ reg.coeff %2.e"%float(A))
print("GS $N=4$ at 3.0 after",round(at_E), "Iterations")
iters=np.linspace(l_lim,h_lim,10)
plt.scatter(y,X,marker=".")
plt.plot(iters,reg.predict(iters[:,None]),label="GS $N=4$")

X=0
y=0
file_number=0


my_data = genfromtxt("2_deltaE.dat", delimiter='\n')
X= my_data[10:]
y = np.linspace(l_lim,h_lim,190)
reg = sklearn.linear_model.LinearRegression().fit(y[:,None], X[:,None]
    )
A=reg.coef_
b=reg.intercept_

at_E=float((b)/A)
print("GS 2 reg.coeff %2.e"%float(A))
print("gs 2 at 3.0 after",round(at_E), "Iterations")
iters=np.linspace(l_lim,h_lim,10)
plt.scatter(y,X,marker="x")
plt.plot(iters,reg.predict(iters[:,None]),label="GS $N=2$")

plt.legend(loc="best")

plt.grid("on")
plt.ylim(0.0005,0.002)
plt.xlabel(r'Training iterations ')#, fontsize=10)
plt.ylabel(r'$(\Delta E)^2 $')#,fontsize=10)
plt.savefig("Regression_GS_inter.pdf")
plt.show()



plt.show()
