import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import matplotlib as mpl
import glob, os

import sklearn.linear_model
from matplotlib import rc
rc('text', usetex=True)
rc('font', size=14)
rc('legend', fontsize=13)
rc('text.latex', preamble=r'\usepackage{lmodern}')


from numpy import genfromtxt

l_lim=11
h_lim=20

file_number=0
for file in glob.glob("0*.dat"):
    my_data = genfromtxt(file, delimiter='\n')
    if file_number==0:
        X= my_data[-10:]

        y = np.linspace(l_lim,h_lim,10)

    else:
        X=np.concatenate((X,  my_data[-10:]))
        y=np.concatenate((y,  np.linspace(l_lim,h_lim,10)))
    file_number+=1
reg = sklearn.linear_model.LinearRegression().fit(y[:,None], X[:,None]
    )
A=reg.coef_
b=reg.intercept_

at_E=float((3-b)/A)
print("BF reg.coeff %2.e"%float(A))
print("BF at 3.0 after",round(at_E), "Iterations")
iters=np.linspace(l_lim,h_lim,10)
plt.scatter(y,X,marker=".")
plt.plot(iters,reg.predict(iters[:,None]),label="Brute Force")

X=0
y=0
file_number=0
for file in glob.glob("1*.dat"):
    my_data = genfromtxt(file, delimiter='\n')
    if file_number==0:
        X= my_data[-10:]

        y = np.linspace(l_lim,h_lim,10)

    else:
        X=np.concatenate((X,  my_data[-10:]))
        y=np.concatenate((y,  np.linspace(l_lim,h_lim,10)))
    file_number+=1
reg = sklearn.linear_model.LinearRegression().fit(y[:,None], X[:,None]
    )
A=reg.coef_
b=reg.intercept_

at_E=float((3-b)/A)
print("IS reg.coeff %2.e"%float(A))
print("IS at 3.0 after",round(at_E), "Iterations")
iters=np.linspace(l_lim,h_lim,10)
plt.scatter(y,X,marker="x")
plt.plot(iters,reg.predict(iters[:,None]),label="Importance Sampling")

plt.legend(loc="best")

plt.grid("on")
plt.xlabel(r'Training iterations ')#, fontsize=10)
plt.ylabel(r'$\bar{E}$')#,fontsize=10)
plt.savefig("Regression_BF_IS.pdf")
plt.show()
X=0
y=0
file_number=0
for file in glob.glob("2*.dat"):
    my_data = genfromtxt(file, delimiter='\n')
    if file_number==0:
        X= my_data[-10:]

        y = np.linspace(l_lim,h_lim,10)

    else:
        X=np.concatenate((X,  my_data[-10:]))
        y=np.concatenate((y,  np.linspace(l_lim,h_lim,10)))
    file_number+=1
reg = sklearn.linear_model.LinearRegression().fit(y[:,None], X[:,None]
    )
A=reg.coef_
b=reg.intercept_

at_E=float((3-b)/A)
print("GS reg.coeff %2.e"%float(A))
print("GS at 3.0 after",round(at_E), "Iterations")
iters=np.linspace(l_lim,h_lim,10)
plt.scatter(y,X,marker="x")
plt.plot(iters,reg.predict(iters[:,None]),label="Gibbs Sampling")



plt.xlabel(r'Training iterations ')#, fontsize=10)
plt.ylabel(r'$\bar{E}$')#,fontsize=10)


plt.legend(loc="best")

plt.savefig("n2_Regression_GS.pdf")
plt.show()
