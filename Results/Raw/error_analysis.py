# Directly copied, then modified, from lecture notes in FYS4411, which is based on the work of Marius Jonsson.

# Common imports
import os

# Where to save the figures and data files

infile = open("Energies.dat",'r')

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt, int64
from numpy.linalg import inv

def block(x):
    # preliminaries

    d = int(log2(len(x)))
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)
    M=zeros(d)
    n=zeros(d,int64)

    # estimate the auto-covariance and variances
    # for each blocking transformation

    for i in arange(0,d):
        n[i] = len(x)
        if n[i]%2 != 0:
            x=x[1::-1]
            n[i]=n[i]-1

        # estimate autocovariance of x
        gamma[i] = 1.0/float(n[i])*sum( (x[0:(n[i]-1)]-mu)*(x[1:n[i]]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::n[i]/2] + x[n[i]/2::n[i]])

    # generate the test observator M_k from the theorem
    if sum(s)==0:
        return mu, s[0]
    for j in range(0,d):
        for k in range(j,d):
            M[j]=M[j]+n[k]*( (n[k]-1) * (s[i]**2 / n[k]**2)  +gamma[i]  )**2  / s[k]**4
        #    M[j]=M[j]+(n[k]*( (n[k]-1) * (s[i]**2 / n[k]**2)  +gamma[i]  )**2) / s[k]**4

    #M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # Chi distribution for 99 percentile
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    return mu, s[k]/n[k]


x = loadtxt(infile)
(mean, var) = block(x)
std = sqrt(var)

import pandas as pd
from pandas import DataFrame

data ={'Mean':[mean], 'STDev':[std]}
frame = pd.DataFrame(data,index=['Values'])
print(frame)
