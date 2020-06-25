import numpy as np

import string as str


import glob, os




table=np.zeros((9,2))


#open for reading

file ="energy_deviation.dat"

f= open(file)
if f.mode == 'r':
    lines = f.readlines()
n=len(lines)
for j in range(n):
    table[j,0]=float(lines[j])


file ="error.dat"
f= open(file)
if f.mode == 'r':
    lines = f.readlines()
n=len(lines)
for j in range(n):
    table[j,1]=float(lines[j])



eta=np.linspace(0.22,0.38,9)
print(eta)
filename="table.txt"
f = open(filename, "w+")
for i in range(8):
    f.write("$ %s $ & "%eta[i])
    for j in range(2):
        if j==1:
            f.write(" $ %.2e F $ \\\ "%table[i,j])
        else:
            f.write(" $ %.2e F &"%table[i,j])
    f.write("\n")
f.close()
