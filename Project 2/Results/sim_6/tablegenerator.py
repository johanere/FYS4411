import numpy as np

import string as str


import glob, os




table=np.zeros((10,6))


#open for reading
file_number=0
for file in glob.glob("*error.dat"):
    print(file)
    f= open(file)
    if f.mode == 'r':
        lines = f.readlines()
    n=len(lines)

    for j in range(n):
        table[j,file_number*2]=float(lines[j])
    file_number+=1

file_number=0
for file in glob.glob("*energy_deviation.dat"):
    print(file)
    f= open(file)
    if f.mode == 'r':
        lines = f.readlines()
    n=len(lines)

    for j in range(n):
        table[j,file_number*2+1]=float(lines[j])
    file_number+=1

print(table)

eta=np.linspace(0.26,0.44,10)
print(eta)
filename="table.txt"
f = open(filename, "w+")
for i in range(10):
    f.write("$ %s $ & "%eta[i])
    for j in range(6):
        if j==5:
            f.write(" $ %.2e F $ \\ "%table[i,j])
        else:
            f.write(" $ %.2e F &"%table[i,j])
    f.write("\n")
f.close()
