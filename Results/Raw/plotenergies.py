import numpy as np
import matplotlib.pyplot as plt


outdata="output_test.txt"


#open for reading
f= open(outdata)
if f.mode == 'r':
    lines = f.readlines()
n=len(lines)

Energies=[]
#extract Energy bins
for j in range(len(lines)):
    hold=lines[j].split()
    for i in range(len(hold)):
        Energies.append(float(hold[i]))
Energies=np.asarray(Energies)

print(Energies)
