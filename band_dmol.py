import numpy as np
import matplotlib.pyplot as plt

filename = 'sp5.txt'
numsep = 113
xsep = []
ysep = np.array([np.zeros(numsep)])
with open(filename, 'r') as fp:
    for i in range(numsep):
        line = fp.readline()
        xsep += [float(line.split()[0])]
        
with open(filename, 'r') as fp:
    linenum = 0
    band = []
    for line in fp:
        linenum += 1
        if linenum % (numsep + 1) == 0:
            ysep = np.append(ysep, [band], axis = 0)
            # print size(band)
            band = []
            continue
        band += [float(line.split()[1])]

ysep = ysep[1:, :]
for i in range(ysep.shape[0]):
    plt.plot(xsep, ysep[i, :])
    
filtered = np.array(filter(lambda l: l[0] >-0.5 and l[0] < 3.0, ysep))
plt.figure(num=None, figsize=(12, 12), dpi=1200, facecolor='w', edgecolor='k')
plt.xlabel('k', fontsize=20)
plt.ylabel('E(k) (eV)', fontsize=20)
for i in range(filtered.shape[0]):
    plt.plot(xsep, filtered[i, :])
plt.savefig('C:/users/yat/desktop/band.jpg', dpi=1200, facecolor='w', edgecolor='k')
