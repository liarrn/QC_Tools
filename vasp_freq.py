from collections import OrderedDict
import re
import StringIO

# visulazation the vibration calculated by vasp with IBRION = 5
inFilename = "OUTCAR"
elem = OrderedDict([(96,  'O'), (98,  'H'), (146,  'Ti')])
nAtoms = 146
x = [0.0 for i in range(nAtoms)]
y = [0.0 for i in range(nAtoms)]
z = [0.0 for i in range(nAtoms)]
dx = [0.0 for i in range(nAtoms)]
dy = [0.0 for i in range(nAtoms)]
dz = [0.0 for i in range(nAtoms)]
i = 0
m = 0

with open(inFilename, 'r') as fp:
    line = fp.readline()
    while (line != ''):
        line = fp.readline()
        if (re.search('meV$', line) != None):
            m = m + 1
            fp.readline()
            for i in range(nAtoms):
                line = fp.readline()
                [x[i], y[i], z[i], dx[i], dy[i], dz[i]] = [float(coor) for coor in line.split()]
                if (i >= nAtoms):
                    print "number of atoms does not match"
                    break
            outFilename = StringIO.StringIO()
            outFilename.write("%d.xyz"%m)
            with open(outFilename.getvalue(), 'w') as fpo:
                fpo.write("%d\n"%nAtoms)
                fpo.write("%d\n"%m)
                for i in range(len(x)):
                    for j in elem.keys():
                        if i < j: 
                            fpo.write("    %s   "%elem[j])
                            break
                        elif(i >= j):
                            continue
                    fpo.write("%.3f   %.3f   %.3f   \n"%(x[i], y[i], z[i]))
                fpo.write("%d\n"%nAtoms)
                fpo.write("%d\n"%m)
                for i in range(len(x)):
                    for j in elem.keys():
                        if i < j: 
                            fpo.write("    %s   "%elem[j])
                            break
                        elif(i >= j):
                            continue
                    fpo.write("%.3f   %.3f   %.3f   \n"%(x[i]+dx[i], y[i]+dy[i], z[i]+dz[i]))