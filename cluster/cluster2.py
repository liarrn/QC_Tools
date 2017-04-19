import numpy as np
import math
from numpy import linalg as LA
from StringIO import StringIO
import matplotlib.pyplot as plt
from matplotlib import mlab


'''
program fro guanxishi4
'''


def ReadXYZe(xyzfile, natoms):
    coor = np.array(np.zeros((1, natoms * 3)))
    energy = np.array(0.0)
    with open(xyzfile, 'r') as xyzp:
        ln = -1
        strcoor = ''
        for line in xyzp:
            ln = ln + 1
            if (ln % (natoms + 2) == 0):
                strcoor = ''
                continue
            elif (ln % (natoms + 2) == 1):
                energy = np.append(energy, float(line.split()[4]))
                continue
            elif (ln % (natoms + 2) == (natoms + 1)):
                strcoor = strcoor + line[4: -1]
                coor = np.concatenate((coor, np.array(np.matrix(strcoor).ravel().tolist())))
                # print coor.shape
                # print np.array(np.matrix(strcoor).ravel().tolist()).shape
                # coor.append(np.matrix(strcoor).getT().ravel().tolist()[0])
            else:
                strcoor = strcoor + line[4: -1] + ';'
    coor = coor[1:, :]
    energy = energy[1:]
    return coor, energy


def ReadXYZ(xyzfile, natoms):
    coor = np.array(np.zeros((1, natoms * 3)))
    # energy = np.array(0.0)
    with open(xyzfile, 'r') as xyzp:
        ln = -1
        strcoor = ''
        for line in xyzp:
            ln = ln + 1
            if (ln % (natoms + 2) == 0):
                strcoor = ''
                continue
            elif (ln % (natoms + 2) == 1):
                # energy = np.append(energy, float(line.split()[4]))
                continue
            elif (ln % (natoms + 2) == (natoms + 1)):
                strcoor = strcoor + line[4: -1]
                coor = np.concatenate((coor, np.array(np.matrix(strcoor).ravel().tolist())))
                # print coor.shape
                # print np.array(np.matrix(strcoor).ravel().tolist()).shape
                # coor.append(np.matrix(strcoor).getT().ravel().tolist()[0])
            else:
                strcoor = strcoor + line[4: -1] + ';'
    coor = coor[1:, :]
    # energy = energy[1:]
    return coor


def rephraseCoor(cluster):
    # reprase a compressed cluster coordinates to a natoms*3 representations
    natoms = cluster.size / 3
    return cluster.reshape((natoms, 3))




def CalcBondLength(coor1, coor2):
    bondlength = math.sqrt(np.power(coor1 - coor2, 2).sum())
    return bondlength


def CalcCN(atnr, cluster, atomRange=-1, rdmax=3.2):
    # the atom number, atnr, should start from 0
    # return coornr, the coordination number
    # nnatoms, array of nearest atoms
    # nnatom, the nearest atom
    # rdmax = 3.5
    # rdmax = 2.96
    coornr = 0  # the coordination number
    nnatoms = np.array([])  # the first nearest atoms
    nnatom = 0  # the nearest atom
    nnbondlength = 100.0
    natoms = cluster.shape[0]

    if type(atomRange) == int and atomRange == -1:
        atomRange = range(natoms)

    for i in atomRange:
        if i == atnr:
            continue
        bondlength = CalcBondLength(cluster[atnr, :], cluster[i, :])
        if bondlength < rdmax:
            coornr += 1
            nnatoms = np.append(nnatoms, [int(i)])
        if bondlength < nnbondlength:
            nnbondlength = bondlength
            nnatom = int(i)
    nnatoms = map(lambda s: int(s), nnatoms)
    nnatoms = np.array(nnatoms)
    coornr = int(coornr)
    return coornr, nnatoms, nnatom


def bondSet(nnatoms, cluster):
    natoms = nnatoms.size
    rdmax = 3.2
    # rdmax = 3.0
    bondlength = 0.0
    bonds = np.array([[0, 0]])
    for i in range(natoms):
        for j in range(i + 1, natoms):
            bondlength = CalcBondLength(cluster[nnatoms[i], :], cluster[nnatoms[j], :])
            if bondlength < rdmax:
                bonds = np.concatenate((bonds, [[nnatoms[i], nnatoms[j]]]))
    bonds = bonds[1:]
    return bonds


def CalcCA(atnr, cluster):
    # Calculated the average cone angle constructed by atnr's surface nearest neighbours
    # the nearest atom, and opposite bond are *** STRICTLY ON THE SURFACE *** !!!
    # return ACA, averaged cone angle
    surfatoms = getSurfAtoms(cluster)
    coornr, nnatoms, nnatom = CalcCN(atnr, cluster, surfatoms, 3.4)
    nnatoms = np.array(nnatoms)
    ACA = 0.0  # averaged cone angle
    for i in range(nnatoms.size):
        nnatom = nnatoms[i]
        tmpatoms = np.delete(nnatoms, np.argwhere(nnatoms == nnatom)[0, 0])
        bonds = bondSet(tmpatoms, cluster)
        # print bonds
        CA = 0.0
        CABond = np.array([0, 0])
        l1 = cluster[nnatom, :] - cluster[atnr, :]
        for i in range(bonds.shape[0]):
            midpoint = (cluster[bonds[i, 0], :] + cluster[bonds[i, 1], :]) / 2
            l2 = midpoint - cluster[atnr, :]
            angle = math.acos(np.dot(l1, l2) / (LA.norm(l1) * LA.norm(l2))) / math.pi * 180
            if angle > CA:
                CA = angle
                CABond = bonds[i, :]
                # print 'surface atom: %d; nnatom: %d; Bond: %d, %d; angle: %.3f'\
                #     %(int(atnr+1), int(nnatom+1), int(bonds[i, 0]+1), int(bonds[i, 1]+1), angle)
                # print l1, l2
                # print '\n'
        for i in range(tmpatoms.size):
            l2 = cluster[tmpatoms[i], :] - cluster[atnr, :]
            angle = math.acos(np.dot(l1, l2) / (LA.norm(l1) * LA.norm(l2))) / math.pi * 180
            if angle > CA:
                CA = angle
                CABond = np.array([tmpatoms[i], tmpatoms[i]])
        ACA += CA
        # print 'surface atom: %d; nnatom: %d; Bond: %d, %d; angle: %.3f'\
        #       % (int(atnr+1), int(nnatom+1), int(CABond[0]+1), int(CABond[1]+1), CA)
    ACA /= nnatoms.size
    return ACA




def getSurfAtoms(cluster):
    # return surface atoms, started with zero
    criterion = 1.6
    surfatoms = np.array([])
    natoms = cluster.shape[0]
    NormOfAccuXY = np.zeros(natoms)
    for i in range(natoms):
        coornr, nnatoms, tmp = CalcCN(i, cluster)
        accuxy = np.zeros((1, 3))
        for n in range(coornr):
            xy = cluster[nnatoms[n], :] - cluster[i, :]
            xy /= LA.norm(xy)
            accuxy += xy
        NormOfAccuXY[i] = LA.norm(accuxy)
        if NormOfAccuXY[i] > criterion:
            surfatoms = np.append(surfatoms, i)
    surfatoms = map(lambda s: int(s), surfatoms)
    surfatoms = np.array(surfatoms)
    return surfatoms
    # return NormOfAccuXY


def CalcAd(atnr, cluster):
    # Calculate the CO, O2 adsorption energy (COad, O2ad) for atnr
    # return COad, O2ad
    # atnr start from 0
    cnmax = 12.0
    COad = 0.0  # CO adsorption energy
    O2ad = 0.0  # O2 adsorption energy
    CNR = 0.0  # reduced CN
    CAR = 0.0  # reduced CA
    CAR = CalcCA(atnr, cluster)
    CAR /= 180.0
    CN, nnatoms, nnatom = CalcCN(atnr, cluster)
    for i in range(nnatoms.size):
        a, b, c = CalcCN(nnatoms[i], cluster)
        CNR += a
    CNR /= cnmax
    if CNR < 3:
        COad = 1.54 * CAR - 1.87
        O2ad = 0.93 * CAR - 1.00
    else:
        COad = 0.17 * CNR - 1.57
        O2ad = 0.2 * CNR - 1.2
    # print 'atnr: %d; reduced CA: %.2f; reduced CN: %.2f'%(atnr+1, CAR, CNR)
    return COad, O2ad


# guanxishi4
def CalcEa(atnr1, atnr2, cluster):
    # Calculate the activation energy for atnr1 and atnr2
    # return the corresponding activation energy, Ea
    # and ((atnrCO, Ead CO), (atnrO2, Ead O2))
    COad1, O2ad1 = CalcAd(atnr1, cluster)
    COad2, O2ad2 = CalcAd(atnr2, cluster)
    Ead1 = abs((COad1-O2ad2)/(COad1+O2ad2))
    Ead2 = abs((COad2-O2ad1)/(COad2+O2ad1))
    Ea1 = 1.07*Ead1-0.207
    Ea2 = 1.07*Ead2-0.207
    if Ea1 < Ea2:
        Ea = Ea1
        Ead = np.array([[atnr1, COad1], [atnr2, O2ad2]])
    else:
        Ea = Ea2
        Ead = np.array([[atnr2, COad2], [atnr1, O2ad1]])
    return Ea, Ead


def CalcLowestEa(cluster):
    # return the lowest activation energy of the cluster
    # return the corresponding bond
    COT = -0.61  # threshold for CO adsorption energy
    surfatoms = getSurfAtoms(cluster)
    bonds = bondSet(surfatoms, cluster)
    Ea, Ead = CalcEa(bonds[0, 0], bonds[0, 1], cluster)
    for i in range(1, bonds.shape[0]):
        tmp, tmpEad = CalcEa(bonds[i, 0], bonds[i, 1], cluster)
        # print bonds[i], tmp
        if (tmp < Ea) and (tmpEad[0, 1] < COT):
            Ea = tmp
            Ead = tmpEad
    return Ea, Ead


if __name__ == '__main__':

    natoms = 20

    # Calculate the CO, O2 adsorption energy
    coors = ReadXYZ('./pyramid.xyz', natoms)
    coor = rephraseCoor(coors[0, :])
    surfatoms = getSurfAtoms(coor)
    for i in range(surfatoms.size):
        COad, O2ad = CalcAd(surfatoms[i], coor)
        print 'atom id: %d; CO Ead: %.3f eV'%(int(surfatoms[i]+1), COad)

    # COad, O2ad = np.zeros(coors.shape), np.zeros(coors.shape)
    # for i in range(coors.shape[0]):
        # coor = rephraseCoor(coors[i, :])
        # natoms = coor.shape[0]
        # surfatoms = getSurfAtoms(coor)
        # for n in range(surfatoms):

