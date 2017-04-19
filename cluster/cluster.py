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


def ball(n=1):
    # return an (n*3) array, with each point within a ball with radius of 1.0
    coor = np.zeros((n, 3))
    for i in range(n):
        while (1):
            tmpCoor = np.random.rand(1, 3) * 2 - 1.0
            dist = np.sqrt(LA.norm(tmpCoor, 2))
            if (dist < 1.0):
                break
        coor[i, :] = tmpCoor
    return coor


def tooclose(atom, cluster, criteria, natoms):
    # return 1 if atom is too close to cluster, return 0 if new atom is not too close
    result = 0
    for i in range(natoms):
        dist = np.sqrt(LA.norm(cluster[i, :] - atom, 2))
        if dist < criteria:
            result = 1
            break
    return result


def toofar(atom, cluster, criteria, natoms):
    # return 1 if atom is too far away from cluster, return 0 if new atom is not too far away
    result = 1
    for i in range(natoms):
        dist = np.sqrt(LA.norm(cluster[i, :] - atom, 2))
        if dist < criteria:
            result = 0
            break
    return result


def printxyz(coor, filename, permission='a', description=''):
    natoms = coor.shape[0]
    with open(filename, permission) as fp:
        fp.write('%d\n%s\n' % (natoms, description))
        for i in range(natoms):
            fp.write('\tAu\t%.3f\t%.3f\t%.3f\n' % (coor[i, 0], coor[i, 1], coor[i, 2]))


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


'''
def CalcCA(atnr, cluster):
    # find the nearest neighbour of atnr atom, and construct a cone angle based on nn, antr
    # the nearest atom, and opposite bond are *** STRICTLY ON THE SURFACE *** !!!
    # return CA, cone angle
    # the nearest atom
    # the opposite bond
    surfatoms = getSurfAtoms(cluster)
    coornr, nnatoms, nnatom = CalcCN(atnr, cluster, surfatoms)
    nnatoms = np.array(nnatoms)
    nnatoms = np.delete(nnatoms, np.argwhere(nnatoms == nnatom)[0, 0])
    bonds = bondSet(nnatoms, cluster)
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
    for i in range(nnatoms.size):
        l2 = cluster[nnatoms[i], :] - cluster[atnr, :]
        angle = math.acos(np.dot(l1, l2) / (LA.norm(l1) * LA.norm(l2))) / math.pi * 180
        if angle > CA:
            CA = angle
            CABond = np.array([nnatoms[i], nnatoms[i]])
    return CA, nnatom, CABond
'''


'''
def getSurfAtoms(cluster):
    # return surface atoms, started with zero
    nn = 8
    # if the coordination number is less or equal than nn, this atom will be considered as surface atom
    result = np.array([])
    natoms = cluster.shape[0]
    cn = np.zeros(natoms)
    a, b = 0, 0
    for i in range(natoms):
        cn[i], a, b = CalcCN(i, cluster)
        if cn[i] <= nn:
                result = np.append(result, i)
    result = map(lambda s: int(s), result)
    return result
 '''


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
    if CNR < 4:
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


'''
#new guanxishi
def CalcEa(atnr1, atnr2, cluster):
    # Calculate the activation energy for atnr1 and atnr2
    # return activation energy
    # and the corresponding activation site, [CO site, O2 site]
    # and reaction energy profile [Ecoad, ETS1, EOOCO, ETS2]
    COad, O2ad = np.zeros(2), np.zeros(2)
    Ecoad, ETS1, EOOCO, ETS2 = np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)
    Ea = np.zeros(2)
    COad[0], O2ad[1] = CalcAd(atnr1, cluster)
    COad[1], O2ad[0] = CalcAd(atnr2, cluster)
    Ecoad = 1.51 * (COad + O2ad) + 0.38  # energy level of coadsorption
    ETS1 = 1.39 * (COad + O2ad) + 0.37  # energy level of first TS
    EOOCO = 2.25 * (COad + O2ad) + 0.95  # energy level of OOCO
    ETS2 = 1.87 * (COad + O2ad) + 0.76  # energy level of second TS
    Ea = ETS2 - Ecoad
    if Ea[0] <= Ea[1]:
        return Ea[0], np.array[atnr1, atnr2], np.array[Ecoad[0], ETS1[0], EOOCO[0], ETS2[0]]
    else:
        return Ea[1], np.array[atnr2, atnr1], np.array[Ecoad[1], ETS1[1], EOOCO[1], ETS2[1]]


def CalcLowestEa(cluster):
    # return the lowest activation energy of the cluster
    # return the corresponding bond
    surfatoms = getSurfAtoms(cluster)
    bonds = bondSet(surfatoms, cluster)
    Ea, ActiveSite, EProfile = CalcEa(bonds[0, 0], bonds[0, 1], cluster)
    for i in range(1, bonds.shape[0]):
        tmpEa, tmpAS, tmpEP = CalcEa(bonds[i, 0], bonds[i, 1], cluster)
        # print bonds[i], tmp
        if tmpEa < Ea:
            Ea = tmpEa
            ActiveSite = tmpAS
            EProfile = tmpEP
    return Ea, ActiveSite, EProfile
'''


def guptaPotential(coor, gradToggle=0):
    r0, A, delta, p, q = 2.8843, 0.2061, 1.79, 10.229, 4.036
    # V = A*SUM(j=1..n)(SUM(i!=j)(exp(-p*(rij/r0-1)))) - delta*SUM(j=1..n)(sqrt(SUM(i!=j)(exp(-2q*(rij/r0-1)))))
    grad = np.zeros(coor.shape)
    natoms = coor.shape[0]
    gupta, gupta_pair, gupta_nonlocal = 0.0, 0.0, 0.0
    dist = np.zeros((natoms, natoms))  # distance matrix
    rho = np.zeros(natoms)  # density array
    rhotmp = 0.0

    for i in range(natoms - 1):
        for j in range(i + 1, natoms):
            dist[i, j] = CalcBondLength(coor[i, :], coor[j, :])
            dist[j, i] = dist[i, j]
    for i in range(natoms - 1):
        for j in range(i + 1, natoms):
            rij = dist[i, j] / r0  # reduced distance
            gupta_pair += math.exp(p * (1 - rij))
            rho_tmp = math.exp(2 * q * (1 - rij))
            rho[i] += rho_tmp
            rho[j] += rho_tmp
    gupta_pair *= 2 * A
    # Note that in orignal Cleri & Rosato paper, they define
    # 2-body term as a sum over all atom pairs, but it is much quicker to
    # evaluate sum over just j>i and double the constant A
    for i in range(natoms):
        rho[i] = math.sqrt(rho[i])
        gupta_nonlocal -= delta * rho[i]
    gupta = gupta_nonlocal + gupta_nonlocal

    # Calculate gradient
    if gradToggle:
        for i in range(natoms - 1):
            for j in range(i + 1, natoms):
                rij = dist[i, j] / r0
                vtmp = (q * delta * (1.0 / rho[i] + 1.0 / rho[j]) * math.exp(2 * q * (1 - rij)) -
                        2 * A * p * math.exp(p * (1 - rij))) / (r0 * rij)
                dx = coor[i, :] - coor[j, :]
                grad[i, :] += vtmp * dx
                grad[j, :] -= vtmp * dx
    return gupta, grad


# def guptaDerivative(coor):
#     grad = np.zeros(coor.shape)
#     r0, A, delta, p, q = 2.884, 0.2061, 1.79, 10.229, 4.036
#     natoms = coor.shape[0]
#     dist = np.zeros((natoms, natoms))
#
#     for i in range(natoms):
#         for j in range(i + 1, natoms):
#             dist[i, j] = CalcBondLength(coor[i, :], coor[j, :])
#             dist[j, i] = dist[i, j]
#
#     for i in range(natoms):
#         for j in range(natoms):
#             if i == j:
#                 continue
#             rij = dist[i, j]
#             grad[i, 0] += -p / r0 * math.exp(-p * (rij / r0 - 1)) * coor[i, 0] / rij * A
#             grad[i, 1] += -p / r0 * math.exp(-p * (rij / r0 - 1)) * coor[i, 1] / rij * A
#             grad[i, 2] += -p / r0 * math.exp(-p * (rij / r0 - 1)) * coor[i, 2] / rij * A
#
#     for i in range(natoms):
#         for j in range(natoms):
#             if i == j:
#                 continue
#             rij = dist[i, j]
#             grad[i, 0] += delta * q / r0 * math.exp(-2 * q * (rij / r0 - 1)) * coor[i, 0] / rij
#             grad[i, 1] += delta * q / r0 * math.exp(-2 * q * (rij / r0 - 1)) * coor[i, 1] / rij
#             grad[i, 2] += delta * q / r0 * math.exp(-2 * q * (rij / r0 - 1)) * coor[i, 2] / rij
#
#     return grad


def steepestDescent(coor):
    encut = 0.002
    alpha = 0.002

    natoms = coor.shape[0]
    energy, energy0 = 0.0, 0.2
    c = coor
    i = 0
    while (energy0 - energy) >= encut:
        energy, grad = guptaPotential(coor, 1)
        energy0 = energy
        print 'iternation: ' + str(i) + '  energy: ' + str(energy)
        printxyz(coor, 'opt.xyz', 'a')
        coor -= alpha * grad
        energy, tmp = guptaPotential(coor)
        i += 1
    print 'iternation: ' + str(i) + '  energy: ' + str(energy)
    printxyz(coor, 'opt.xyz', 'a')
    return coor, energy


def randLJ(natoms):
    rdmin = 2.7
    rdmax = 2.9
    # length = (natoms*4.2)**(1/3)
    radius = 2.0 + 4.0 * natoms ** (1.0 / 3)

    dis = 0.0
    coor = np.zeros((natoms, 3))
    coor[0, :] = ball() * radius
    fail = 0
    for i in range(1, natoms):
        if fail == 1:
            break
        print 'generating the %dth atom' % i
        notOK = 1
        iter = 0
        while (notOK == 1):
            iter += 1
            if iter > 100:
                fail = 1
                break
            tmpCoor = ball() * radius
            # notOK = 0 if i == 1 else (tooclose(tmpCoor, coor, rdmin, i - 1) \
            # or toofar(tmpCoor, coor, rdmax, i - 1))
            notOK = 0 if i == 1 else tooclose(tmpCoor, coor, rdmin, i - 1)
            coor[i, :] = tmpCoor
            printxyz(coor, 'GEOMETRY.xyz', 'a')
            # print 'tooclose %d'%tooclose(tmpCoor, coor, rdmin, i-1)
            # print coor
            # print tmpCoor
        coor[i, :] = tmpCoor
    if fail != 1:
        printxyz(coor, 'GEOMETRY.xyz', 'a')
    printxyz(coor, 'GEOMETRY.xyz', 'a')


if __name__ == '__main__':
    # for i in range(100):
    # randLJ(20)
    # print i

    # cd D:\research\Au55\rand\genClusters-py
    # from rand import *

    a, b = 0, 0
    natoms = 55

    '''
    coors, energy = ReadXYZe('Au_cluster/T3.xyz', natoms)
    lowestEnergy = energy[0]
    highestEnergy = lowestEnergy + 0.5
    coors = coors[energy < highestEnergy, :]
    energy = energy[energy < highestEnergy]

    nrconf = coors.shape[0]
    Ea = np.zeros(nrconf)  # the lowest activation energy for each conf
    Ead = np.zeros((nrconf, 2, 2))  # the corresponding atoms and their adsorption energy for CO and O2
    for i in range(nrconf):
        print 'Calculating the %dth configuration'%i
        coor = rephraseCoor(coors[i, :])
        Ea[i], Ead[i, :] = CalcLowestEa(coor)
    index = Ea.argsort()
    for i in range(nrconf):
        print 'Writing the %dth configuration'%i
        coor = rephraseCoor(coors[index[i], :])
        des = 'lowest Ea: %.2f eV\tCO adsorption site: %d: Ead: %.2f eV\tO2 adsorption site: %d: Ead: %.2f eV\tEnergy: %.2f eV'\
              %(Ea[index[i]], Ead[index[i], 0, 0]+1, Ead[index[i], 0, 1], Ead[index[i], 1, 0]+1, Ead[index[i], 1, 1], energy[index[i]])
        printxyz(coor, './Au_cluster/active.xyz', description=des)
    '''



    # output the bond length distribution
    # bonds = np.array([])
    # for i in range(coors.shape[0]):
    #     coor = rephraseCoor(coors[i, :])
    #     for m in range(natoms):
    #         for n in range(m+1, natoms):
    #             bonds = np.append(bonds, CalcBondLength(coor[m, :], coor[n, :]))
    # np.savetxt('bonds.txt', bonds)
    # print bonds.shape
    # plt.xticks(np.arange(2, 12, 0.5))
    # plt.hist(bonds, bins=1000)

    # output the coordination number distribution
    # rdmax = 3.5
    # coornr = 0
    # coornrs = np.array([])
    # for i in range(coors.shape[0]):
    #     coor = rephraseCoor(coors[i, :])
    #     natoms = coor.shape[0]
    #     for n in range(natoms):
    #         coornr, a, b = CalcCN(n, coor)
    #         coornrs = np.append(coornrs, coornr)
    # np.savetxt('coornrs.txt', coornrs)
    # plt.hist(coornrs, bins=100)
    # print coornrs.shape

    # output the norm of accumulated xy distribution
    # NormOfAccuXY = np.array([])
    # for i in range(coors.shape[0]):
    #     coor = rephraseCoor(coors[i, :])
    #     natoms = coor.shape[0]
    #     NormOfAccuXY = np.append(NormOfAccuXY, getSurfAtoms(coor))
    # np.savetxt('NormOfAccuXY.txt', NormOfAccuXY)

    # output the surface atoms
    # for i in range(coors.shape[0]):
    #     coor = rephraseCoor(coors[i, :])
    #     natoms = coor.shape[0]
    #     surfatoms = getSurfAtoms(coor)
    #     with open('surf.xyz', 'a') as fp:
    #         fp.write('%d\n'%(int(natoms)))
    #         fp.write('%d\n'%(int(i)))
    #         for n in range(natoms):
    #             if n in surfatoms:
    #                 fp.write('\tAu\t%.3f\t%.3f\t%.3f\n'%(coor[n, 0], coor[n, 1], coor[n, 2]))
    #             else:
    #                 fp.write('\tCu\t%.3f\t%.3f\t%.3f\n'%(coor[n, 0], coor[n, 1], coor[n, 2]))

    # Calculation the cone angle
    # coor = rephraseCoor(coors[0, :])
    # surfatoms = getSurfAtoms(coor)
    # for i in range(surfatoms.shape[0]):
    #     CA, nnatom, CABond = CalcCA(surfatoms[i], coor)
    #     print 'surface atom: %d; nnatom: %d; CA Bond: %d, %d; CA: %.3f\n'\
    #         %(int(surfatoms[i]+1), int(nnatom+1), int(CABond[0]+1), int(CABond[1]+1), CA)


    # coor = rephraseCoor(coors[0, :])
    # surfatoms = getSurfAtoms(coor)
    # print surfatoms
    # for i in range(surfatoms.shape[0]):
    #     ACA = CalcCA(surfatoms[i], coor)
    #     print 'surface atom: %d; averaged CA: %.3f' % (int(surfatoms[i]+1), ACA)
    # CA, nnatom, CABond = CalcCA(1, coor)

    # Calculation the coordination number
    # coor = rephraseCoor(coors[0, :])
    # surfatoms = getSurfAtoms(coor)
    # for i in range(surfatoms.shape[0]):
    #     CN, nnatoms, nnatom = CalcCN(surfatoms[i], coor)
    #     print 'surface atom: %d; CN: %d; nnatoms: '\
    #         %(int(surfatoms[i]+1), int(CN))
    #     print nnatoms
    #     print '\n'

    # Calculate the CO, O2 adsorption energy
    # coor = rephraseCoor(coors[0, :])
    # surfatoms = getSurfAtoms(coor)
    # for i in range(surfatoms.size):
    #     COad, O2ad = CalcAd(surfatoms[i], coor)
    #     print 'atom id: %d; CO Ead: %.3f eV; O2 Ead: %.3f eV    '%(int(surfatoms[i]+1), COad, O2ad)

    # COad, O2ad = np.zeros(coors.shape), np.zeros(coors.shape)
    # for i in range(coors.shape[0]):
    #     coor = rephraseCoor(coors[i, :])
    #     natoms = coor.shape[0]
    #     surfatoms = getSurfAtoms(coor)
    #     for n in range(surfatoms):


    # Calculate the lowest Ea for a cluster
    coors = ReadXYZ('./Au_cluster/T3.xyz', natoms)
    coor = rephraseCoor(coors)
    Ea, Ead = CalcLowestEa(coor)
    print 'lowest Ea: %.2f eV\nCO adsorption site: %d: Ead: %.2f eV\nO2 adsorption site: %d: Ead: %.2f eV'\
          %(Ea, Ead[0, 0]+1, Ead[0, 1], Ead[1, 0]+1, Ead[1, 1])

    # print coor
    # potential, grad = guptaPotential(coor, 1)
    # print potential
    # print grad
    # grad = guptaDerivative(coor)
    # print 'Gupta potential: ' + str(potential)
    # print 'Gupta grad: '
    # print grad

    # coor, energy = steepestDescent(coor)
    # print coor
    # printxyz(coor, 'opt.xyz', 'w')

    # (a, b, c) = CalcCN(28, coor)
    # print a
    # print b
    # print 'nearest neighbour: ' + str(c)
    # (a, b) = CalcCA(28, coor)
    # print a
    # print 'opposite bond: ' + str(b)

    # print c
