
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
