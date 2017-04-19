import numpy as np


def toofar(atom ,coor, criteria, natoms):
    #  return True if atom is too far away from cluster
    #  return False if new atom is not too far away
    is_toofar = True
    for i in range(natoms):
        dist = np.sqrt(np.sum((coor[i, :] - atom)**2))
        if dist <= criteria:
            is_toofar = False
            break
    return is_toofar

def tooclose(atom ,coor, criteria, natoms):
    #  return True if atom is too close away from cluster
    #  return False if new atom is not too close
    is_tooclose = False
    for i in range(natoms):
        dist = np.sqrt(np.sum((coor[i, :] - atom)**2))
        if dist <= criteria:
            is_tooclose = True
            break
    return is_tooclose

def ball():
    # return coor in a uniform ball
    coor = np.random.rand(3) * 2 - 1
    dist = np.sqrt(np.sum(coor**2))
    while dist > 1.0:
        coor = np.random.rand(3) * 2 - 1
        dist = np.sqrt(np.sum(coor**2))
    return coor

def randcluster(natoms, rdmin, rdmax):
    # return coor of rand cluster
    length = (natoms*5.0)**(1/3.0)

    coor = np.zeros([natoms, 3])
    coor[0, : ] = ball()*length
    is_fail = False
    for i in range(1, natoms):
        if is_fail == True:
            break
        is_satisfied = False
        iteration = 0
        newcoor = ball() * length
        while is_satisfied == False:
            iteration += 1
            if iteration > 10000:
                is_fail = True
                break
            is_satisfied = True
            newcoor = ball() * length
            if tooclose(newcoor, coor, rdmin, i) or toofar(newcoor, coor, rdmax, i):
                is_satisfied = False
        coor[i, :] = newcoor
    if is_fail == True:
        return -1
    else:
        return coor
