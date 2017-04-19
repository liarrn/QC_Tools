

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

