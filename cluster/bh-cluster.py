import numpy as np
from numpy import linalg as LA
from itertools import islice
import re
from subprocess import call
import math
import copy


def coor2xyz(elems, coor, des=''):
    xyz = ''
    natoms = coor.shape[0]
    xyz += '%d\n'%(natoms)
    xyz += des.rstrip()
    xyz += '\n'
    for i in range(natoms):
        xyz += '%-5s%15.9f%15.9f%15.9f\n'%(elems[i], coor[i][0], coor[i][1], coor[i][2])
    return xyz

def car2coor(carfile):
    carfile = [i.split() for i in carfile[4:-2]]
    elem = np.array([i[0] for i in carfile])
    coor = np.array([map(float, i[1: 4]) for i in carfile])
    return elem, coor

def xyz2coor(xyzfile):
    #     input xyz string list, each element in the list is a single line in xyz file
    #     output elem array and coor array with shape natoms*3
    natoms = int(xyzfile[0])
    xyzfile = xyzfile[2: ]
    xyzfile.sort()
    xyzfile = [i.split() for i in xyzfile]
    elems = [i[0] for i in xyzfile]
    coors = [map(float, i[1: ]) for i in xyzfile]
    elems = np.array(elems)
    coors = np.array(coors)
    return elems, coors

def coor2car(elems, coors):
    header = '!BIOSYM archive 3\nPBC=OFF\nMaterials Studio Generated CAR File\n'
    header += '!DATE     Apr 26 20:39:18 2016\n'
    end = 'end\nend\n'
    blanks = '                     '
    natoms = coors.shape[0]
    cartesian = ''
    for i in range(natoms):
        cartesian += "%-5s%15.9f%15.9f%15.9f%s%-2s%8.3f\n"%(elems[i], coors[i, 0], coors[i, 1], coors[i, 2], blanks, elems[i], 0.000)
    return header+cartesian+end

def xyz2car(xyzfile):
    #     input xyz string list, each element in the list is a single line in xyz file
    #     output car file for dmol calculations
    elems, coors = xyz2coor(xyzfile)
    header = '!BIOSYM archive 3\nPBC=OFF\nMaterials Studio Generated CAR File\n'
    header += '!DATE     Apr 26 20:39:18 2016\n'
    end = 'end\nend\n'
    blanks = '                     '
    natoms = coors.shape[0]
    cartesian = ''
    for i in range(natoms):
        cartesian += "%-5s%15.9f%15.9f%15.9f%s%-2s%8.3f\n"%(elems[i], coors[i, 0], coors[i, 1], coors[i, 2], blanks, elems[i], 0.000)
    return header+cartesian+end

def readxyzfile(fp, natoms):
    # read xyzfile
    # input file pointer, number of atoms
    # output: array of xyzfiles
    xyzfiles = []
    while True:
        next_n_lines = list(islice(fp, natoms+2))
        if not next_n_lines:
            break
        xyzfiles.append(next_n_lines)
    xyzfiles = np.array(xyzfiles)
    return xyzfiles
    
    
def writexyz(filename, write_type, elems, coor, des=''):
    with open(filename, write_type) as xyz_out_fp:
        xyz = coor2xyz(elems, coor, des)
        xyz_out_fp.write(xyz)
        
def writeline(filename, write_type, des):
    with open(filename, write_type) as fp:
        fp.write(des)


def extract_energy_dmol_opt(project_name):
    '''
    extract final energy from dmol3's geometry optimization calculation
    return energy in eV and elems, coor, force
    if not success, return -1, -1, -1, -1
    '''
    bohr2ang = 0.529177208
    ha_bohr2ev_ang = 51.42208619083232 
    ha2ev = 27.211396
    with open(project_name+'.outmol', 'r') as fp:
        outmol = fp.readlines()
    length = len(outmol)
    for opt_line in list(reversed(range(length))):
        if 'opt==' in outmol[opt_line]:
            break
    if opt_line == 0:
        return -1, -1, -1, -1
    energy = float(outmol[opt_line].split()[2])
    energy *= ha2ev
    for final_line in list(reversed(range(length))):
        if 'DERIVATIVES' in outmol[final_line]:
            break
    final_line_end = final_line
    while outmol[final_line_end] != '\n':
        final_line_end += 1
    final = copy.copy(outmol[final_line+2: final_line_end-1])
    final = [line.split()[1: ] for line in final]
    elems = np.array([line[0] for line in final])
    coor = np.array([map(float, line[1: 4]) for line in final])
    coor *= bohr2ang
    force = np.array([map(float, line[4: ]) for line in final])
    force *= ha_bohr2ev_ang
    return energy, elems, coor, force
    

def extract_energy(calculator, jobtype, project_name):
    '''
    input calculator, jobtype, filename
    calculator = enum('dmol3')
    jobtype = enum('sp', 'opt')
    if calc_params['calc_type'] is 'sp', then return the converged scf energy or -1 if scf failed
    if calc_params['calc_type'] is 'opt', then return energy in eV and elems, coor, force, or -1, -1, -1, -1 if failed
    '''
    if calculator == 'dmol3' and jobtype == 'sp':
        return extract_energy_dmol_sp(project_name)
    if calculator == 'dmol3' and jobtype == 'opt':
        return extract_energy_dmol_opt(project_name)
		
def call_dmol3(calc_params, elems, coor):
    '''
    input calc_params, elems, coor
    if calc_params['calc_type'] is 'sp', then return the converged scf energy or -1 if scf failed
    if calc_params['calc_type'] is 'opt', then return energy in eV and elems, coor, force, or -1, -1, -1, -1 if failed
    '''
    project_name = calc_params['project_name']
    calculator = 'dmol3'
    calc_call = calc_params['calc_call']
    calc_type = calc_params['calc_type']
    
    with open(project_name+'.car', 'w') as carfp:
        carfp.write(coor2car(elems, coor))
    call(calc_call)
    return extract_energy(calculator, calc_type, project_name)
	
def BoltzmannFactor(en1, en2, Temp, unit):
    # unit = enum('ha', 'ev')
    BoltzmannConstant = 1.38064852E-23 # in unit of Joule/Klein
    Joule2Hartree = 229371044869059200 # Joul2Hartree ha per joule
    Joule2EV = 6241506479963234000 # Joule2EV ev per joule
    if unit == 'ha':
        k = BoltzmannConstant * Joule2Hartree
    elif unit == 'ev':
        k = BoltzmannConstant * Joule2EV
    factor = math.exp(-(en1-en2)/(k*Temp))
    return factor
	
def move_cluster_randall(coor, params):
    natoms = coor.shape[0]
    step_size = params['step_size']
    moved = coor + (np.random.rand(natoms, 3)*2-1) * step_size
    return moved

def move_cluster(coor, params):
    # input: coordination, method, params
    # output: displaced coor
    # method = enum('randall)
    if params['method'] == 'randall':
        return move_cluster_randall(coor, params)
    if params['method'] == 'randone':
        return move_cluster_randone(coor, params)
		
def basin_hopping_dmol3(xyzfilename, calc_params, mc_params, move_params):
    calculator = 'dmol3'
    temp = mc_params['temp']
    num_steps = mc_params['num_steps']
    project_name = calc_params['project_name']
    natoms = calc_params['natoms']
    xyz_out = project_name+'-xyz.xyz'
    energy_outfile_name = project_name+'-energy.out'
    
    with open(xyzfilename, 'r') as xyzfp:
        xyz = readxyzfile(xyzfp, natoms)
    curr_elems, curr_coor = xyz2coor(xyz[0])
    curr_energy, curr_elems, curr_coor, curr_force = call_dmol3(calc_params, curr_elems, curr_coor)
    curr_adsorbed_energy = curr_energy + calc_CO_tot(curr_coor)
    
    writexyz(xyz_out, 'w', curr_elems, curr_coor, '1th %.9f eV'%curr_energy)
    writeline(energy_outfile_name, 'w', 'current_energy(eV)    next_energy(eV)    current_adsorbed_energy(eV)    next_adsorbed_energy(eV)\n')
    writeline(energy_outfile_name, 'a', '%.9f    %.9f\n'%(curr_energy, curr_energy, curr_adsorbed_energy, curr_adsorbed_energy))
    
    for i in range(num_steps):
        next_coor = move_cluster(curr_coor, move_params)
        next_elems = copy.copy(curr_elems)
        next_energy, next_elems, next_coor, next_force = call_dmol3(calc_params, next_elems, next_coor)
        
        # if scf run is not sucessful, continue to next iteration
        if next_energy == -1:
            continue
        next_adsorbed_energy = next_energy + calc_CO_tot(next_coor)
        writexyz(xyz_out, 'a', next_elems, next_coor, '%dth %.9f eV'%(i+2, next_energy))
        writeline(energy_outfile_name, 'a', '%.9f    %.9f\n'%(curr_energy, next_energy, curr_adsorbed_energy, next_adsorbed_energy))
        
        is_move = False
        if next_adsorbed_energy<= curr_adsorbed_energy:
            is_move = True
        else:
            bf = BoltzmannFactor(next_adsorbed_energy, curr_adsorbed_energy, temp, 'ev')
            rn = np.random.rand()
            if rn < bf:
                is_move = True
        if is_move == True:
            curr_energy = next_energy
            curr_adsorbed_energy = next_adsorbed_energy
            curr_coor = next_coor
            curr_elems = next_elems

def basin_hopping(xyzfilename, calc_params, mc_params, move_params):
    if calc_params['calculator'] == 'dmol3':
        basin_hopping_dmol3(xyzfilename, calc_params, mc_params, move_params)
        





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


def in_plane(coor):
    criterion = 0.3
    if coor.shape[0] <= 3:
        return True
    norm = np.cross(coor[1, :] - coor[0, :], coor[2, : ]- coor[0, :])
    norm /= LA.norm(norm)
    is_in_plane = True
    for i in range(3, coor.shape[0]):
        vec = coor[i, :] - coor[0, :]
        vec /= LA.norm(vec)
#         print np.dot(vec, norm)
        if np.abs(np.dot(vec, norm)) > criterion:
            is_in_plane = False
            break
    return is_in_plane



def getSurfAtoms(cluster):
    # return surface atoms, started with zero
    criterion = 1.82
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
        else:
            if in_plane(cluster[nnatoms]):
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
        
def calc_CO_tot(coor):
    surfatoms = getSurfAtoms(coor)
    COad = np.array([CalcAd(i, coor)[0] for i in surfatoms])
    return sum(COad)
        
        
        
        
        
        
        
        
        
        
        
        
if __name__ == '__main__':
    xyzfilename = './au18-anion.xyz'
    calc_params = {}
    calc_params['project_name'] = 'au18-anion'
    calc_params['natoms'] = 18
    calc_params['calculator'] = 'dmol3'
    calc_params['calc_call'] = './RunDMol3.sh -np 24 ' + calc_params['project_name']
    calc_params['calc_call'] = calc_params['calc_call'].split()
    calc_params['calc_type'] = 'opt'
    
    mc_params = {}
    mc_params['temp'] = 7500
    mc_params['num_steps'] = 10000
    
    move_params={'method': 'randall', 'step_size': 1.2}
    
    basin_hopping_dmol3(xyzfilename, calc_params, mc_params, move_params)
    