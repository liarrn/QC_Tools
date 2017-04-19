import numpy as np
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

def coor2car(elems, coors, lattice_params={'PBC': False}):
    PBC = lattice_params['PBC']
    header = '!BIOSYM archive 3\n'
    if PBC == True:
        header += 'PBC=ON\n'
    elif PBC == False:
        header += 'PBC=OFF\n'
    header += 'Materials Studio Generated CAR File\n'
    header += '!DATE     Apr 26 20:39:18 2016\n'
    if PBC == True:
        header += 'PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f (P1)\n'\
        %(lattice_params['a'], lattice_params['b'], lattice_params['c'], lattice_params['b2c'], lattice_params['a2c'], lattice_params['a2b'])
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
        
def extract_energy_dmol_sp(project_name):
    '''
    extract final energy from dmol3's static point calculation
    return energy in eV
    if not success, return -1
    '''
    ha2ev = 27.211396
    fp = open(project_name+'.outmol', 'r')
    outmol = fp.read()
    fp.close()
    is_success = re.findall("successfully", outmol)
    if is_success != []:
        is_success = True
    else:
        is_success = False
        return -1
    outmol = outmol.splitlines()
    Ef_lines = []
    for line in outmol:
        if 'Ef' in line:
            Ef_lines.append(line)
    energy = float(Ef_lines[-1].split()[1][:-2])
    energy *= ha2ev
    return energy

# def extract_energy_dmol_opt(project_name):
#     '''
#     extract final energy from dmol3's geometry optimization calculation
#     return energy in eV and elems, coor, force
#     if not success, return -1, -1, -1, -1
#     '''
#     bohr2ang = 0.529177208
#     ha_bohr2ev_ang = 51.42208619083232 
#     ha2ev = 27.211396
#     with open(project_name+'.outmol', 'r') as fp:
#         outmol = fp.readlines()
#     length = len(outmol)
#     for opt_line in list(reversed(range(length))):
#         if 'opt==' in outmol[opt_line]:
#             break
#     if opt_line == 0:
#         return -1, -1, -1, -1
#     energy = float(outmol[opt_line].split()[2])
#     energy *= ha2ev
#     for final_line in list(reversed(range(length))):
#         if 'DERIVATIVES' in outmol[final_line]:
#             break
#     final_line_end = final_line
#     while outmol[final_line_end] != '\n':
#         final_line_end += 1
#     final = copy.copy(outmol[final_line+2: final_line_end-1])
#     final = [line.split()[1: ] for line in final]
#     elems = np.array([line[0] for line in final])
#     coor = np.array([map(float, line[1: 4]) for line in final])
#     coor *= bohr2ang
#     force = np.array([map(float, line[4: ]) for line in final])
#     force *= ha_bohr2ev_ang
#     return energy, elems, coor, force



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
        if 'Ef' in outmol[opt_line]:
            break
    if opt_line == 0:
        return -1, -1, -1, -1
    energy = float(outmol[opt_line].split()[1][:-2])
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
    
def call_dmol3(calc_params, elems, coor, lattice_param = {'PBC': False}):
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
        carfp.write(coor2car(elems, coor, lattice_param))
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
    
def elem_index(elems, target_elem):
    result = []
    for i in range(len(elems)):
        if elems[i] == target_elem:
            result.append(i)
    return result

# elems = ['Au', 'Au', 'Au', 'Au', 'Au', 'Au', 'S', 'S', 'S', 'S', 'S', 'S', 'H', 'H', 'H', 'H', 'H', 'H']
# elem_index(elems, 'H')


def nearest_hydrogen(elems, coors, sulfur):
    hydrogen = elem_index(elems, 'H')
    shortest_dist = 100
    for h in hydrogen:
        dist = np.sqrt(sum((coors[sulfur]-coors[h])**2))
        if dist < shortest_dist:
            shortest_dist = dist
            nearest_hydrogen = h
    return nearest_hydrogen
    
def move_cluster_randall(coor, params):
    natoms = coor.shape[0]
    step_size = params['step_size']
    moved = coor + (np.random.rand(natoms, 3)*2-1) * step_size
    return moved
     
def move_cluster_all(coor, params):
    natoms = coor.shape[0]
    step_size = params['step_size']
    moved = coor + (np.tile(np.random.rand(3), [natoms, 1])*2-1) * step_size
    return moved

def move_cluster_randone(coor, params):
    natoms = coor.shape[0]
    step_size = params['step_size']
    chosen_atom = int(np.random.rand()*10000) % natoms
    moved = copy.copy(coor)
    moved[chosen_atom] += (np.random.rand(3)*2-1) * step_size
    return moved

def move_cluster(coor, params):
    # input: coordination, method, params
    # output: displaced coor
    # method = enum('randall)
    if params['method'] == 'randall':
        return move_cluster_randall(coor, params)
    if params['method'] == 'randone':
        return move_cluster_randone(coor, params)
    

def move_cluster_thiolate(elems, coor, params):
    au = elem_index(elems, 'Au')
    sulfur_hydrogen_pair = [[s, nearest_hydrogen(elems, coor, s)] for s in elem_index(elems, 'S')]
#     print au
#     print sulfur_hydrogen_pair
    moved = copy.copy(coor)
    # random move au atoms
    moved[au] = move_cluster_randall(coor[au], params)
    for i in range(len(sulfur_hydrogen_pair)):
        moved[sulfur_hydrogen_pair[i]] = move_cluster_all(coor[sulfur_hydrogen_pair[i]], params)
    return elems, moved
    
def basin_hopping_dmol3(xyzfilename, calc_params, mc_params, move_params):
    calculator = 'dmol3'
    temp = mc_params['temp']
    num_steps = mc_params['num_steps']
    project_name = calc_params['project_name']
    natoms = calc_params['natoms']
    xyz_out = project_name+'-xyz.xyz'
    xyz_force_out = project_name+'-xyz-force.xyz'
    energy_outfile_name = project_name+'-energy.out'
    
    with open(xyzfilename, 'r') as xyzfp:
        xyz = readxyzfile(xyzfp, natoms)
    curr_elems, curr_coor = xyz2coor(xyz[0])
    curr_energy, curr_elems, curr_coor, curr_force = call_dmol3(calc_params, curr_elems, curr_coor)
    
    writexyz(xyz_out, 'w', curr_elems, curr_coor, '1th %.9f eV'%curr_energy)
    writeline(energy_outfile_name, 'w', 'current_energy(eV)    next_energy(eV)\n')
    writeline(energy_outfile_name, 'a', '%.9f    %.9f\n'%(curr_energy, curr_energy))
    
    for i in range(num_steps):
        next_elems, next_coor = move_cluster_thiolate(curr_elems, curr_coor, move_params)
        # next_elems = copy.copy(curr_elems)
        next_energy, next_elems, next_coor, next_force = call_dmol3(calc_params, next_elems, next_coor)
        
        # if scf run is not sucessful, continue to next iteration
        if next_energy == -1:
            continue
        
        writexyz(xyz_out, 'a', next_elems, next_coor, '%dth %.9f eV'%(i+2, next_energy))
        writeline(energy_outfile_name, 'a', '%.9f    %.9f\n'%(curr_energy, next_energy))
        
        is_move = False
        if next_energy<= curr_energy:
            is_move = True
        else:
            bf = BoltzmannFactor(next_energy, curr_energy, temp, 'ev')
            rn = np.random.rand()
            if rn < bf:
                is_move = True
        if is_move == True:
            curr_energy = next_energy
            curr_coor = next_coor
            curr_elems = next_elems

def basin_hopping(xyzfilename, calc_params, mc_params, move_params):
    if calc_params['calculator'] == 'dmol3':
        basin_hopping_dmol3(xyzfilename, calc_params, mc_params, move_params)
        
        
if __name__ == '__main__':
    xyzfilename = './Au7S6H6.xyz'
    calc_params = {}
    calc_params['project_name'] = 'Au7S6H6'
    calc_params['natoms'] = 19
    calc_params['calculator'] = 'dmol3'
    calc_params['calc_call'] = './RunDMol3.sh -np 24 ' + calc_params['project_name']
    calc_params['calc_call'] = calc_params['calc_call'].split()
    calc_params['calc_type'] = 'opt'
    
    mc_params = {}
    mc_params['temp'] = 7500
    mc_params['num_steps'] = 10000
    
    move_params={'method': 'randall', 'step_size': 0.7}
    
    basin_hopping_dmol3(xyzfilename, calc_params, mc_params, move_params)
    
