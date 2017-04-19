import numpy as np
from itertools import islice
import re
from subprocess import call
import math
import copy

def coor2xsf(elems, coor, lattice_param={'PBC': False}, force=False, des='#'):
    natoms = elems.shape[0]
    header = des.rstrip()+'\n\n'
    if lattice_param['PBC'] == False:
        header += 'ATOMS\n'
    elif lattice_param['PBC'] == True:
        lattice_param['a2b'] = lattice_param['a2b']/180.0*np.pi
        lattice_param['b2c'] = lattice_param['b2c']/180.0*np.pi
        lattice_param['a2c'] = lattice_param['a2c']/180.0*np.pi
        va = np.array([lattice_param['a'], 0, 0])
        vb = np.array([lattice_param['b']*np.cos(lattice_param['a2b']), lattice_param['b']*np.sin(lattice_param['a2b']), 0])
        vc_x= np.cos(lattice_param['a2c'])
        vc_y = (np.cos(lattice_param['b2c'])-np.cos(lattice_param['a2c'])*np.cos(lattice_param['a2b']))/np.sin(lattice_param['a2b'])
        vc_z = np.sqrt(1-vc_x**2-vc_y**2)
        vc = np.array([vc_x, vc_y, vc_z])*lattice_param['c']
        header += 'CRYSTAL\nPRIMVEC\n'
        header += '%15.9f%15.9f%15.9f\n'%(va[0], va[1], va[2])
        header += '%15.9f%15.9f%15.9f\n'%(vb[0], vb[1], vb[2])
        header += '%15.9f%15.9f%15.9f\n'%(vc[0], vc[1], vc[2])
        header += 'PRIMCOORD\n'
        header+= '%d    %d\n'%(coor.shape[0], 1)
    content = ''
    if type(force)!=np.ndarray:
        for i in range(natoms):
            content += '%-4s%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n'\
            %(elems[i], coor[i, 0], coor[i, 1], coor[i, 2], 0.0, 0.0, 0.0)
    else:
        for i in range(natoms):
            content += '%-4s%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n'\
            %(elems[i], coor[i, 0], coor[i, 1], coor[i, 2], force[i, 0], force[i, 1], force[i, 2])
    return header+content

def coor2xyzforce(elems, coor, force, des=''):
    xyz = ''
    natoms = coor.shape[0]
    xyz += '%d\n'%(natoms)
    xyz += des.rstrip()
    xyz += '\n'
    for i in range(natoms):
        xyz += '%-5s%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n'%(elems[i], coor[i][0], coor[i][1], coor[i][2], force[i][0], force[i][1], force[i][2])
    return xyz

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
	
def writexyzforce(filename, write_type, elems, coor, force, des=''):
    with open(filename, write_type) as xyz_out_fp:
        xyz = coor2xyzforce(elems, coor, force, des)
        xyz_out_fp.write(xyz)
    
    
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
	
	
def genfcc(lattice_param):
    '''
    input lattice_param, return corresponding fcc coordinates
    '''
    a = lattice_param['a']
    b = lattice_param['b']
    c = lattice_param['c']
    a2b = lattice_param['a2b']/180.0*np.pi
    b2c = lattice_param['b2c']/180.0*np.pi
    a2c = lattice_param['a2c']/180.0*np.pi
    
    va = np.array([a, 0, 0])
    vb = np.array([b*np.cos(a2b), b*np.sin(a2b), 0])
    vc_x= np.cos(a2c)
    vc_y = (np.cos(b2c)-np.cos(a2c)*np.cos(a2b))/np.sin(a2b)
    vc_z = np.sqrt(1-vc_x**2-vc_y**2)
    vc = np.array([vc_x, vc_y, vc_z])*c
#     print va, vb, vc
    
    coor = np.zeros([4, 3])
    coor[1, :] = (vb+vc)/2.0
    coor[2, :] = (va+vc)/2.0
    coor[3, :] = (va+vb)/2.0
    return coor

# def lattice_constant_scaling(lattice_param, a_range, b_range, c_range):
#     # return list of coor and list of lattice_param
#     coor = np.array([genfcc(lattice_param)])
#     lattice_params = np.array([lattice_param])
#     curr_lattice = copy.copy(lattice_param)
#     for a in a_range+lattice_param['a']:
#         curr_lattice['a'] = a
#         for b in b_range+lattice_param['b']:
#             curr_lattice['b'] = b
#             for c in c_range+lattice_param['c']:
#                 curr_lattice['c'] = c
#                 coor = np.concatenate([coor, np.array([genfcc(curr_lattice)])])
#                 lattice_params = np.concatenate([lattice_params, np.array([copy.copy(curr_lattice)])])
# #                 print curr_lattice
# #                 print lattice_params
#     coor = coor[1: ]
#     lattice_params = lattice_params[1: ]
#     return coor, lattice_params

# lattice_param = {'a': 4.4564, 'b': 4.4564, 'c': 4.4564, 'a2b': 90, 'b2c': 90, 'a2c': 90}
# a_range = np.arange(-0.6, 0.7, 0.1)
# b_range = np.arange(-0.6, 0.7, 0.1)
# c_range = np.arange(-0.6, 0.7, 0.1)
# coor, lattice_param = lattice_constant_scaling(lattice_param, a_range, b_range, c_range)


# def lattice_constant_scaling(lattice_param, scaling_params):
#     # return list of coor and list of lattice_param
#     coor = np.array([genfcc(lattice_param)])
#     lattice_params = np.array([lattice_param])
#     curr_lattice = copy.copy(lattice_param)
#     for a in np.arange(scaling_params['a_min'], scaling_params['a_max']+scaling_params['a_step'], scaling_params['a_step']):
#         curr_lattice['a'] = a
#         for b in np.arange(a, scaling_params['b_max']+scaling_params['b_step'], scaling_params['b_step']):
#             curr_lattice['b'] = b
#             for c in np.arange(b, scaling_params['c_max']+scaling_params['c_step'], scaling_params['c_step']):
#                 curr_lattice['c'] = c
#                 coor = np.concatenate([coor, np.array([genfcc(curr_lattice)])])
#                 lattice_params = np.concatenate([lattice_params, np.array([copy.copy(curr_lattice)])])
#     coor = coor[1: ]
#     lattice_params = lattice_params[1: ]
#     return coor, lattice_params

# lattice_param = {'a': 4.4564, 'b': 4.4564, 'c': 4.4564, 'a2b': 90, 'b2c': 90, 'a2c': 90}
# scaling_params = {'a_min': 3.8564, 'a_max': 5.0564, 'a_step': 0.1, 
#                   'b_min': 3.8564, 'b_max': 5.0564, 'b_step': 0.1,
#                   'c_min': 3.8564, 'c_max': 5.0564, 'c_step': 0.1}
# coor, lattice_param = lattice_constant_scaling(lattice_param, scaling_params)


def lattice_constant_scaling(lattice_param, scaling_params):
    # return list of coor and list of lattice_param
    coor = np.array([genfcc(lattice_param)])
    lattice_params = np.array([lattice_param])
    curr_lattice = copy.copy(lattice_param)
    for scaling in scaling_params['scaling']:
        curr_lattice['a'] = lattice_param['a'] + scaling
        curr_lattice['b'] = lattice_param['b'] + scaling
        curr_lattice['c'] = lattice_param['c'] + scaling
        coor = np.concatenate([coor, np.array([genfcc(curr_lattice)])])
        lattice_params = np.concatenate([lattice_params, np.array([copy.copy(curr_lattice)])])
    coor = coor[1: ]
    lattice_params = lattice_params[1: ]
    return coor, lattice_params

# lattice_param = {'a': 4.4564, 'b': 4.4564, 'c': 4.4564, 'a2b': 90, 'b2c': 90, 'a2c': 90}
# scaling_params = {'scaling': np.arange(-0.6, 0.7, 0.1)}
# coor, lattice_param = lattice_constant_scaling(lattice_param, scaling_params)
# print lattice_param.shape

def monoclinic_strain(lattice_param, mstrain_params):
    # return list of coor and list of lattice_param
    coor = np.array([genfcc(lattice_param)])
    lattice_params = np.array([lattice_param])
    curr_lattice = copy.copy(lattice_param)
    c_origin = lattice_param['c']
    for disp in np.arange(mstrain_params['left'], mstrain_params['right']+mstrain_params['step'], mstrain_params['step']):
        curr_lattice['a2c'] = np.arctan(c_origin/disp)/np.pi*180
        if curr_lattice['a2c'] < 0:
            curr_lattice['a2c'] = 180+curr_lattice['a2c']
        curr_lattice['c'] = np.sqrt(c_origin**2+disp**2)
        coor = np.concatenate([coor, np.array([genfcc(curr_lattice)])])
        lattice_params = np.concatenate([lattice_params, np.array([copy.copy(curr_lattice)])])
    coor = coor[1: ]
    lattice_params = lattice_params[1: ]
    return coor, lattice_params

# lattice_param = {'a': 4.4564, 'b': 4.4564, 'c': 4.4564, 'a2b': 90, 'b2c': 90, 'a2c': 90}
# mstrain_params = {'left': -1.6, 'right': 1.6, 'step': 0.2}
# coor, lattice_param = monoclinic_strain(lattice_param, mstrain_params)
# print lattice_param.shape

def orthorhombic_strain(lattice_param, ostrain_params):
    # return list of coor and list of lattice_param
    coor = np.array([genfcc(lattice_param)])
    lattice_params = np.array([lattice_param])
    curr_lattice = copy.copy(lattice_param)
    ac = lattice_param['a']*lattice_param['c']
    for strain in ostrain_params['strain']:
        curr_lattice['a'] = lattice_param['a'] + strain
        curr_lattice['c'] = ac / curr_lattice['a']
        coor = np.concatenate([coor, np.array([genfcc(curr_lattice)])])
        lattice_params = np.concatenate([lattice_params, np.array([copy.copy(curr_lattice)])])
    coor = coor[1: ]
    lattice_params = lattice_params[1: ]
    return coor, lattice_params

# lattice_param = {'a': 4.4564, 'b': 4.4564, 'c': 4.4564, 'a2b': 90, 'b2c': 90, 'a2c': 90}
# ostrain_params = {'strain': np.arange(-0.6, 0.7, 0.1)}
# coor, lattice_param = orthorhombic_strain(lattice_param, ostrain_params)
# print lattice_param.shape


if __name__ == '__main__':
    lattice_param = {'a': 4.4564, 'b': 4.4564, 'c': 4.4564, 'a2b': 90, 'b2c': 90, 'a2c': 90, 'PBC': True}
    
    scaling_params = {'scaling': np.arange(-0.6, 0.7, 0.1)}
    mstrain_params = {'left': -1.6, 'right': 1.6, 'step': 0.2}
    ostrain_params = {'strain': np.arange(-0.6, 0.7, 0.1)}
    
    coor_scaling, lattice_param_scaling = lattice_constant_scaling(lattice_param, scaling_params)
    
    all_lattice = copy.copy(lattice_param_scaling)
    all_coor = copy.copy(coor_scaling)
    
    for i in range(lattice_param_scaling.shape[0]):
        coor_monoclinic_strain, lattice_param_monoclinic_strain = monoclinic_strain(lattice_param_scaling[i], mstrain_params)
        all_coor = np.concatenate([all_coor, coor_monoclinic_strain])
        all_lattice = np.concatenate([all_lattice, lattice_param_monoclinic_strain])
        
    for i in range(lattice_param_scaling.shape[0]):
        coor_orthorhombic_strain, lattice_param_orthorhombic_strain = orthorhombic_strain(lattice_param_scaling[i], ostrain_params)
        all_coor = np.concatenate([all_coor, coor_orthorhombic_strain])
        all_lattice = np.concatenate([all_lattice, lattice_param_orthorhombic_strain])  
        
    prefix = '4atom_bulk_'
    calc_params = {}
    calc_params['calculator'] = 'dmol3'
    calc_params['calc_type'] = 'opt'
    calc_params['natoms'] = 4
    elems = np.array(['Au']*4)
    for i in range(403):
        movecall = 'cp opt.input %s%d.input'%(prefix, i)
        call(movecall.split())
        calc_params['project_name'] = '%s%d'%(prefix, i)
        calc_params['calc_call'] = './RunDMol3.sh -np 24 ' + calc_params['project_name']
        calc_params['calc_call'] = calc_params['calc_call'].split()
        energy, elems, coor, force = call_dmol3(calc_params, elems, all_coor[i], all_lattice[i])
        # print energy, elems, coor, force
        xsf = coor2xsf(elems, coor, all_lattice[i], force, des='# total energy = %.9f eV'%energy)
		with open(calc_params['project_name']+'.xsf', 'w') as fp:
            fp.write(xsf)
