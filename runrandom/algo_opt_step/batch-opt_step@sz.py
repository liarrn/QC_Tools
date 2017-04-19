import numpy as np
from itertools import islice
import re
from subprocess import call
import math
import copy
from subprocess import check_output
import time

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
		

def extract_energy_dmol_opt(project_name):
    '''
    extract final energy from dmol3's geometry optimization calculation
    return energy in eV and elems, coor, force
    if not success, return -1, -1, -1, -1
    '''
    def get_coor(deri_line, opt_line):
        energy = float(outmol[opt_line].split()[2])
        energy *= ha2ev
        end_line = deri_line
        while outmol[end_line] != '\n':
            end_line += 1
        final = copy.copy(outmol[deri_line+2: end_line-1])
        final = [line.split()[1: ] for line in final]
        elems = np.array([line[0] for line in final])
        coor = np.array([map(float, line[1: 4]) for line in final])
        coor *= bohr2ang
        force = np.array([map(float, line[4: ]) for line in final])
        force *= ha_bohr2ev_ang
        result = {'elems': elems, 'coor': coor, 'force': force, 'energy': energy}
        return result
    
    bohr2ang = 0.529177208
    ha_bohr2ev_ang = 51.42208619083232 
    ha2ev = 27.211396
    with open(project_name+'.outmol', 'r') as fp:
        outmol = fp.readlines()
    opt_lines = []
    deri_lines = []
    for ln in range(len(outmol)):
        if 'opt==' in outmol[ln]:
            opt_lines.append(ln)
        elif 'DERIVATIVES' in outmol[ln]:
            deri_lines.append(ln)
    opt_lines = opt_lines[2: ]
    result = [get_coor(deri_lines[i], opt_lines[i]) for i in range(len(opt_lines))]
    return result
    

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
    np = calc_params['np']
    calculator = 'dmol3'
    calc_type = calc_params['calc_type']
    
    with open(project_name+'.car', 'w') as carfp:
        carfp.write(coor2car(elems, coor))
    run_Dmol(project_name, np)
    return extract_energy(calculator, calc_type, project_name)

def subjob(sub_command):
    sub_command = sub_command.split()
    subout = check_output(sub_command)
    jobid = int(subout.split()[1][1:-1])
    return jobid

def is_running(jobid):
    bjobs_command = 'bjobs'.split()
    bjobs_out = check_output(bjobs_command).split('\n')[1:]
    bjobs_out = [line.split() for line in bjobs_out]
    all_jobid = [int(line[0]) for line in bjobs_out if len(line)>1]
    running = jobid in all_jobid
    return running

def wait_until_finish(jobid):
    while True:
        time.sleep(30)
        # print 'running'
        running = is_running(jobid)
        if running == False:
            return

def run_Dmol(project_name, np=36):
    lsf = 'APP_NAME=intelw_exc\n'
    lsf += 'NP=%d\n'%np
    lsf += 'NP_PER_NODE=12\n'
    lsf += 'RUN="RAW"\n'
    lsf += 'filename=%s\n'%project_name
    lsf += '/home-yw/Soft/msi/MS70/MaterialsStudio7.0/etc/DMol3/bin/RunDMol3.sh -np $NP $filename\n'
    with open('MS70.lsf', 'w') as fp:
        fp.write(lsf)
    sub_command = 'bsub MS70.lsf'
    jobid = subjob(sub_command)
    # print 'jobid = %d'%jobid
    wait_until_finish(jobid)
    # print 'finish'
    

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

def runrandom(calc_params):
    # calculator = enum('dmol3')
    project_name = calc_params['project_name']
    natoms = calc_params['natoms']
    xyz_force_out = project_name+'-xyz-force.xyz'
    xyz_out = project_name+'-xyz.xyz'
    # create output file, if existed, wash them out
    open(xyz_force_out, 'w').close()
    open(xyz_out, 'w').close()
    for step in range(calc_params['num_steps']):
        coor = randcluster(natoms, calc_params['rdmin'], calc_params['rdmax'])
        if type(coor) == int:
            continue
        elems = ['Au'] * natoms
        next = call_dmol3(calc_params, elems, coor)
        next_energy = next[-1]['energy']
        # next_energy, next_elems, next_coor, next_force = call_dmol3(calc_params, next_elems, next_coor)
        
        for opt_step in range(len(next)):
            writexyz(xyz_out, 'a', next[opt_step]['elems'], next[opt_step]['coor'], '%dth %.9f eV'%(step+1, next[opt_step]['energy']))
            writexyzforce(xyz_force_out, 'a', next[opt_step]['elems'], next[opt_step]['coor'], next[opt_step]['force'], '%dth %.9f eV'%(step+1, next[opt_step]['energy']))
        # writeline(energy_outfile_name, 'a', '%.9f    %.9f\n'%(curr_energy, next_energy))

		
if __name__ == '__main__':
    calc_params = {}
    calc_params['num_steps'] = 1000
    calc_params['calculator'] = 'dmol3'
    calc_params['rdmin'] = 2.6
    calc_params['rdmax'] = 3.2
    calc_params['calc_type'] = 'opt'
    
    for natoms in range(8, 9):
        movecall = 'cp opt.input au%d.input'%natoms
        call(movecall.split())
        calc_params['project_name'] = 'au%d'%natoms
        calc_params['natoms'] = natoms
        calc_params['np'] = 24
        runrandom(calc_params)
