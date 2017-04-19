import numpy as np
import os

def read_fix_list(filename, n):
    with open(filename, 'r') as fp:
        xyz = fp.readlines()
    xyz = xyz[8:]
    read_z = [float(i.split()[2])  for i in xyz]
    fix_list = np.unique(read_z)
    return fix_list[:n]

def fix_atom_selected(source_file, target_file, fix_list):
    fixed_atom_num = 0
    row_sign = 0
    selected = 0
    with open(source_file, mode = 'r') as initial_file:
        with open(target_file, mode = 'w') as finished_file:
            for line in initial_file:
                if row_sign < 7:
                    finished_file.write(line)
                elif row_sign == 7:
                    if list(line.strip())[0] == 's' or list(line.strip())[0] == 'S':
                        selected = 1
                    else:
                        selected = 0
                        finished_file.write('selective\n')
                    finished_file.write(line)
                else:
                    if selected:
                        finished_file.write(line)
                    else:
                        finished_file.write(line.rstrip())
                        z=float(line.split()[2])
                        if z in fix_list:
                            finished_file.write(' F F F \n')
                            fixed_atom_num +=1
                        else:
                            finished_file.write(' T T T \n')
                row_sign += 1
            print 'there are %d atoms fixed' % fixed_atom_num
        return target_file

def get_vaspfile(infile):
    if infile != '':
        return infile
    filename_list = os.listdir('./')
    for filename in filename_list:
        if '.vasp' in filename:
            return filename
 
if __name__ == '__main__':
    num_layer_fixed = 3
    input_file = ''
    input_file = get_vaspfile(input_file)
    print 'input file: %s'%input_file
    output_file = 'POSCAR'
    a = fix_atom_selected(input_file, output_file, read_fix_list(input_file, num_layer_fixed))