import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

def DOS(filename, LDOS = False, atom_list=[]):
    with open(filename, mode = 'r') as fp:
        dos = fp.readlines()
        num_of_state = int(dos[5].split()[2])
        fermi_energy = float(dos[5].split()[3])
        total_dos_block = dos[6 : 6+num_of_state]
        total_dos_spin_up = np.array([ map(float, i.split())  for i in total_dos_block])[:, 1]
        total_dos_spin_down = np.array([ map(float, i.split())  for i in total_dos_block])[:, 2]
        dos_energy = np.array([ map(float, i.split())  for i in total_dos_block])[:, 0] - fermi_energy
        total_dos = total_dos_spin_up + total_dos_spin_down
        if not LDOS:
            return dos_energy, total_dos_spin_up, total_dos_spin_down, total_dos
        else:
            dos_sum = np.zeros(num_of_state)
            dos_s = np.zeros(num_of_state)
            dos_p = np.zeros(num_of_state)
            dos_d = np.zeros(num_of_state)
            for atom_num in atom_list:
                dos_atom_num = dos[6+atom_num*(num_of_state+1) : 5+(atom_num+1)*(num_of_state+1)]
                dos_atom_num = np.array([ map(float, i.split())  for i in dos_atom_num])
                dos_sum = dos_sum +  np.sum(dos_atom_num[:, 1:], axis = 1)
                dos_s = dos_s + np.sum(dos_atom_num[:, 1:3], axis = 1)
                dos_p = dos_p + np.sum(dos_atom_num[:, 3:9], axis = 1)
                dos_d = dos_d + np.sum(dos_atom_num[:, 9:19], axis = 1)
            return dos_energy, total_dos_spin_up, total_dos_spin_down, total_dos, dos_sum, dos_s, dos_p, dos_d

if __name__ == '__main__':
	dos_energy, total_dos_spin_up, total_dos_spin_down, total_dos, dos_sum, dos_s, dos_p, dos_d=  DOS('DOSCAR', True, range(1,5))
	plt.plot(dos_energy, total_dos, 'ro-', dos_energy, dos_sum, 'bo-')