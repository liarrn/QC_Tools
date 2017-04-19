import os

def vasp_complete_check(directory):
    completed = [] #  have ' reached required accuracy' in OUTCAR
    completed_energy = []
    pending = [] # have INCAR in directory, don't have OUTCAR
    running = [] # have INCAR and OUTCAR, don't have  ' reached required accuracy'  in OUTCAR
    for dirpath, dirnames, filenames in os.walk(directory):
        if 'INCAR' in filenames:
            if 'OUTCAR' in filenames:
                with open(dirpath+'/OUTCAR', 'r') as fp:
                    out = fp.readlines()
                is_converged = False
                for nl in list(reversed(range(len(out)))):
                    if 'reached required accuracy' in out[nl]:
                        is_converged = True
                        break
                if is_converged == True:
                    completed.append(dirpath)
                    for nl in list(reversed(range(len(out)))):
                        if 'free' in out[nl]:
                            completed_energy.append(float(out[nl].split()[-2]))
                            break
                else:
                    running.append(dirpath)
            else:
                pending.append(dirpath)
    print 'COMPLETED: '
    for i in range(len(completed)):
        print '    %-30sfree energy: %.4f eV'%(completed[i], completed_energy[i])
    print 'PENDING (or have not been submited): '
    for i in pending:
        print '    '+i
    print 'RUNNING (or failed): '
    for i in running:
        print '    '+i
        
if __name__ == '__main__':
    vasp_complete_check('./')
