from subprocess import check_output
import time

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
    

            
if __name__ == '__main__':
    project_name = 'au4'
    np = 12
    run_Dmol(project_name, np)

    
    
    
    
    