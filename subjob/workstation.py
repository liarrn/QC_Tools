import paramiko
import time
import math
import numpy as np

##########################WORKSTATION##########################
ncores = 48
idle_criteria = 10
ipadd = '10.15.1.30'
username='liyd'
password='12315'
my_jobs = ['cd /vol-th/home/Chgy/liyd/AuNT/zigzag/30/2/5.52/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 302f RunVASP.5.4.1.sh', \
          'cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/2/9.3/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 332f RunVASP.5.4.1.sh', \
          'cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/3/9.31/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 333f RunVASP.5.4.1.sh']

ssh = paramiko.client.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(ipadd, username=username, password=password)
while True:
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command('sar -u 10 1')
    cpuinfo = ssh_stdout.readlines()
    idle_cpu = int(np.round(ncores*float(cpuinfo[-1].split()[-1]))/100.0)
    print idle_cpu, len(my_jobs)
    if idle_cpu < idle_criteria:
        continue
    if len(my_jobs) == 0:
        break
    a, b, c = ssh.exec_command(my_jobs[0])
    my_jobs = my_jobs[1: ]