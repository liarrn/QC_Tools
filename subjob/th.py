import paramiko
import time

# my_jobs = [['cd /vol-th/home/Chgy/liyd/AuNT/zigzag/20/3/5.52/opt_fine', 'yhbatch -N 6 -n 36 -p TH_NET -J 203f RunVASP.5.4.1.sh'], \
#           ['cd /vol-th/home/Chgy/liyd/AuNT/zigzag/30/1/5.5/opt_fine', 'yhbatch -N 6 -n 36 -p TH_NET -J 301f RunVASP.5.4.1.sh'], \
#           ['cd /vol-th/home/Chgy/liyd/AuNT/zigzag/30/2/5.52/opt_fine', 'yhbatch -N 6 -n 36 -p TH_NET -J 302f RunVASP.5.4.1.sh'], \
#           ['cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/2/9.3/opt_fine', 'yhbatch -N 6 -n 36 -p TH_NET -J 332f RunVASP.5.4.1.sh'], \
#           ['cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/3/9.31/opt_fine', 'yhbatch -N 6 -n 36 -p TH_NET -J 333f RunVASP.5.4.1.sh']]

# my_jobs = ['cd /vol-th/home/Chgy/liyd/AuNT/zigzag/30/2/5.52/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 302f RunVASP.5.4.1.sh', \
          # 'cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/2/9.3/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 332f RunVASP.5.4.1.sh', \
          # 'cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/3/9.31/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 333f RunVASP.5.4.1.sh']

# my_jobs = ['cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/2/9.3/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 332f RunVASP.5.4.1.sh', \
          # 'cd /vol-th/home/Chgy/liyd/AuNT/armchair/33/3/9.31/opt_fine; yhbatch -N 6 -n 36 -p TH_NET -J 333f RunVASP.5.4.1.sh']
          
# my_jobs = ['cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Cu-H2/100cell111; yhbatch -N 2 -n 24 -p TH_NET -J H100 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Cu-H2/110cell111; yhbatch -N 2 -n 24 -p TH_NET -J H110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Cu-H2/111cell111; yhbatch -N 2 -n 24 -p TH_NET -J H111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Cu-N2/111cell111; yhbatch -N 2 -n 24 -p TH_NET -J N111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Cu-N2/110cell111; yhbatch -N 2 -n 24 -p TH_NET -J N110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Cu-N2/100cell111; yhbatch -N 2 -n 24 -p TH_NET -J N100 RunVASP.5.4.1.sh', \
                        # 'cd /vol-th/home/Chgy/liyd/AuNT/zigzag/20/1/5.49/phonon/disp-008; yhbatch -N 2 -n 24 -p TH_NET -J 008 RunVASP.5.4.1.sh', \
                        # 'cd /vol-th/home/Chgy/liyd/AuNT/zigzag/20/1/5.49/phonon/disp-009; yhbatch -N 2 -n 24 -p TH_NET -J 009 RunVASP.5.4.1.sh', \
                        # 'cd /vol-th/home/Chgy/liyd/AuNT/zigzag/20/1/5.49/phonon/disp-011; yhbatch -N 2 -n 24 -p TH_NET -J 011 RunVASP.5.4.1.sh']

# my_jobs = ['cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/100; yhbatch -N 1 -n 12 -p TH_NET -J CO100 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/1002; yhbatch -N 1 -n 12 -p TH_NET -J CO1002 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/110; yhbatch -N 1 -n 12 -p TH_NET -J CO110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/1102; yhbatch -N 1 -n 12 -p TH_NET -J CO1102 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/111; yhbatch -N 1 -n 12 -p TH_NET -J CO111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/1112; yhbatch -N 1 -n 12 -p TH_NET -J CO1112 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/100; yhbatch -N 1 -n 12 -p TH_NET -J H2100 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/1002; yhbatch -N 1 -n 12 -p TH_NET -J H21002 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/110; yhbatch -N 1 -n 12 -p TH_NET -J H2110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/1102; yhbatch -N 1 -n 12 -p TH_NET -J H21102 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/111; yhbatch -N 1 -n 12 -p TH_NET -J H2111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/1112; yhbatch -N 1 -n 12 -p TH_NET -J H21112 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-N2/100; yhbatch -N 1 -n 12 -p TH_NET -J N2100 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-N2/1002; yhbatch -N 1 -n 12 -p TH_NET -J N21002 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-N2/110; yhbatch -N 1 -n 12 -p TH_NET -J N2110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-N2/1102; yhbatch -N 1 -n 12 -p TH_NET -J N21102 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-N2/111; yhbatch -N 1 -n 12 -p TH_NET -J N2111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-N2/1112; yhbatch -N 1 -n 12 -p TH_NET -J N21112 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-NO/100; yhbatch -N 1 -n 12 -p TH_NET -J NO100 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-NO/1002; yhbatch -N 1 -n 12 -p TH_NET -J NO1002 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-NO/110; yhbatch -N 1 -n 12 -p TH_NET -J NO110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-NO/1102; yhbatch -N 1 -n 12 -p TH_NET -J NO1102 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-NO/111; yhbatch -N 1 -n 12 -p TH_NET -J NO111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-NO/1112; yhbatch -N 1 -n 12 -p TH_NET -J NO1112 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/100; yhbatch -N 1 -n 12 -p TH_NET -J O2100 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/1002; yhbatch -N 1 -n 12 -p TH_NET -J O21002 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/110; yhbatch -N 1 -n 12 -p TH_NET -J O2110 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/1102; yhbatch -N 1 -n 12 -p TH_NET -J O21102 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/111; yhbatch -N 1 -n 12 -p TH_NET -J O2111 RunVASP.5.4.1.sh', \
                        # 'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/1112; yhbatch -N 1 -n 12 -p TH_NET -J O21112 RunVASP.5.4.1.sh', \
                        # ]

                        
my_jobs = [\
                        'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/1112; yhbatch -N 1 -n 12 -p TH_NET -J O21112 RunVASP.5.4.1.sh', \
                        'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-O2/111; yhbatch -N 1 -n 12 -p TH_NET -J O2111 RunVASP.5.4.1.sh', \
                        'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-CO/111; yhbatch -N 1 -n 12 -p TH_NET -J CO111 RunVASP.5.4.1.sh', \
                        'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/111; yhbatch -N 1 -n 12 -p TH_NET -J H2111 RunVASP.5.4.1.sh', \
                        'cd  /vol-th/home/Chgy/liyd/MJ_TH_test/Pd/Pd-H2/110; yhbatch -N 1 -n 12 -p TH_NET -J H2110 RunVASP.5.4.1.sh', \
                        ]
          
# phonon_root = '/vol-th/home/Chgy/liyd/AuNT/zigzag/20/1/5.49/phonon/'
# for i in range(3, 21):
    # job = ['cd '+phonon_root+'disp-%03d; yhbatch -N 2 -n 24 -p TH_NET -J %03d RunVASP.5.4.1.sh'%(i, i)]
    # my_jobs.extend(job)
logfile = open('th.log', 'w')    
logfile.close()
max_jobs = 30
ssh = paramiko.client.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('TH-1A-LN2', username='Chgy', password='@DjojW0Lt@()')
while True:
    time.sleep(30)
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command('yhqueue')
    # print ssh_stdout.read()
    curr_jobs = ssh_stdout.readlines()
    curr_jobs = curr_jobs[1: ]
    print len(curr_jobs), len(my_jobs)
    if len(curr_jobs) == max_jobs:
        continue
    for i in range(max_jobs-len(curr_jobs)):
        if len(my_jobs) == 0:
            break
        a, b, c = ssh.exec_command(my_jobs[0])
        info = 'std out: ' + b.read().rstrip() + '\nstd err: ' + c.read().rstrip()+'\n'
        print info
        logfile = open('th.log', 'a')    
        logfile.write(my_jobs[0]+'\n')
        logfile.write(info+'\n')
        logfile.close()
        my_jobs = my_jobs[1: ]
    