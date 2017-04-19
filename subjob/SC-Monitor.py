import paramiko
import time

TH1 = paramiko.client.SSHClient()
TH1.set_missing_host_key_policy(paramiko.AutoAddPolicy())
TH1.connect('TH-1A-LN2', username='Chgy', password='@DjojW0Lt@()')
while True:
    
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command('yhqueue')
    # print ssh_stdout.read()
    curr_jobs = ssh_stdout.readlines()
    curr_jobs = [line.split() for line in curr_jobs[1: ]]
    time.sleep(30)
    