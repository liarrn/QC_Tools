import paramiko
import time

k = paramiko.RSAKey.from_private_key_file("./nsfc2015_450.id")
c = paramiko.SSHClient()
c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
print "connecting"
c.connect( hostname = "172.16.22.11", username = "nsfc2015_450", pkey = k)

ssh_stdin, ssh_stdout, ssh_stderr = c.exec_command('ssh ln7; yhqueue')
curr_jobs = ssh_stdout.readlines()
print curr_jobs