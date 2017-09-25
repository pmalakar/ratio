import os
import sys
import time
from subprocess import *

nodes = [sys.argv[1]]

def runcmd (node, iter):

  script = './run.sh ' + str(iter) 
  cmd = 'sbatch -N '+ str(node) + ' -t ' + sys.argv[2] + ' ' + script
  print 'Executing ' + cmd
  jobid = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
  print 'Jobid : ' + jobid
  
  time.sleep(10)
  while True:
    #cmd = 'squeue --job ' + jobid.strip() + ' | grep preeti | awk \'{print $1}\''
    cmd = 'squeue | grep preeti | grep " ' + str(node)  + ' " | awk \'{print $1}\''
    jobrun = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    if jobrun == '':
      break
    time.sleep(30)

  return jobid.strip()

for iter in range (1, 3):
 for node in nodes:
    print '\nStarting on ' + str(node) + ' nodes' #+ str(rank) + ' ranks per node'
    jobid = runcmd(node, iter)

