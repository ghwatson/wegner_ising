import subprocess
import numpy as np
import multiprocessing
import os
import sys

#Parameter space for jobs
e_i = 0.
e_f = 0.5
e_steps = 6
tf_i = 1.
tf_f = 2.
tf_steps = 11

L = 10

def worker(data):
    #unpack worker data.
    e,tf = data;

    #call the c++ program that will perform the job
    script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    wegner_path = script_path + "/wegner {0} {1} {2}".format(L,e,tf)
    args = "bin/bar -c somefile.xml -d text.txt -r aString -f anotherString".split()
    cmd = ["/bin/bash",wegner_path]
    cmd = wegner_path
    subprocess.call(cmd, shell=True)


data = [(e,tf) for tf in np.linspace(tf_i,tf_f,tf_steps) for e in np.linspace(e_i,e_f,e_steps)]
print data

#worker([0.1,0.2])
# Calculate correlators using a pool of workers.
pool_size = multiprocessing.cpu_count() - 1
pool = multiprocessing.Pool(processes=pool_size)
#pool_outputs = pool.map(worker,data)
pool.map(worker, data)
pool.close()
pool.join()
