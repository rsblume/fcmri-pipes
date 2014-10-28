#!/usr/bin/python

from subprocess import call

call("3dmaskdump -mask automask_d2_001_TS.1+orig -noijk cleanTS_001_TS.1_taskonly+orig > TaskTS_001_TS.1.txt",shell=True)
