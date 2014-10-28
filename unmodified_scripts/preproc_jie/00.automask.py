#!/usr/bin/python

from subprocess import call

call("3dAutomask -dilate 2 -prefix automask_d2_001_TS.1 TS.1+orig[48]",shell=True)
