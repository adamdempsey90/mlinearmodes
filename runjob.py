#!/usr/bin/python
import os
from subprocess import call
import eigen
import sys



argc = len(sys.argv)

if argc < 4:
    print "Not Enough Arguments"


dirname = str(sys.argv[1])

if dirname == '.' or dirname == './' or dirname == ' ':
    dirname = None


sourcedir = str(sys.argv[2])

optfile = str(sys.argv[3])

paramsfile = str(sys.argv[4])


if argc == 6:
    planetsfile = str(sys.argv[5])
else:
    planetsfile = 'None'

print 'Using the files:\n\tOptions File: %s\n\tParameter file: %s\n\tPlanets File: %s' %(optfile, paramsfile, planetsfile)
print 'Source directory: %s' % sourcedir



if dirname != None:
    print "Making directory " + dirname
    call(['mkdir',dirname])
    call(['cp', optfile, paramsfile,dirname])
    if planetsfile != 'None':
        call(['cp',planetsfile,dirname])

print 'Compiling code'

eigen.code_compile(optfile,sourcedir)

if dirname != None:
    call(['mv','eigen',dirname])
    os.chdir(dirname)

#call(['./eigen', '-i','params.in'])
