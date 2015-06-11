#!/usr/bin/env/python

from sys import argv
from subprocess import call


def create_defines_file(optfile='params.opt',sourcedir='src/'):
	if sourcedir[-1] != '/':
		sourcedir += '/'

	profiles = {'EXPDECAY_PROF':'src/profiles/expdecay_profile.c', \
				'MKLIN_PROF':'src/profiles/mklin_profile.c', \
				'POWERLAW_PROF':'src/profiles/powerlaw_profile.c', \
				'KUZMIN_PROF':'src/profiles/kuzmin.c', \
				'GAUSSIANRING_PROF':'src/profiles/gaussian_ring.c', \
				'GAUSSIANBUMP_PROF':'src/profiles/gaussian_bump.c', \
				'INNERTAPER_PROF':'src/profiles/innertaper_profile.c', \
				'PAPALOIZOU_PROF': 'src/profiles/papaloizou_profile.c'}

	with open(optfile,'r') as f:
		temp = [x.split('+') for x in f.readlines()]
		defs=[]
		for x in temp:
			if x[0] == '' and '#' not in x:
				def_str = x[-1].split('\n')[0]
				if '_PROF' in def_str:
					copy_profile( profiles[def_str] ,sourcedir)

				defs.append(def_str)

	with open(sourcedir + 'defines.h','w') as g:
		if defs != []:
			for x in defs:
				g.write('#define ' + x + '\n')
			if 'BAROTROPIC' not in defs:
				if 'POLYTROPE' in defs:
					g.write('#define BAROTROPIC\n')

		else:
			g.write('\n')


	print 'Created the defines.h file:'
	call(['cat',sourcedir+'defines.h'])
	return defs

def copy_profile(fname,sourcedir='src/'):
	if type(fname) == str:
		print 'Copying ', fname, ' into ' + sourcedir+'profiles.c'
		call(['cp',fname,sourcedir + 'profiles.c'])
	else:
		print fname, ' is not a valid file!'
	return



if __name__ == "__main__":

	if len(argv) == 2:
		if str(argv[1])[-4:] == '.opt':
			optfile = str(argv[1])
			sourcedir = 'src/'
		else:
			sourcedir = str(argv[1])
			if sourcedir[-1] != '/':
				sourcedir += '/'
			optfile = 'params.opt'

	elif len(argv) == 3:

		if str(argv[1])[-4:] == '.opt':
			optfile = str(argv[1])
			sourcedir = str(argv[2])
		else:
			sourcedir = str(argv[1])
			if str(argv[2])[-4:] == '.opt':
				optfile = str(argv[2])
			else:
				print 'Could not find options file. Using default.'
				optfile = 'params.opt'

	else:
		print 'Using default options file and source directory'
		optfile = 'params.opt'
		sourcedir = 'src/'

	if sourcedir[-1] != '/':
		sourcedir += '/'

	print 'Using Source Directory as: %s\nUsing Options File: %s' %(sourcedir,optfile)

	defs = create_defines_file(optfile,sourcedir)
