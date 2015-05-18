#!/usr/bin/env/python

from sys import argv
from subprocess import call


def create_defines_file():
	profiles = {'EXPDECAY_PROF':'src/profiles/expdecay_profile.c', \
				'MKLIN_PROF':'src/profiles/mklin_profile.c', \
				'POWERLAW_PROF':'src/profiles/powerlaw_profile.c', \
				'KUZMIN_PROF':'src/profiles/kuzmin.c', \
				'GAUSSIANRING_PROF':'src/profiles/gaussian_ring.c', \
				'GAUSSIANBUMP_PROF':'src/profiles/gaussian_bump.c', \
				'INNERTAPER_PROF':'src/profiles/innertaper_profile.c', \
				'PAPALOIZOU_PROF': 'src/profiles/papaloizou_profile.c'}
	
	with open('params.opt','r') as f:
		temp = [x.split('+') for x in f.readlines()]
		defs=[]
		for x in temp:
			if x[0] == '' and '#' not in x:
				def_str = x[-1].split('\n')[0]
				if '_PROF' in def_str:
					copy_profile( profiles[ def_str] )
					
				defs.append(def_str)
				
	with open('src/defines.h','w') as g:
		if defs != []:
			for x in defs:
				g.write('#define ' + x + '\n')
			if 'BAROTROPIC' not in defs:
				if 'POLYTROPE' in defs:
					g.write('#define BAROTROPIC\n')
					
		else:
			g.write('\n')
	
	
	
				
	return defs			

def copy_profile(fname):
	if type(fname) == str:
		print 'Copying ', fname, ' into src/profiles.c'
		call(['cp',fname,'src/profiles.c'])
	else:
		print fname, ' is not a valid file!'
	return

defs = create_defines_file()	



		
