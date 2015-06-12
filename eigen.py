from subprocess import call
import os
from time import time
import config

def load_params(fname='params.in'):
	params={}


	int_keys = ['nr','Nplanets','np']
	str_keys = ['outputname']

	with open(fname,'r') as f:
		for line in f.readlines():
			if '#' not in line:
				sline = line.split('=')
				key = sline[0].strip()
				val = sline[-1].strip()
				if key in int_keys:
					params[key] = int(val)
				elif key in str_keys:
					params[key] = str(val)
				else:
					params[key] = float(val)


	return params

def load_defines():
	defs=[]
	with open('params.opt','r') as f:
		for line in f.readlines():
		    if '#' not in line and '+' in line:
		        defs.append((line.split('+')[1]).strip())
	return defs


def dump_params(params,fname='params.in'):
	skeys = ['nr','ri','ro','mdisk','rs','h0','sig_ind', \
			'flare_ind','alpha_s','alpha_b','np','gam', \
			'beta','tol','Nplanets','outputname']

	lines = [k + ' = ' + str(params[k]) for k in skeys]
	lines.insert(0,'# Input parameters for matrixeigenvalue code')
	with open(fname,"w") as f:
		f.write('\n'.join(lines))

	return

def param_dict(prof=None):
	params={}
	params['h0'] = .05
	params['rs'] = .1
	params['mdisk'] = .04
	params['np'] = 8
	params['sig_ind'] = -1.5
	params['alpha_s'] = 0
	params['alpha_b'] = 0
	params['ri'] = .4
	params['ro'] = 10
	params['nr'] = 400
	params['np'] = 8
	params['flare_ind'] = 0
	params['beta'] = 0
	params['gam'] = 2
	params['tol'] = 1e-8
	params['Nplanets'] = 0
	params['outputname'] = 'results'
	if prof=='yoram':
		params = {'alpha_b': 0, \
		 'alpha_s': 0, \
 		 'beta': 0, \
		 'flare_ind': 0, \
		 'gam': 2, \
		 'h0': 0.05, \
 		 'mdisk': 0.01, \
		 'np': 8, \
		 'nr': 400, \
		 'ri': 0.1, \
		 'ro': 10, \
		 'rs': 0.1, \
		 'sig_ind': -1.5, \
		 'tol': 1e-08}

	if prof=='tremaine':
		params['ri'] = .0009
#		params['ri'] = 0.00674
		params['ro'] = 148.41
		params['rs'] = .01
		params['mdisk'] = 1

	if prof in ['papa_a','papa_b']:
		params['ri'] = 1.1
		params['ro'] = 99.9
		if prof == 'papa_a':
			params['mdisk'] = .04
		else:
			params['mdisk'] = .004

	if prof == 'narrow ring':
		params['ri'] =  0.95122942
		params['ro'] =  1.0512711
		params['mdisk']= 1
		params['rs']=.01
		params['sig_ind'] = .01


	return params

def run_code(params, defines = None):

	if defines != None:
		add_defines(defines)
		call(['./compile'])

	callstr = ['./eigen']
	callstr.append(str(params['nr']))
	callstr.append(str(params['ri']))
	callstr.append(str(params['ro']))
	callstr.append(str(params['mdisk']))
	callstr.append(str(params['rs']))
	callstr.append(str(params['h0']))
	callstr.append(str(params['sig_ind']))
	callstr.append(str(params['flare_ind']))
	callstr.append(str(params['alpha_s']))
	callstr.append(str(params['alpha_b']))
	callstr.append(str(params['np']))
	callstr.append(str(params['gam']))
	callstr.append(str(params['beta']))
	callstr.append(str(params['tol']))
	callstr.append(str(params['Nplanets']))
	callstr.append(str(params['outputname']))

	res = call(callstr)
	if res != 0:
		print '\n\nProblem with '
		print params
		print ' '.join(callstr)
		return -1
	else:
		tic = time()
		fname = params['outputname'] + '.hdf5'
		print 'Loading from file ' + fname
		fld = Field(fname)
		toc = time()
		print 'Loading time: %.4f seconds' % (toc - tic)


	return fld

def add_option(opt,recompile=False):
	if type(opt) != list:
		opt = [opt]
	with open('params.opt','r+') as f:
		lines = f.readlines()
		f.seek(0)
		for i,line in enumerate(lines):
			for key in opt:
				if key in line:
					if '#' in line:
						lines[i] = lines[i].split('#')[-1]
		f.write(''.join(lines))

	if recompile:
		code_compile()

def remove_option(opt,recompile=False):
	if type(opt) != list:
		opt = [opt]
	with open('params.opt','r+') as f:
		lines = f.readlines()
		f.seek(0)
		for i,line in enumerate(lines):
			for key in opt:
				if key in line:
					if '#' not in line:
						lines[i] = '#' + lines[i]
		f.write(''.join(lines))

	if recompile:
		code_compile()

def add_remove_option(add_opts,rm_opts,recompile=False):
	if type(add_opts) != list:
		add_opts = [add_opts]
	if type(rm_opts) != list:
		rm_opts = [rm_opts]

	with open('params.opt','r+') as f:
		lines = f.readlines()
		f.seek(0)
		for i,line in enumerate(lines):
			for key in add_opts:
				if key in line:
					if '#' in line:
						lines[i] = lines[i].split('#')[-1]
			for key in rm_opts:
				if key in line:
					if '#' not in line:
						lines[i] = '#' + lines[i]
		f.write(''.join(lines))

	if recompile:
		code_compile()


def code_compile(optfile='params.opt',sourcedir='src/'):
	cdir = os.getcwd()

	config.create_defines_file(optfile,sourcedir)
	os.chdir(sourcedir)
	os.chdir('../')
	call(['make'])
	os.chdir(cdir)



def set_profile(prof):

	allowed_profs = [ 'EXPDECAY', 'POWER', 'MLIN','KUZMIN','RING','USER','INNERTAPER']

	if prof not in allowed_profs:
		print 'Not a valid profile name'
		print 'Choose from'
		print prof

	return
