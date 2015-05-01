from subprocess import call
from time import time

def param_dict(prof=None):
	params=dict()
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
			
	return params
	
def run_code(params, defines = None):
	tic = time()
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

	res = call(callstr)
	if res != 0:
		print '\n\nProblem with '
		print params
		print ' '.join(callstr)
		return -1
	else:
		fld = Field(params)
	
	toc = time()
#	print 'Final Running Time: %.4f seconds' % (tic - toc)
	return fld
	
def add_defines(defines_list):
	if type(defines_list) != list:
		defines_list = [defines_list]
	with open("defines.h","r+") as f:
		lines = f.readlines()
		for i,line in enumerate(lines):
			for val in defines_list:
				if val in line:
					if line[0:2] == '//':
						lines[i] = line[2:]
		f.seek(0,0)
		f.writelines(lines)
	return			
def remove_defines(defines_list):
	if type(defines_list) != list:
		defines_list = [defines_list]
	with open("defines.h","r+") as f:
		lines = f.readlines()
		
		for i,line in enumerate(lines):
			for val in defines_list:
				if val in line:
					if line[0:2] != '//':
						lines[i] = '//' + line
		f.seek(0,0)
		f.writelines(lines)
	return

def set_profile(prof):

	allowed_profs = [ 'EXPDECAY', 'POWER', 'MLIN','KUZMIN','RING','USER','INNERTAPER']
	
	if prof not in allowed_profs:
		print 'Not a valid profile name'
		print 'Choose from'
		print prof
	
	return
	


