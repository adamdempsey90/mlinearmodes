from subprocess import call

def param_dict():
	params=dict()
	params['h0'] = .05
	params['rs'] = .1
	params['mdisk'] = .04
	params['np'] = 8
	params['sig_ind'] = -1.5
	params['alpha'] = 0
	params['ri'] = .1
	params['ro'] = 10	
	params['nr'] = 512
	params['np'] = 8
	params['flare_ind'] = 0
	params['beta'] = 0
	params['gam'] = 2
	
	return params
	
def run_code(params, defines = None):
	if defines != None:
		add_defines(defines)
		call(['./compile'])
	
	callstr = ['./a.out']
	callstr.append(str(params['nr'])) 
	callstr.append(str(params['ri'])) 
	callstr.append(str(params['ro'])) 
	callstr.append(str(params['mdisk'])) 
	callstr.append(str(params['rs'])) 
	callstr.append(str(params['h0'])) 
	callstr.append(str(params['sig_ind'])) 
	callstr.append(str(params['flare_ind'])) 
	callstr.append(str(params['alpha'])) 
	callstr.append(str(params['np'])) 
	callstr.append(str(params['gam'])) 
	callstr.append(str(params['beta'])) 

	res = call(callstr)
	if res != 0:
		print '\n\nProblem with '
		print params
		return -1
	else:
		fld = Field()
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