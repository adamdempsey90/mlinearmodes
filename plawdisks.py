from numpy import pi


def sigma_from_machq_power_law(moq,h,ro,mu,delta):
	return h*h*moq*ro**(1+mu-delta) / pi

def dump_params(params,fname='params.in'):
	skeys = ['nr','ri','ro','mdisk','rs','h0','sig_ind', \
			'flare_ind','alpha_s','alpha_b','np','gam', \
			'beta','tol','Nplanets','outputname']

	lines = [k + ' = ' + str(params[k]) for k in skeys]
	lines.insert(0,'# Input parameters for matrixeigenvalue code')
	with open(fname,"w") as f:
		f.write('\n'.join(lines))

	return

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




if __name__ == "__main__":
    mu_vals = [-1.75,-1.5,-1,0]
    f_vals = [0,.125,.25,.5]
    moq = [.003,.01,.03,.1,.3,1,3,10,30]


    params = load_params()



    flines = []

    for i,mu in enumerate(mu_vals):
        for j,f in enumerate(f_vals):
            for k,m in enumerate(moq):
                sig = sigma_from_machq_power_law(m,params['h0'],params['ro'],m,2*f-1)
                fname = 'pl'+str(i)+str(j)+str(k)
                params['sig_ind'] = mu
                params['flare_ind'] = f
                params['sig_0'] = sig
                params['outputname'] = fname
                fname += '.in'
                flines.append('\t'.join([str(mu),str(f),str(m),fname]))
                dump_params(params,fname)

    with open('keys.txt','w') as f:
        f.write('\n'.join(flines))
