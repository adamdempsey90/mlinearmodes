
class Field():
	def __init__(self):
		dat=loadtxt('globals.dat')
		emat=loadtxt('eigen.dat')
		self.matrix=loadtxt('matrix.dat')
		self.matrix = self.matrix[:,::2] + 1j*self.matrix[:,1::2]
		self.lr = dat[:,0]
		self.r = dat[:,1]
		self.omega = dat[:,2]
		self.c2 = dat[:,3]
		self.sigma = dat[:,4]
		self.hor = dat[:,5]
		self.soft = dat[:,6]
		self.dldc2 = dat[:,7]
		self.dlds = dat[:,8]
		self.kappa2 = dat[:,9]
		self.kappa = sqrt(self.kappa2)
		self.d2lds = dat[:,10]
		self.dlr = diff(self.lr)[0]
		self.nr = len(self.r)
		evals = emat[:,0] + 1j*emat[:,1]
		emat = emat[:,2:]
		evecs = emat[:,::2] + 1j*emat[:,1::2]
				
		inds = argsort(evals)
		
		self.evals = evals[inds]
		self.evecs = evecs[inds,:]
		self.edict = { ev:self.evecs[i,:] for i,ev in enumerate(self.evals)}
		self.sigp = { ev:-gradient(self.sigma*self.evecs[i,:],self.dlr) for i,ev in enumerate(self.evals)}
		self.vrp = { ev:self.edict[ev]*self.r*(self.omega-ev)*1j for ev in self.evals }
		self.vpp = { ev:1j*self.vrp[ev]/2 for ev in self.evals }
	
#		self.bvort = self.kappa2/(2*self.omega)
#		self.vort = {ev:-.5*gradient(self.omega*self.edict[ev],self.dlr) for ev in self.evals}
#		self.vortens = {ev: (self.vort[ev]/self.sigma  - self.sigp[ev]*self.bvort/(self.sigma**2)) for ev in self.evals}
	
	
	def plot(self,q,logr=False,logy=False):
		arglist = vars(self).keys()
		if logr:
			r = self.lr
			xstr = '$\ln r$'
		else:
			r = self.r
			xstr = '$r$'
		
		
		if q not in arglist:
			print 'Bad argument'
			return
		
		dat = getattr(self,q)
		

		if q=='matrix':
			ry,rx = meshgrid(r,r)
			figure()
			pcolormesh(rx,ry,real(dat))
			colorbar()
			xlabel(xstr,fontsize='large')
			ylabel(xstr,fontsize='large')
			title('Re('+q+')',fontsize='large')
			
			
			figure()
			pcolormesh(rx,ry,imag(dat))
			colorbar()
			xlabel(xstr,fontsize='large')
			ylabel(xstr,fontsize='large')
			title('Im('+q+')',fontsize='large')
			
		elif q=='evals':
		
			fig,ax=subplots()
			
			ax.set_xscale('symlog')
			ax.plot(real(self.evals),imag(self.evals),'x')
			ax.set_xlabel('$\\Omega_p$',fontsize='large')
			ax.set_ylabel('$\\gamma$',fontsize='large')
			ax.set_title('Eigenvalues',fontsize='large') 
			
			if imag(fld.evals).max() < 1e-2:
				ax.set_ylim(-1,1)
		else:
			if logy:
				dat = log(dat)
				tstr = '$\ln ($' + q + ')'
			else:
				tstr = q
				
			figure()
			plot(r,dat)
			title(tstr,fontsize='large')
			xlabel(xstr,fontsize='large')
		
		
		return
	
	def plotuvs(self,ev,logr=False,logy=False):
		
		tstr = '$\\Omega_p$ = %.2e,\t$\\gamma$ = %.2e' % (real(ev), imag(ev))
		fig,(axu,axv,axs) = subplots(3,sharex=True)
		
		axu.plot(self.r,real(self.vrp[ev]),'-k',label=r'Re(u)')
		axu.plot(self.r,imag(self.vrp[ev]),'--k',label=r'Im(u)')
		
		axv.plot(self.r,real(self.vpp[ev]),'-k',label=r'Re(v)')
		axv.plot(self.r,imag(self.vpp[ev]),'--k',label=r'Im(v)')
		
		axs.plot(self.r,real(self.sigp[ev]),'-k',label=r'Re($\sigma$)')
		axs.plot(self.r,imag(self.sigp[ev]),'--k',label=r'Im($\sigma$)')
		
		
		axu.set_ylabel('u',fontsize='large')
		axv.set_ylabel('v',fontsize='large')
		axs.set_ylabel('$\sigma$',fontsize='large')

		
		axu.set_title(tstr)
		axs.set_xlabel('r',fontsize='large')
		
		axu.legend(loc='best')
		axv.legend(loc='best')
		axs.legend(loc='best')
		
		
		
		if logr:
			axu.set_xscale('symlog')
			axv.set_xscale('symlog')
			axs.set_xscale('symlog')
		
		if logy:
			axu.set_yscale('symlog')
			axv.set_yscale('symlog')
			axs.set_yscale('symlog')
		
		fig.subplots_adjust(hspace=0)
		
		return
		
	
	def plotmode(self,ev,logr=False,renormalize=False,scale=0):
		if logr:
			r = self.lr
			xstr = '$\ln r$'
		else:
			r = self.r
			xstr = '$r$'
			
			
		figure()
		xlabel(xstr,fontsize='large')
		ylabel('$e(r)$',fontsize='large')
		if type(ev) == list or type(ev)==numpy.ndarray:
			for x in ev:
				dat = copy(self.edict[x])
				if renormalize:
					dat /= self.sigma
				if scale != 0:
					if scale == 'max':
						dat /= dat.max()
					else:
						dat *= scale/dat[0]
				
				plot(r,real(dat),'-k')
				plot(r,imag(dat),'--k')
		else:
			dat = copy(self.edict[ev])
			if renormalize:
					dat /= self.sigma
			if scale != 0:
				if scale == 'max':
					dat /= dat.max()
				else:
					dat *= scale/dat[0]
			
			plot(r,real(dat),'-k',label=('$\\Omega_p = %.2e$' % real(ev)))
			plot(r,imag(dat),'--k',label=('$\\gamma = %.2e$' % imag(ev)))
			legend(loc='best')
		
		return
	
	# def convert_real(self,ev,Nphi=500):
# 		phi = linspace(-pi,pi,Nphi)
# 		sigmatot = array([self.sigma[i] + 2*real(self.sigp[ev][i]*exp(phi)) for i in range(len(self.r))])
# 		
# 		rr,pp = meshgrid(phi,self.r)
# 		x = rr*cos(pp)
# 		y = rr*sin(pp)
# 		figure()
# 		pcolormesh(x,y,log10(sigmatot))
# 		colorbar()
# 		title('$\\Sigma$')
# 		
# 		return
	def predicted_k(self,ev,sg=True,bt2=False,logr=False,logy=False):
	
		a = self.c2/(self.r*self.r)
		if sg:
			b = -2*pi*self.sigma/self.r
		else:
			b = 0
			
		if bt2:
			c = self.kappa**2 - (self.omega - ev)**2
		else:
			c = -2*self.omega*(self.kappa - self.omega - ev)
			
	
		krp = -b + sqrt(b*b - 4*a*c)
		krm = -b - sqrt(b*b - 4*a*c)
		krp /= (2*a)
		krm /= (2*a)
		
				
		

		emode = self.edict[ev]
		
		dedr = gradient(emode,self.dlr)
		

		
		figure()
		plot(self.r,real(1j*krp*emode),'-r')
		plot(self.r,imag(1j*krp*emode),'--r')
		plot(self.r,real(1j*krm*emode),'-b')
		plot(self.r,imag(1j*krm*emode),'--b')
		plot(self.r,real(dedr),'-k')
		plot(self.r,imag(dedr),'-k')
		
		if logy:
			yscale('symlog')
		
		if logr:
			xscale('symlog')
			
		return
	
	def output_mode(self,ev,filename):
		with open(filename,'w') as f:
			for i in range(self.nr):
				line = (self.lr[i],self.r[i],real(self.vrp[ev][i]),imag(self.vrp[ev][i]),real(self.vpp[ev][i]),imag(self.vpp[ev][i]),real(self.sigp[ev][i]),imag(self.sigp[ev][i]))
				linestr = '%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n' % line
				f.write(linestr)
		
		return
			
		
		