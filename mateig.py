
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
		self.kappa = dat[:,9]
		self.d2lds = dat[:,10]
		
		evals = emat[:,0] + 1j*emat[:,1]
		emat = emat[:,2:]
		evecs = emat[:,::2] + 1j*emat[:,1::2]
				
		inds = argsort(evals)
		
		self.evals = evals[inds]
		self.evecs = evecs[inds,:]
		self.edict = { ev:self.evecs[i,:] for i,ev in enumerate(self.evals)}
	
	def plot(self,q,logr=False):
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
			
		else:
			figure()
			plot(r,dat)
			title(q,fontsize='large')
			xlabel(xstr,fontsize='large')
		
		
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
					dat *= scale/dat[0]
				
				plot(r,real(dat),'-k')
				plot(r,imag(dat),'--k')
		else:
			dat = copy(self.edict[ev])
			if renormalize:
					dat /= self.sigma
			if scale != 0:
				dat *= scale/dat[0]
			
			plot(r,real(dat),'-k',label=('$\\Omega_p = %.2e$' % real(ev)))
			plot(r,imag(dat),'--k',label=('$\\gamma = %.2e$' % imag(ev)))
			legend(loc='best')
		
		return
		
			
		
		