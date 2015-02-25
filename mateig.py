from scipy.integrate import cumtrapz
from subprocess import call

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
		self.dldom = dat[:,11]
		self.d2dom = dat[:,12]
		self.nu = dat[:,13]
		self.dldnu = dat[:,14]
		self.dlr = diff(self.lr)[0]
		self.nr = len(self.r)
		evals = emat[:,0] + 1j*emat[:,1]
		emat = emat[:,2:]
		evecs = emat[:,::2] + 1j*emat[:,1::2]
				
		inds = argsort(evals)
		
		self.Q = self.kappa * self.omega/(pi*self.sigma)
		
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
		
	
	def plotmode(self,ev,logr=False,logy=False,renormalize=False,scale=0):
		if logr:
			r = self.lr
			xstr = '$\ln r$'
		else:
			r = self.r
			xstr = '$r$'
			
		
		fig,(axex,axey,axe,axw) = subplots(4,1,sharex='col')
		axw.set_xlabel(xstr,fontsize='large')
		if logy:
			axe.set_ylabel('$ \ln e(r)$',fontsize='large')
		else:
			axe.set_ylabel('$e(r)$',fontsize='large')
		axex.set_ylabel('$e_x(r)$',fontsize='large')
		axey.set_ylabel('$e_y(r)$',fontsize='large')
		axw.set_ylabel('$ | \omega(r) |/ \pi $',fontsize='large')
#		axw.set_ylim((-.9,1.1))
		
		
#	 	fig,(axe,axw) = subplots(2,1,sharex='col')
# 		axw.set_xlabel(xstr,fontsize='large')
# 		axe.set_ylabel('$e(r)$',fontsize='large')
# 		axw.set_ylabel('$ | \omega(r) |/ \pi $',fontsize='large')
		
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
				
				axex.plot(r,real(dat))
				axey.plot(r,imag(dat))
				axe.plot(r,abs(dat))
				axw.plot(r,abs(angle(self.edict[x])/pi))
		else:
			dat = copy(self.edict[ev])
			if renormalize:
					dat /= self.sigma
			if scale != 0:
				if scale == 'max':
					dat /= dat.max()
				else:
					dat *= scale/dat[0]
			
			axex.plot(r,real(dat),'-k',label=('$\\Omega_p = %.2e$' % real(ev)))
			axey.plot(r,imag(dat),'-k',label=('$\\gamma = %.2e$' % imag(ev)))
			axex.legend(loc='best')
			axey.legend(loc='best')
			axe.plot(r,abs(dat),'-k')
			axw.plot(r,abs(angle(self.edict[ev])/pi),'-k')
		
		subplots_adjust(hspace=.1)
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
	
	
		if bt2:
			kt2 = self.kappa2 -(self.omega - ev)**2
		else:
			kt2 = 2*self.omega*(self.kappa-self.omega - ev)
			
		kt2 /= self.c2
		
		if sg:
			kc = pi*fld.sigma/fld.c2
		else:
			kc = 0
			
		kp = kc + sqrt(kc*kc + kt2)
		km = kc - sqrt(kc*kc + kt2)
		
		
		kpr = cumtrapz(kp,x=self.r,initial=0)
		kmr = cumtrapz(km,x=self.r,initial=0)
		
		
	
# 		a = self.c2/(self.r*self.r)
# 		if sg:
# 			b = -2*pi*self.sigma/self.r
# 		else:
# 			b = 0
# 			
# 		if bt2:
# 			c = self.kappa**2 - (self.omega - ev)**2
# 		else:
# 			c = -2*self.omega*(self.kappa - self.omega - ev)
# 			
# 	
# 		krp = -b + sqrt(b*b - 4*a*c)
# 		krm = -b - sqrt(b*b - 4*a*c)
# 		krp /= (2*a)
# 		krm /= (2*a)
# 		
# 				
# 		
		
		emode = self.edict[ev]
		
		dedr = gradient(emode,self.dlr)
		dabsedr = gradient(abs(emode),self.dlr)
		
#		emodep = dabsedr*exp(1j*krp) + 1j*krp*emode
#		emodem = dabsedr*exp(1j*krm) + 1j*krm*emode
		
		
		
#		emodep = 1j*kp*emode*self.r
#		emodem = 1j*km*emode*self.r


#		emodep =  emode[0] * exp(1j*kp*(self.r-self.r[0]))
#		emodem =  emode[0]*exp(1j*km*(self.r-self.r[0]))

		emodep = emode[0] * exp(1j*kpr)	
		emodem = emode[0] * exp(1j*kmr)	
#		emodep = abs(emode) * exp(1j*kp*self.r)
#		emodem = abs(emode) * exp(1j*km*self.r)
			
		figure()
		plot(self.r,kp,'-r',self.r,km,'-b',self.r,imag(kp),'--r',self.r,imag(km),'--b')
	
	
		figure()
		plot(self.r, emodep, '-r',self.r, emodem,'-b')
		plot(self.r,emode,'-k')
		
# 		figure()
# 		plot(self.r,real(emodep),'-r')
# #		plot(self.r,imag(emodep),'--r')
# 		plot(self.r,real(emodem),'-b')
# #		plot(self.r,imag(emodem),'--b')
# 		plot(self.r,real(dedr),'-k')
# #		plot(self.r,imag(dedr),'-k')
# 		
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

	def draw_ellipse(self,ev,num_ellipse,(xc,yc)=(0,0),Nph=500):
		pgrid = linspace(0,2*pi,Nph)
		# (Np , Nr)		
		
		emode = self.edict[ev]
		
		u = self.r*(self.omega - ev)*1j*emode
		v = .5*1j*u
		
		x = zeros((Nph,self.nr))
		y = zeros((Nph,self.nr))
		for i in range(Nph):
			# semi latus rectum
			l = self.r * (self.omega*self.r  + real(v*exp(-1j*pgrid[i])))
			l0 = self.r*self.r * self.omega
			p = self.r*(l/l0)*(l/l0)
			theta = pgrid[i] - angle(emode)
			a = p/(1-abs(emode)**2)
			b = p/sqrt(1-abs(emode)**2)
			x[i,:] = xc +a*cos(theta) 
			y[i,:] = yc + b*sin(theta)
			
			
		figure();	
		size_x = abs(x).max()
		size_y = abs(y).max()
		xlim((-size_x,size_x))
		ylim((-size_y,size_y))
		plot(0,0,'k*',markersize=10)
		for i in range(self.nr)[::self.nr/num_ellipse]:
			plot(x[:,i],y[:,i],'-k')
				
		return 
		
	def argand(self,logx=False,logy=False):
	# Plots the Argand diagram for the pattern speed and growthrates
	
		figure()
		plot(real(self.evals),imag(self.evals),'x')
		
		xscale('symlog',linthreshx=1e-8)
		yscale('symlog',linthreshy=1e-8)
		
		xlabel('$\\Omega_p$',fontsize='large')
		ylabel('$\\gamma$',fontsize='large')
		
		return


def argand_compare(flds,labelstr=None,tstr=None,linrange=(1e-6,1e-6)):
		
	if labelstr == None:
		labelstr = [ str(i) for i in range(len(flds)) ]
	

	figure()

	if tstr != None:
		title(tstr,fontsize='large')
		
			
	point_type = ['x', 'o', 'v', 'h', 's', 'D', '*', '^','+','p','8','<','>'] 
	for i,fld in enumerate(flds):
		plot(real(fld.evals),imag(fld.evals),point_type[i],label=labelstr[i])
	
	xscale('symlog',linthreshx=linrange[0])
	yscale('symlog',linthreshy=linrange[1])
	
	xlabel('$\\Omega_p$',fontsize='large')
	ylabel('$\\gamma$',fontsize='large')
	legend(loc='best')
	xlim(-1e-3,1e-3)
	ylim(-1,1)
	grid('on')
	return
	
	
def compare(flds, evnum, nustr=None,logr=False,scale=0,logy=False):
	
	if type(evnum) != list and type(evnum) != array: 
		evstr = ['$\\Omega_p$ = %.1e + %.1e i ' % (real(fld.evals[evnum]),imag(fld.evals[evnum])) for fld in flds]
	else:
		evstr = ['$\\Omega_p$ = %.1e + %.1e i ' % (real(fld.evals[evnum[i]]),imag(fld.evals[evnum[i]])) for i,fld in enumerate(flds)]
	
	
	if nustr != None:
		fldstr = [nustr[i] + ', ' + evstr[i] for i in range(len(flds))]
	else:
		fldstr = evstr
	
	if logr:
		r = [fld.lr for fld in flds]
		xstr = '$\ln r$'
	else:
		r = [fld.r for fld in flds]
		xstr = '$r$'
	
	
		
	fig,(axex,axey,axe,axw) = subplots(4,1,sharex='col')
	axw.set_xlabel(xstr,fontsize='large')
	if logy:
		axe.set_ylabel('$ \ln e(r)$',fontsize='large')
	else:
		axe.set_ylabel('$e(r)$',fontsize='large')
	axex.set_ylabel('$e_x(r)$',fontsize='large')
	axey.set_ylabel('$e_y(r)$',fontsize='large')
	axw.set_ylabel('$ | \omega(r) |/ \pi $',fontsize='large')
#	axw.set_ylim((-.9,1.1))
	for i,fld in enumerate(flds):
		if type(evnum) != list and type(evnum) != array:
			ev = fld.evals[evnum]
		else:
			ev = fld.evals[evnum[i]]
		dat = copy(fld.edict[ev])
		
		if scale != 0:
			if scale == 'max':
				dat /= dat.max()
			else:
				dat *= scale/dat[0]
		
		if fld == '':
			axex.plot(r[i],real(dat))
		else:
			axex.plot(r[i],real(dat),label=fldstr[i])
		axey.plot(r[i],imag(dat))
		if logy:
			axe.plot(r[i],log(abs(dat)))
		else:
			axe.plot(r[i],abs(dat))
		axw.plot(r[i],abs(angle(fld.edict[ev])/pi))
	

	if len(flds) > 3:
		axex.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=len(flds)/2, mode="expand", borderaxespad=0.)
	else:
		axex.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=len(flds), mode="expand", borderaxespad=0.)
	
	subplots_adjust(hspace=.1)
	return
	

def run_masses(mvals,alpha):
	nr = 200
	ri = 1
	ro = 100
	rs = .1 
	h0 = .05
	beta = -1.5
	flare = 0
	np = 8
	
	flds=dict()
	
	call(["./compile"])
	for i,mass in enumerate(mvals):
		print '\n\nWorking on Mdisk = ', mass
		callstr = ['./a.out',str(nr),str(ri),str(ro),str(mass),str(rs),str(h0),str(beta),str(flare),str(alpha),str(np)]
		res = call(callstr)
		if res != 0:
			print '\n\nProblem with ',mass
			return -1
		flds[mass] = Field()
	
	argand_compare([ flds[m] for m in mvals],labelstr=[ '$M_D = %.2e $' % m for m in mvals ],linrange=(1e-6,1e-9),tstr='$\\alpha=%.1e$' % alpha)
		
	return flds

def run_batch(ri,ro,nr):
	
	h0 = .05
	rs = .1
	mdisk = .004
	np = 8
	
	alpha_vals = array([0,1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
	beta_vals = array([ -1.5, -1, -.75,-.5, 0])
	flare_vals = array([0,.25,.5,.6])
	
	aa,bb = meshgrid(alpha_vals,beta_vals)
	
	flds = dict()
	
	
	
	(ni,nj) = aa.shape
	
	call(["./compile"])
	
	tot = 0
	for i in range(len(alpha_vals)):
		for j in range(len(beta_vals)):
			for k in range(len(flare_vals)):
				alpha = alpha_vals[i]
				beta = beta_vals[j]
				flare = flare_vals[k]
			
				callstr = ['./a.out',str(nr),str(ri),str(ro),str(mdisk),str(rs),str(h0),str(beta),str(flare),str(alpha),str(np)]
				print '\n\n\n'
				print tot, alpha, beta, flare
				print callstr
				print '\n\n\n'
				tot += 1
				res=call(callstr)
				if res != 0:
					print '\n\nProblem with', (alpha,beta,flare)
					return alpha_vals,beta_vals,flare_vals,-1
				flds[(alpha,beta,flare)] = Field()
	
	
	
	
	return alpha_vals, beta_vals, flare_vals, flds
	
	
def zeromode(flds,alpha_vals,beta_vals,flare_vals):
	
	
	na = len(alpha_vals)
	nb = len(beta_vals)
	nf = len(flare_vals)
	
	zm = []
	ev = []
	for beta in beta_vals:
		zm.append( [ flds[(a,beta,0)] for a in alpha_vals[1:] ] )
		ev.append( array([ flds[(a,beta,0)].evals[-2] for a in alpha_vals[1:] ]) )
	
	point_type = ['x', 'o', 'v', 'h', 's', 'D', '*', '^','+','p','8','<','>'] 
	linestyle = [ '-'+p for p in point_type]
	
	
	figure()
	for i in range(nb):
		loglog(alpha_vals[1:],abs(imag(ev[i])),point_type[i], label='$\\beta $= %.2f' % beta_vals[i])
	
	xlabel('$\\alpha$',fontsize='large')
	ylabel('|$\\gamma$|',fontsize='large')
	
#	xscale('symlog')
#	yscale('symlog')
	
	legend(loc='best')
	ylim(1e-12,1e-4)
	
#	print [ log10(abs(imag(y))) for y in ev ]
	
#	fits=[ polyfit(log10(aa[0,1:]),log10(abs(imag(y))),1) for y in ev ]
	
	
	zm=[]
	ev=[]
	for alpha in alpha_vals[1:]:
		zm.append([ flds[(alpha,beta,0)] for beta in beta_vals])
		ev.append( array([ flds[(alpha,beta,0)].evals[-2] for beta in beta_vals]))
	
	figure()
	for i,alpha in enumerate(alpha_vals[1:]):
		semilogy(beta_vals,abs(imag(ev[i])),point_type[i],label='$\\alpha$ = %.2e' % alpha)
	
	xlabel('$\\beta$',fontsize='large')
	ylabel('|$\\gamma$|',fontsize='large')
	legend(loc='best')
	
	
	
# alpha variation

	# Choose flare  = 0 slice


	ev_ab = [ array([ abs(imag(flds[ (alpha,beta,0) ].evals[-2])) for alpha in alpha_vals[1:]]) for beta in beta_vals]
	
	  

	# Choose beta = -1.5 slice

	ev_af = [ array([ abs(imag(flds[ (alpha,-1.5,flare) ].evals[-2])) for alpha in alpha_vals[1:]]) for flare in flare_vals]

# both 
# 		bb,ff = meshgrid(beta_vals,flare_vals)
# 		bbff=[ zip(bb[i,:],ff[i,:]) for i in range(bb.shape[0]) ]
# 		
# 		g_alpha = zeros(bb.shape)
# 		
# 		for i in range(bb.shape[0]):
# 			for j in range(bb.shape[1]):
# 				g_alpha[i,j] = abs(imag(flds[ (
		

# beta variation

	# Choose flare = 0 slice 
	ev_ba =  [ array([ abs(imag(flds[ (alpha,beta,0) ].evals[-2])) for beta in beta_vals]) for alpha in alpha_vals[1:]]
	
	# Choose alpha = 1e-3 slice
	ev_bf =  [ array([ abs(imag(flds[ (1e-3,beta,flare) ].evals[-2])) for beta in beta_vals]) for flare in flare_vals]
	

# delta variation
	
	# Choose beta = -1.5 slice
	ev_fa =  [ array([ abs(imag(flds[ (alpha,-1.5,flare) ].evals[-2])) for flare in flare_vals]) for alpha in alpha_vals[1:]]
	
	# Choose alpha = 1e-3 slice
	ev_fb =  [ array([ abs(imag(flds[ (1e-3,beta,flare) ].evals[-2])) for flare in flare_vals]) for beta in beta_vals]

	
	fig, ((ax_af,ax_bf,ax_fb),(ax_ab,ax_ba,ax_fa)) = subplots(2,3,sharex='col', sharey='row')
	
	ax_ab.set_xlabel('$\\alpha$',fontsize='large')
	ax_ba.set_xlabel('$\\beta$',fontsize='large')
	ax_fa.set_xlabel('$f$',fontsize='large')
	
	ax_af.set_ylabel('$|\\gamma|$',fontsize='large')
	ax_ab.set_ylabel('$|\\gamma|$',fontsize='large')
	
	ax_af.set_title('$\\beta = -1.5$',fontsize='large')
	ax_bf.set_title('$\\alpha = 1.0e-03$',fontsize='large')
	ax_fb.set_title('$\\alpha = 1.0e-03$',fontsize='large')
	
	ax_ab.set_title('$f = 0$',fontsize='large')
	ax_ba.set_title('$f  = 0$',fontsize='large')
	ax_fa.set_title('$\\beta = -1.5$',fontsize='large')
	
	
	for i,beta in enumerate(beta_vals):
		ax_ab.loglog(alpha_vals[1:],ev_ab[i],linestyle[i],label='$\\beta $= %.2f' % beta)
		ax_fb.semilogy(flare_vals, ev_fb[i],linestyle[i],label='$\\beta $= %.2f' % beta)
	
	for i, alpha in enumerate(alpha_vals[1:]):
		ax_ba.semilogy(beta_vals, ev_ba[i],linestyle[i],label='$\\alpha $= %.1e' % alpha)
		ax_fa.semilogy(flare_vals, ev_fa[i],linestyle[i],label='$\\alpha $= %.1e' % alpha)
	
	for i, flare in enumerate(flare_vals):
		ax_af.loglog(alpha_vals[1:],ev_af[i],linestyle[i], label='$f$ = %.2f' % flare)
		ax_bf.semilogy(beta_vals,ev_bf[i],linestyle[i], label='$f$ = %.2f' % flare)
	
	ax_af.legend(loc='best')
	ax_bf.legend(loc='best')
	ax_fb.legend(loc='best')
	
	ax_ab.legend(loc='best')
	ax_ba.legend(loc='best')
	ax_fa.legend(loc='best')
	
	ax_ab.set_xlim((.5*min(alpha_vals[1:]),1.5*max(alpha_vals[1:])))
	ax_af.set_xlim((.5*min(alpha_vals[1:]),1.5*max(alpha_vals[1:])))
	
	ax_ba.set_xlim(( min(beta_vals)-.5, max(beta_vals)+.5))
	ax_bf.set_xlim(( min(beta_vals)-.5, max(beta_vals)+.5))
	
	ax_fa.set_xlim(( min(flare_vals)-.5, max(flare_vals)+.5))
	ax_fb.set_xlim(( min(flare_vals)-.5, max(flare_vals)+.5))

	subplots_adjust(wspace=0,hspace=.1)

# 	color plots 

# fix alpha = 1e-3
		
	bb_fb,ff_fb = meshgrid(beta_vals,flare_vals)
	
	bbff=[ zip(bb_fb[i,:],ff_fb[i,:]) for i in range(bb_fb.shape[0]) ]
	
	g_fb = zeros(bb_fb.shape)
	
	for i in range(bb_fb.shape[0]):
		for j in range(bb_fb.shape[1]):
			g_fb[i,j] = abs(imag(flds[ (1e-3,bbff[i][j][0],bbff[i][j][1])].evals[-2]))


#	fix flare = 0
		
	bb_ba,aa_ba = meshgrid(beta_vals,alpha_vals[1:])
	
	bbaa=[ zip(bb_ba[i,:],aa_ba[i,:]) for i in range(bb_ba.shape[0]) ]
	
	g_ba = zeros(bb_ba.shape)
	
	for i in range(bb_ba.shape[0]):
		for j in range(bb_ba.shape[1]):
			g_ba[i,j] = abs(imag(flds[ (bbaa[i][j][1],bbaa[i][j][0],0)].evals[-2]))

# fix beta = -1.5

	ff_fa,aa_fa = meshgrid(flare_vals,alpha_vals[1:])
	
	ffaa=[ zip(ff_fa[i,:],aa_fa[i,:]) for i in range(ff_fa.shape[0]) ]
	
	g_fa = zeros(ff_fa.shape)
	
	for i in range(ff_fa.shape[0]):
		for j in range(ff_fa.shape[1]):
			g_fa[i,j] = abs(imag(flds[ (ffaa[i][j][1],-1.5,ffaa[i][j][0])].evals[-2]))


	
	fig1, (ax_fb,ax_ba,ax_fa) = subplots(1,3)
	
	ax_fb.set_title('$\\alpha$ = 1.0e-03, $\log_{10} | \\gamma |$',fontsize='large')
	ax_ba.set_title('$f$ =0, $\log_{10} | \\gamma |$',fontsize='large')
	ax_fa.set_title('$\\beta $= -1.5, $\log_{10} | \\gamma |$',fontsize='large')
	
	ax_fb.set_xlabel('$\\beta$',fontsize='large')
	ax_fb.set_ylabel('$f$',fontsize='large')
	
	ax_ba.set_xlabel('$\\beta$',fontsize='large')
	ax_ba.set_ylabel('$\log_{10} \\alpha$',fontsize='large')
	
	ax_fa.set_xlabel('$f$',fontsize='large')
	ax_fa.set_ylabel('$\log_{10} \\alpha$',fontsize='large')
	
	
	ax_fb.pcolormesh(bb_fb,ff_fb,log10(g_fb))
	ax_ba.pcolormesh(bb_ba,log10(aa_ba),log10(g_ba))
	ax_fa.pcolormesh(ff_fa,log10(aa_fa),log10(g_fa))

	
	
	
	
	return 
	
		
	
	
	
	
	
	
	
		
