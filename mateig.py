from scipy.integrate import cumtrapz
from scipy.special import ellipe,ellipk
from subprocess import call
from copy import deepcopy
import pickle

class Mode():
	def __init__(self,ev,emode,(r,dlr,omega,sigma)):
		self.ev = ev
		self.emode = emode
		self.nodes = self.calc_nodes(emode.real)
		self.sig  = -gradient(sigma*emode,dlr)
		self.u = emode*r*(omega-ev)*1j
		self.v = 1j*self.u/2


	def calc_nodes(self,y):
		overlap = sign(y[1:]) - sign(y[:-1])
		return len(overlap[overlap != 0])

class Planet():
	def __init__(self,dat,pot):
		self.num = dat[0]
		self.a = dat[1]
		self.mp = dat[2]
		self.hill = dat[3]
		self.wp = dat[4]
		self.p0 = pot[0]
		self.p1 = pot[1]

class Field():
	def __init__(self,params):

		dat=loadtxt('globals.dat')
		emat=loadtxt('eigen.dat')
		self.params = deepcopy(params)
		self.defines = load_defines()
		self.matrix=loadtxt('matrix.dat')
		self.matrix = self.matrix[:,::2] + 1j*self.matrix[:,1::2]
		self.D = loadtxt('D.dat')
		self.D2 = loadtxt('D2.dat')


		self.lr = dat[:,0]
		self.r = dat[:,1]
		self.omega = dat[:,2]
		self.c2 = dat[:,3]
		self.sigma = dat[:,4]
		self.hor = dat[:,5]
		self.pres = dat[:,6]
		self.temp = dat[:,7]
		self.wp = dat[:,8]
		self.dldc2 = dat[:,9]
		self.dlds = dat[:,10]
		self.dldpres = dat[:,11]
		self.dldtemp = dat[:,12]
		self.d2lds = dat[:,13]
		self.d2ldpres = dat[:,14]
		self.d2ldtemp = dat[:,15]

		self.soft = self.params['rs'] * self.hor*self.r

		self.nu = self.params['alpha_s'] * self.c2/self.omega
		self.dlr = diff(self.lr)[0]
		self.nr = len(self.r)


		self.npl = self.params['Nplanets']
		if self.npl > 0:
			pdat = loadtxt('planet_summary.dat')
			ppot = loadtxt('planet_pots.dat')
			ppot = [ (ppot[:,4*i]+1j*ppot[:,4*i+1],ppot[:,4*i+2]+1j*ppot[:,4*i+3]) for i in range(self.npl)]
			if self.npl > 1:
				self.planets = [Planet(line,pp) for line,pp in zip(pdat,ppot)]
				self.plr = array([pl.a for pl in self.planets])
			else:
				self.planets = Planet(pdat,ppot[0])
				self.plr = self.planets.a

		evals = emat[:,0] + 1j*emat[:,1]
		emat = emat[:,2:]
		evecs = emat[:,::2] + 1j*emat[:,1::2]
		if self.npl > 0:
			self.pevecs = evecs[:,-self.npl:]
			evecs = evecs[:,:-self.npl]
			self.evecs = evecs
		# if npl != 0:
		# 	self.plevals = evals[-npl:]
		# 	evals = evals[:-npl]
		# 	evecs = evecs[:-npl,:-npl]
		# else:
		# 	self.plevals = 0





		inds = argsort(evals)

		self.Q = self.omega * sqrt(self.temp)/(pi*self.sigma)
		self.Mach = self.r*self.omega/sqrt(self.temp)
		self.mdisk = trapz(self.r*self.sigma,x=self.r)*2*pi

		(self.Veffp,self.Veffm,self.K2eff,self.Psi) = calculate_veff(self)

		self.evals = evals[inds]
		self.evecs = evecs[inds,:]

		self.evdict = {}
		self.edict= {}
		self.pedict = {}
		self.sigp = {}
		self.vrp = {}
		self.vpp = {}
		self.vortp = {}
		self.vortensp = {}
		self.modes = {}

		for i,ev in enumerate(self.evals):
			ee = self.evecs[i,:]
			self.modes[ev] = Mode(ev,ee,(self.r,self.dlr,self.omega,self.sigma))
			self.edict[ev] = ee
			if self.npl > 0:
				self.pedict[ev] = self.pevecs[i,:]
			self.sigp[ev] = -dot(self.D,self.sigma*ee)
			self.vrp[ev] = 1j*self.r*self.omega*ee
			self.vpp[ev] = .5*1j*self.vrp[ev]
			self.vortp[ev] =

#			if abs(ev) != 0:
#				self.evdict[self.nodes(self.evecs[i,:])] = ev
#				self.edict[self.nodes(self.evecs[i,:])] = self.evecs[i,:]
# 				self.sigp[ev]  = -gradient(self.sigma*self.evecs[i,:],self.dlr)
# 				self.sigp[self.nodes(self.evecs[i,:])] = self.sigp[ev]
# 				self.vrp[ev] = self.edict[ev]*self.r*(self.omega-ev)*1j
# 				self.vrp[self.nodes(self.evecs[i,:])] = self.vrp[ev]
# 				self.vpp[ev] = 1j*self.vrp[ev]/2
# 				self.vpp[self.nodes(self.evecs[i,:])] = self.vpp[ev]

#		self.edict = { ev:self.evecs[i,:] for i,ev in enumerate(self.evals)}
#		self.sigp = { ev:-gradient(self.sigma*self.evecs[i,:],self.dlr) for i,ev in enumerate(self.evals)}
#		self.vrp = { ev:self.edict[ev]*self.r*(self.omega-ev)*1j for ev in self.evals }
#		self.vpp = { ev:1j*self.vrp[ev]/2 for ev in self.evals }

#		self.bvort = self.kappa2/(2*self.omega)
#		self.vort = {ev:-.5*gradient(self.omega*self.edict[ev],self.dlr) for ev in self.evals}
#		self.vortens = {ev: (self.vort[ev]/self.sigma  - self.sigp[ev]*self.bvort/(self.sigma**2)) for ev in self.evals}

	def loglog(self,q):
		self.plot(q,True,True)
		return
	def semilogx(self,q):
		self.plot(q,True,False)
		return
	def semilogy(self,q):
		self.plot(q,False,True)

		return
	def plot(self,q,logr=False,logy=False):
		arglist = vars(self).keys()
		if logr:
			r = log10(self.r)
			xstr = '$\log_{10} r$'
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
				dat = log10(dat)
				tstr = '$\log_{10} ($' + q + ')'
			else:
				tstr = q

			figure()
			plot(r,dat)
			title(tstr,fontsize='large')
			xlabel(xstr,fontsize='large')


		return
	def omegap_plot(self,ev,logr=True):

		ilr = self.r[sign(self.wp-ev.real)[1:] - sign(self.wp-ev.real)[:-1] != 0]
		num_ilr = len(ilr)
		min_wp = self.wp.min()
		max_wp = self.wp.max()

		fig, ax1=subplots()
		ax1.set_xlabel('$r$',fontsize='large')
		ax1.set_title('$\\Omega_p = %.2e + %.2ei$' % (ev.real, ev.imag),fontsize='large')
		if logr:
			ax1.semilogx(self.r,self.wp,'-k')
			ax1.semilogx(self.r,ones(self.r.shape)*ev.real,'--')
			for lr in ilr:
				ax1.semilogx([lr,lr],[min_wp,max_wp],'--k',label='$r_{ilr}=%.2e$'%lr)


		else:
			ax1.plot(self.r,self.wp,'-k')
			ax1.plot(self.r,ones(self.r.shape)*ev.real,'--')
			for lr in ilr:
				ax1.plot([lr,lr],[min_wp,max_wp],'--k',label='$r_{ilr}=%.2e$'%lr)

		ax1.legend(loc='best')
		ax2 = ax1.twinx()

		if logr:
			ax2.semilogx(self.r, self.edict[ev].real,'-r')
		else:
			ax2.plot(self.r,self.edict[ev].real,'-r')


	def Veff_plot(self,ev,signk=1,logr=True,plot_lr=False):

		ilr = self.r[sign(self.wp-ev.real)[1:] - sign(self.wp-ev.real)[:-1] != 0]

		cor = self.r[sign(self.omega-ev.real)[1:] - sign(self.omega-ev.real)[:-1] != 0]

		for rc in cor:
			print 'Corotation at %.2e'%rc

		for lr in ilr:
			print 'Inner Linblad resonance at %.2e'%lr

		if signk == 1:
			Veff = self.Veffp
		else:
			Veff = self.Veffm


		num_ilr = len(ilr)
		min_veff = Veff.min()
		max_veff = Veff.max()

		tp = self.r[sign(Veff.real-ev.real)[1:] - sign(Veff.real-ev.real)[:-1] != 0]
		tpi = self.r[sign(Veff.imag-ev.imag)[1:] - sign(Veff.imag-ev.imag)[:-1] != 0]


		fig, (ax1,axk)=subplots(1,2,figsize=(20,10))
		ax1.set_xlabel('$r$',fontsize='large')
		ax1.set_ylabel('$V_{eff}$,    $\\Omega_p$',fontsize='large')
		ax1.set_title('$\\Omega_p = %.2e + %.2ei$' % (ev.real, ev.imag),fontsize='large')
		if logr:
			ax1.semilogx(self.r,Veff.real,'-k')
			ax1.semilogx(self.r,Veff.imag,'--k')




		else:
			ax1.plot(self.r,Veff.real,'-k')
			ax1.plot(self.r,Veff.imag,'--k')

		ax1.axhline(ev.real,color='b',linestyle='-')
		ax1.axhline(ev.imag,color='b',linestyle='--')
		if plot_lr:
			for lr in ilr:
				ax1.axvline(lr,color='k',linestyle='--',label='$r_{ilr}=%.2e$'%lr)
		for t in tp:
			ax1.axvline(t,color='r',linestyle='-',label='$r_{tp}=%.2e$'%t)
		for t in tpi:
			ax1.axvline(t,color='r',linestyle='--',label='$r_{tp}=%.2e$'%t)
		for rc in cor:
			ax1.axvline(rc,color='k',linewidth=2,label='$r_{co}=%.2e$'%rc)



		ax1.legend(loc='best')
		ax2 = ax1.twinx()

		if logr:
			ax2.semilogx(self.r, self.edict[ev].real,'-m')
			ax2.semilogx(self.r, self.edict[ev].imag,'--m')
		else:
			ax2.plot(self.r,self.edict[ev].real,'-m')
			ax2.plot(self.r,self.edict[ev].imag,'--m')

		axk.set_xlabel('$r$',fontsize='large')
		axk.set_title('$\\Omega_p = %.2e + %.2ei$' % (ev.real, ev.imag),fontsize='large')

		val = self.K2eff*(ev - Veff)
		if logr:
			axk.semilogx(self.r,val.real,'-g')
			axk.semilogx(self.r,val.imag,'--g')

		else:
			axk.plot(self.r,val.real,'-g')
			axk.plot(self.r,val.imag,'--g')


		if plot_lr:
			for lr in ilr:
				axk.axvline(lr,color='k',linestyle='--',label='$r_{ilr}=%.2e$'%lr)
		for t in tp:
			axk.axvline(t,color='r',linestyle='-',label='$r_{tp}=%.2e$'%t)
		for t in tpi:
			axk.axvline(t,color='r',linestyle='--',label='$r_{tp}=%.2e$'%t)
		for rc in cor:
			ax1.axvline(rc,color='k',linewidth=2,label='$r_{co}=%.2e$'%rc)


		axk.legend(loc='best')

	def plotreal(self,ev,q,total=False,logz=False, logx=False,logy=False,rlims=None,plot_contours=True):

		if rlims == None:
			r = self.r
			nr = self.nr
			ind = range(nr)
			sigma = self.sigma
			evec  = self.edict[ev]
		else:
			ind = self.r <= rlims
			r = self.r[ind]
			nr = len(r)
			sigma = self.sigma[ind]
			evec = self.edict[ev]
			evec = evec[ind]





		phi = linspace(0,2*pi,6*nr)
		rr,pp = meshgrid(r,phi)
		xx = rr*cos(pp)
		yy = rr*sin(pp)


		if xx.shape[0] == nr:
			ind = 0
		else:
			ind = 1



		sigp = -gradient(sigma*evec,self.dlr)
		ss = zeros(xx.shape)

		for i in range(nr):
			if ind == 0:
				ss[i,:] = real(sigp[i]*exp(1j*phi))
				if total:
					ss[i,:] += sigma[i]
			else:
				ss[:,i] = real(sigp[i]*exp(1j*phi))
				if total:
					ss[:,i] += sigma[i]




		figure();
		if logz:
			pcolormesh(xx,yy,log10(ss),cmap='hot'); colorbar()
			if plot_contours:
				contour(xx,yy,log10(ss))
		else:
			pcolormesh(xx,yy,ss,cmap='hot'); colorbar()
			if plot_contours:
				contour(xx,yy,ss)




		return


	def plotuvs(self,ev=None,node=None,logr=False,logy=False):

		if node == None:
			tstr = '$\\Omega_p$ = %.2e,\t$\\gamma$ = %.2e' % (real(ev), imag(ev))
		else:
			ev = self.evdict[node]
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


	def plotmode(self,ev=None,node=None,logr=False,logy=False,renormalize=False,scale=0):
		if logr:
			r = log10(self.r)
			xstr = '$\log_{10} r$'
			if self.npl > 0:
				plr = log10(self.plr)
		else:
			r = self.r
			xstr = '$r$'
			if self.npl > 0:
				plr = self.plr





		fig,(axex,axey,axe,axw) = subplots(4,1,sharex='col')
		axw.set_xlabel(xstr,fontsize='large')
		if logy:
			axe.set_ylabel('$ \log_{10} e(r)$',fontsize='large')
		else:
			axe.set_ylabel('$e(r)$',fontsize='large')
		axex.set_ylabel('$e_x(r)$',fontsize='large')
		axey.set_ylabel('$e_y(r)$',fontsize='large')
		axw.set_ylabel('$ | \omega(r) |/ \pi $',fontsize='large')
#		axw.set_ylim((-.9,1.1))


		if ev != None:
			if type(ev)==list or type(ev)==numpy.ndarray:
				keys = ev
			else:
				keys = [ev]

		elif node != None:
			if type(node)==list or type(node)==numpy.ndarray:
				keys = node
			else:
				keys = [node]

		else:
			print 'No mode specified, plotting the zero node mode'
			keys = [0]
			if self.edict[0] == 0:
				print 'There is no zero mode'
				return


		for x in keys:
			dat = copy(self.edict[x])

			if self.npl > 0:
				pdat = copy(self.pedict[x])
			if renormalize:
				dat /= self.sigma
			if scale != 0:
				if scale == 'max':
					dat /= dat.max()
					if self.npl > 0:
						pdat /= dat.max()
				else:
					dat *= scale/dat[0]
					if self.npl > 0:
						pdat *= scale/dat[0]

			axex.plot(r,real(dat))
			axey.plot(r,imag(dat))
			axe.plot(r,abs(dat))
			axw.plot(r,angle(self.edict[x])/pi)
			axex.set_title('$\\Omega_p = %.2e + %.2ei$' % (real(x),imag(x)))
			if self.npl > 0:
				axex.plot(plr,real(pdat),'-ro')
				axey.plot(plr,imag(pdat),'-ro')
				axe.plot(plr,abs(pdat),'-ro')
				axw.plot(plr,angle(self.pedict[x])/pi,'-ro')


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

#	def predicted_k(self,ev,sg=True,bt2=False,logr=False,logy=False):
#
#
# # 		if bt2:
# # 			kt2 = self.kappa2 -(self.omega - ev)**2
# # 		else:
# # 			kt2 = 2*self.omega*(self.kappa-self.omega - ev)
# #
# # 		kt2 /= self.c2
# #
# # 		if sg:
# # 			kc = pi*fld.sigma/fld.c2
# # 		else:
# # 			kc = 0
# #
# # 		kp = kc + sqrt(kc*kc + kt2)
# # 		km = kc - sqrt(kc*kc + kt2)
# #
# #
# # 		kpr = cumtrapz(kp,x=self.r,initial=0)
# # 		kmr = cumtrapz(km,x=self.r,initial=0)
# #
# #
#
# 		kc = pi*self.sigma/self.c2
# 		v = (ev.real - self.omega)/self.kappa
#
# 		k_short_trail = kc*(1 + sqrt( 1- (self.Q**2*(1-v**2)))
# 		k_long_trail = kc*(1 - sqrt( 1- (self.Q**2*(1-v**2)))
#
# 		k_short_lead = - k_short_trail
# 		k_long_lead = - k_long_trail
#
# 		kp = k_short_trail
# 		km = k_short_lead
#
# # 		a = self.c2/(self.r*self.r)
# # 		if sg:
# # 			b = -2*pi*self.sigma/self.r
# # 		else:
# # 			b = 0
# #
# # 		if bt2:
# # 			c = self.kappa**2 - (self.omega - ev)**2
# # 		else:
# # 			c = -2*self.omega*(self.kappa - self.omega - ev)
# #
# #
# # 		krp = -b + sqrt(b*b - 4*a*c)
# # 		krm = -b - sqrt(b*b - 4*a*c)
# # 		krp /= (2*a)
# # 		krm /= (2*a)
# #
# #
# #
#
# 		emode = self.edict[ev]
#
# 		dedr = gradient(emode,self.dlr)
# 		dabsedr = gradient(abs(emode),self.dlr)
#
# #		emodep = dabsedr*exp(1j*krp) + 1j*krp*emode
# #		emodem = dabsedr*exp(1j*krm) + 1j*krm*emode
#
#
#
# #		emodep = 1j*kp*emode*self.r
# #		emodem = 1j*km*emode*self.r
#
#
# #		emodep =  emode[0] * exp(1j*kp*(self.r-self.r[0]))
# #		emodem =  emode[0]*exp(1j*km*(self.r-self.r[0]))
#
# 		emodep = emode[0] * exp(1j*kpr)
# 		emodem = emode[0] * exp(1j*kmr)
# #		emodep = abs(emode) * exp(1j*kp*self.r)
# #		emodem = abs(emode) * exp(1j*km*self.r)
#
# 		figure()
# 		plot(self.r,kp,'-r',self.r,km,'-b',self.r,imag(kp),'--r',self.r,imag(km),'--b')
#
#
# 		figure()
# 		plot(self.r, emodep, '-r',self.r, emodem,'-b')
# 		plot(self.r,emode,'-k')
#
# # 		figure()
# # 		plot(self.r,real(emodep),'-r')
# # #		plot(self.r,imag(emodep),'--r')
# # 		plot(self.r,real(emodem),'-b')
# # #		plot(self.r,imag(emodem),'--b')
# # 		plot(self.r,real(dedr),'-k')
# # #		plot(self.r,imag(dedr),'-k')
# #
# 		if logy:
# 			yscale('symlog')
#
# 		if logr:
# 			xscale('symlog')
#
# 		return

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

		emode = copy(self.edict[ev])

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

	def argand(self,logx=False,logy=False,linrange=(1e-6,1e-6)):
	# Plots the Argand diagram for the pattern speed and growthrates

		figure()
		plot(real(self.evals),imag(self.evals),'x')

		xscale('symlog',linthreshx=linrange[0])
		yscale('symlog',linthreshy=linrange[1])
		xlim(-1e-1,1e-1)
		ylim(-1,1)
		grid('on')
		xlabel('$\\Omega_p$',fontsize='large')
		ylabel('$\\gamma$',fontsize='large')

	def nodes(self,ev):
		ex =(self.edict[ev]).real
		overlap = sign(ex[1:]) - sign(ex[:-1])
		return (len(overlap[overlap != 0]),ev)

	def sort_nodes(self,show_val=False):
		nodes_list = array([ self.nodes(x)[0] for x in self.evals])
		ev_list = array([ self.nodes(x)[1] for x in self.evals])
		ind=argsort(nodes_list)

		indx_list = [list(self.evals).index(x) for x in ev_list[ind]]
		sorted_node_list = nodes_list[ind]

		result = [ (nodes_list[i],ev_list[i]) for i in ind]
		sorted_res = [(sorted_node_list[i],indx_list[i]) for i in range(len(nodes_list))]
# 		if show_val==True:
# 			print result
# 		else:
# 			for x in sorted_res:
# 				print x

		return sorted_res
	def find_node(self,num):
		for x in self.evals:
			if self.nodes(x) == num:
				return x
		print 'Could not find the mode with %d nodes' % num
		return -1

	def generate_kernels(self):
		rp,rr = meshgrid(fld.r,fld.r)

		b = fld.params['rs'] * fld.params['h0'] * pow(rp,fld.params['flare_ind']+1)


		kval = 4*rr*rp/(b**2 + (rr+rp)**2)

		ee = ellipe(kval)
		ek = ellipk(kval)


		kern0 = 2*(b**2 - rr**2 + rp**2)*ee - 2*(b**2 + (rr-rp)**2)*ek
		kern0 /= (rr * (b**2 + (rr-rp)**2)*sqrt(b**2 + (rr+rp)**2))


		kern02 = 4*rr*(b**6 + b**4 * (-2*rr**2 + 3*rp**2) + (-rr**2 * rp + rp**3)**2 + b**2 * (-3*rr**4 - 4*rr**2 *rp**2 + 3*rp**4))*ee
		kern02 += -4*rr*(b**2 + (rr-rp)**2)*(b**4 + (rr**2-rp**2)**2 + b**2 * (rr**2 + 2*rp**2))*ek
		kern02 *= pow(b**2 + (rr-rp)**2,-2) * pow(b**2 + (rr+rp)**2,-1.5)


		kern1 = -2*(b**4 + 2*rr**4 - 3*rr**2 * rp**2 + rp**4 + b**2*(3*rr**2 + 2*rp**2))*ee
		kern1 += 2*(b**2 + (rr-rp)**2)*(b**2 + 2*rr**2 + rp**2)*ek
		kern1 /= ((b**2 + (rr-rp)**2)*rp*sqrt(b**2 + (rr+rp)**2))


		kern0dat = loadtxt('kernel0.dat')
		kern1dat = loadtxt('kernel.dat')
		kern02dat = loadtxt('kernel02.dat')

		err0 = abs( (kern0 - kern0dat)/kern0 )
		err1 = abs( (kern1 - kern1dat)/kern1 )
		err02 = abs( (kern02 - kern02dat)/kern02 )


		figure(); imshow(log10(err0)); colorbar(); title('Kernel 0')
		figure(); imshow(log10(err02)); colorbar(); title('Kernel 02')
		figure(); imshow(log10(err1)); colorbar(); title('Kernel 1')


		omg2 = zeros(fld.r.shape)
		kapg2 = zeros(fld.r.shape)
		phi1 = zeros(fld.r.shape)

		for i in range(fld.nr):
			omg2[i] = - trapz(fld.r**2*fld.sigma*kern0[i,:],x=fld.lr)/fld.r[i]
			kapg2[i] = - trapz(fld.r**2*fld.sigma*kern02[i,:],x=fld.lr)/(fld.r[i]**3)
#			phi1[i] = trapz(fld.r**2*fld.sigma*kern1[i,:],x=fld.lr)

		dat=loadtxt('omegacorrecs.dat')
		omg2dat = dat[:,3]
		kapg2dat = dat[:,-1]

		# fakedat = loadtxt('fakepotential.dat')
# 		phi1dat= fakedat[:,1]
#
# 		errph = abs((phi1-phi1dat)/phi1);

		errom = abs( (omg2 - omg2dat)/omg2)
		errkap = abs((kapg2 - kapg2dat)/kapg2)

		figure(); loglog(fld.r,errom,label='$ \\frac{ \\Delta \\Omega^2}{\\Omega^2}$')
		loglog(fld.r,errkap,label='$\\frac{ \\Delta \\kappa^2}{ \\kappa^2} $')
		xlabel('r',fontsize='large');
		legend(loc='best')

		figure(); semilogx(fld.r,omg2,'-',fld.r,omg2dat,'--');
		xlabel('r',fontsize='large');
		ylabel('$\\Omega^2$',fontsize='large');

		figure(); semilogx(fld.r,kapg2,'-',fld.r,kapg2dat,'--');
		xlabel('r',fontsize='large');
		ylabel('$\\kappa^2$',fontsize='large');

# 		figure(); loglog(fld.r,errph)
# 		xlabel('r',fontsize='large');
# 		title('$\\Phi$')
#
		return kern0, kern02, kern1, err0,err02,err1, omg2,kapg2, errom,errkap


def load_coeffs(r):
	dat = loadtxt('coeffs.dat')
	A = dat[:,0] + 1j*dat[:,1]
	B = dat[:,2] + 1j*dat[:,3]
	C = dat[:,4] + 1j*dat[:,5]

	fig,(ax1,ax2,ax3) = subplots(3,1,sharex='col')

	ax1.semilogx(r,A.real,'-b',r,A.imag,'-r')
	ax2.semilogx(r,B.real,'-b',r,B.imag,'-r')
	ax3.semilogx(r,C.real,'-b',r,C.imag,'-r')

	ax3.set_xlabel('$r$',fontsize='large')
	ax1.set_ylabel('$A$',fontsize='large')
	ax2.set_ylabel('$B$',fontsize='large')
	ax3.set_ylabel('$C$',fontsize='large')

	return A,B,C

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
	xlim(-1e-1,1e-1)
	ylim(-1,1)
	grid('on')
	return


def compare(flds, evnum, labels=None,logr=False,scale=0,logy=False,just_e=False):

	if type(evnum) != list and type(evnum) != array:
		evstr = ['$\\Omega_p$ = %.1e + %.1e i ' % (real(fld.evals[evnum]),imag(fld.evals[evnum])) for fld in flds]
	else:
		evstr = ['$\\Omega_p$ = %.1e + %.1e i ' % (real(fld.evals[evnum[i]]),imag(fld.evals[evnum[i]])) for i,fld in enumerate(flds)]


	if labels != None:
		fldstr = [labels[i] + ', ' + evstr[i] for i in range(len(flds))]
	else:
		fldstr = evstr

	if logr:
		r = [log10(fld.r) for fld in flds]
		xstr = '$\log_{10} r$'
	else:
		r = [fld.r for fld in flds]
		xstr = '$r$'


	if just_e:
		figure()
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

			if logr:
				if logy:
					loglog(fld.r,abs(dat,label=fldstr[i]))
				else:
					semilogx(fld.r,abs(dat),label=fldstr[i])
			else:
				if logy:
					semilogy(fld.r,abs(dat),label=fldstr[i])
				else:
					plot(fld.r,abs(dat),label=fldstr[i])
			ylabel('$|e|$', fontsize='large')
			xlabel('$r$', fontsize='large')

		if len(flds) > 3:
			legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=len(flds)/2, mode="expand", borderaxespad=0.)
		else:
			legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=len(flds), mode="expand", borderaxespad=0.)
	else:

		fig,(axex,axey,axe,axw) = subplots(4,1,sharex='col')
		axw.set_xlabel(xstr,fontsize='large')
		if logy:
			axe.set_ylabel('$ \log_{10} e(r)$',fontsize='large')
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
				axe.plot(r[i],log10(abs(dat)))
			else:
				axe.plot(r[i],abs(dat))
			axw.plot(r[i],angle(fld.edict[ev])/pi)

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

def lin_profile(q,r=linspace(.4,4,200)):
	r1=1
	r2=2
	epsilon=.1
	deltar = 5
	h=.05
	Qout= 2

	s = 1.5 + .5*q

	c = h*r**(-.5*q)
	c2 = c*c
	omk = r**(-1.5)
	H = c/omk


	delta1 = deltar*h*pow(r1,1.5-.5*q)
	delta2 = deltar*h*pow(r2,1.5-.5*q)

	bump = lambda r: (.5*(1-epsilon)*(1+tanh((r-r1)/delta1))+epsilon)*(.5*(1-epsilon)*(1-tanh((r-r2)/delta2))+epsilon)


#	sig_r2 = h*r2**(-.5*q) * pow(r2,-1.5) / (pi*Qout)

#	sigref = sig_r2 * pow(r2,s) / bump(r2)

	sigref = h * pow(r2,s - .5*q - 1.5)/(pi*Qout * bump(r2))

	sigma = sigref * bump(r) * pow(r,-s)

	Q = c*omk/(pi*sigma)

	sig_norm = sigref*bump(r1)


	fig,(axt,axm,axb) = subplots(3,1,sharex='col')

	axt.set_title('$q = %.1f$' % q)
	axt.plot(r,sigma/sig_norm)
	axb.set_xlabel('$r$')
	axt.set_ylabel('$\Sigma/\Sigma(R_0)$')
	axb.set_xlim(.4,4)

	axb.plot(r,Q)
	axb.set_ylabel('$Q$')
	axm.plot(r,c)
	axm.set_ylabel('$c^2$')

	figure()
	plot(r,bump(r))
	ylabel('Bump Function')
	xlabel('$r$')


	figure()
	plot(r,-1.5*omk/sigma)


	return r, sigma, c2, Q, H


def gamma_test(gamma_vals,flare_ind):


	h0 = .05
	rs = .1
	mdisk = .004
	np = 8
	sig_ind = -1.5
	alpha = 0
	ri = .4
	ro = 4
	nr = 512
	beta_cool = 0



	flds = dict()


	call(["./compile"])

	tot = 0
	for i in range(len(gamma_vals)):
		gam = gamma_vals[i]

		callstr = ['./a.out',str(nr),str(ri),str(ro),str(mdisk),str(rs),str(h0),str(sig_ind),str(flare_ind),str(alpha),str(np),str(gam),str(beta_cool)]
		print '\n\n\n'
		print tot, gam
		print callstr
		print '\n\n\n'
		tot += 1
		res=call(callstr)
		if res != 0:
			print '\n\nProblem with', (alpha,beta,flare)
			return alpha_vals,beta_vals,flare_vals,-1
		flds[gam] = Field()




	return gamma_vals, flds

def cooling_test(beta_vals):
	if type(beta_vals) == list:
		beta_vals = array(beta_vals)
	if type(beta_vals) == float:
		beta_vals = array([beta_vals])

	h0 = .05
	rs = .1
	mdisk = .04
	np = 8
	sig_ind = -1.5
	alpha = 0
	ri = .1
	ro = 10
	nr = 512
	gam = 2
	flare_ind = 0;



	nsgflds_low = dict()
	sgflds_hi = dict()

	nsgflds_low = dict()
	sgflds_hi = dict()




	call(["./compile"])

	tot = 0
	for i,beta in enumerate(beta_vals.round(2)):
		callstr = ['./a.out',str(nr),str(ri),str(ro),str(mdisk),str(rs),str(h0),str(sig_ind),str(flare_ind),str(alpha),str(np),str(gam),str(beta)]
		print '\n\n\n'
		print tot, beta
		print callstr
		print '\n\n\n'
		tot += 1
		res=call(callstr)
		if res != 0:
			print '\n\nProblem with', (alpha,beta,flare)
			return -1
		flds[beta] = Field()
	return beta_vals.round(2), flds


def predicted_rate(fld,params):
	mu = params['sig_ind']
	delta = 2*params['flare_ind'] - 1
	beta = params['beta']
	gamma = params['gam']
	l = params['ro'] - params['ri']
	k = pi/l


	if gamma != 1:
		fac = 1.+1j*beta/(1.+ beta**2)


		C = fld.temp*(gamma -1)/(2 * fld.omega * fld.r**2)

		A = fld.omega - fld.kappa + C*(2*mu + fac*(2 + mu + delta)*delta)

		B = C * (2 + mu + fac*(gamma * delta + (gamma -1) *( mu + 2)))

		C = C*( 1 + fac*(gamma -1))
		dedr = delta*beta * ( beta - 1j*gamma)/float((gamma**2 + beta**2))


	else:
		print 'Isothermal'
		C = fld.c2/(2*fld.omega*fld.r*fld.r)
		A = fld.omega - fld.kappa + C*mu*(2 + delta)
		B = C*(2 + delta + mu)
		dedr = 0

#	omp = A[-1] + dedr * (B[-1] + C[-1]*dedr)
#	omp = A[-1] + B[-1] * dedr + C[-1] * -1

	omp = A[-1] + (B[-1] + C[-1])*dedr - C[-1]*k*k*fld.r[-1]**2

	return omp

def fit_zn_power_law(fld):
	emode = copy(abs(fld.edict[fld.evals[-3]]))
	r  = copy(fld.r)
	lr = log10(r)

	lower_bound = lr[0] + .25*(lr[-1] - lr[0])
	upper_bound = lr[-1] - .25*(lr[-1] - lr[0])

	ind = (lr >= lower_bound) & (lr <= upper_bound)

	s = polyfit(lr[ind],log10(emode[ind]),1)

	figure()
	loglog(r,emode,'--k')
	loglog(r[ind],pow(10,s[0]*lr[ind] + s[1]),'-k',label='$ e \sim r^{%.4f}$' % s[0])
	xlabel('r',fontsize='large')
	ylabel('|e|',fontsize='large')
	legend(loc='best')

	return s

def predicted_omegap(params,eos):
	h = params['h0']
	mu = params['sig_ind']
	f = params['flare_ind']
	ro = params['ro']
	ri = params['ri']
	beta = params['beta']
	gamma = params['gam']
	d = 2*f -1
	k2 = (pi/log10(ro/ri))**2



	if eos == 'Barotropic':
		return .5*h**2*(mu - k2)*pow(ro,.5*(4*f -3))
	elif eos == 'Isothermal':
		return .5*h**2*(2*mu - k2 - f*(2*f-1+mu))*pow(ro,.5*(4*f-3))
	elif eos == 'Cooling':
		return gamma*.5*h**2*pow(ro,.5*(4*f-3)) * (-2*f*(mu + d)+2*mu+d*(2+mu+d)-gamma*k2)
	else:
		return 0

def kuzmin_veff(r,gam):
	return 3*(r**4 * ( 5 - 4*gam + 3*gam*gam)+r**2 * (6-8*gam) + 1)/(8*sqrt(r)*(1+r**2)**(.5*(3*gam+1)))

def load_defines():
	with open('src/defines.h') as f:
		lines = f.readlines()
	return [line.split()[1] for line in lines]

def calculate_veff(fld,signk=1):
	T = fld.temp
	om = fld.omega
	r = fld.r
	S = fld.sigma

	d = fld.dldtemp
	m = fld.dlds
	mp = fld.d2lds
	dp = fld.d2lds
	g = fld.params['gam']
	if 'SELFGRAVITY' in fld.defines:
		moq = fld.Mach/fld.Q
	else:
		moq = 0

	if 'BAROTROPIC' in fld.defines:
		veffp = fld.wp- (T/(8*om*r*r))*( (1+d-m)*(3+d-m)+2*(dp-mp))
		veffm = veffp +  (T/(8*om*r*r))*4*moq*(moq - 3*1j)
		veffp = veffp +  (T/(8*om*r*r))*4*moq*(moq + 3*1j)
		k2eff = -2*om/T
		psi = (sqrt(S/(r**3 * T)) + exp(-1j*signk*moq)) * S
	elif 'ISOTHERMAL' in fld.defines:
		veffp = fld.wp-(T/(8*om*r*r))*((m-3)*(m-1)-2*m)
		veffm = veffp + (T/(8*om*r*r))*4*moq*(moq - (3-d)*1j)
		veffp = veffp + (T/(8*om*r*r))*4*moq*(moq + (3-d)*1j)

		k2eff = -2*om/T
		psi = (sqrt(S/r**3) + exp(-1j*signk*moq)) * S
	elif 'COOLING' in fld.defines:
		veffp = fld.wp - (T/(8*om*r*r))*(g*(1+d+m)*(3+d+m)-4*(2+d)*(d+m)+2*(g-2)*(dp+mp))
		veffm = veffp + (T/(8*om*r*r))*4*moq*(moq - 3*1j)
		veffp = veffp + (T/(8*om*r*r))*4*moq*(moq + 3*1j)
		k2eff = -2*om/(g*T)
		psi = 1/sqrt(g*r**3*S*T) + exp(-1j*signk*moq)
	else:
		return 'No good EOS, outputinng precession rate'
		veffp = fld.wp
		veffm = fld.wp
		k2eff = 2*fld.omega/fld.temp

	return veffp,veffm,k2eff,psi

def save_Field(obj, filename):
	print 'Saving Field to ' + filename
	with open(filename, 'wb') as output:
		pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_Field(filename):
	print 'Loading Field from ' + filename
	with open(filename,'rb') as output:
		return pickle.load(output)


def create_power_law_plots(fld,moq,soft):
	cut_moq = lambda fld,m: ([fld[(r,m)] for r in soft], ['$r_s = %.2e$'%r for r in soft])
	cut_rs = lambda fld,r: ([fld[(r,m)] for m in moq], ['$M/Q_{outer} = %.2e$'%m for m in moq])

	figr, axr = subplots(2,3,sharex='col',sharey='row',figsize=(20,15))
	for ax in axr[1,:]:
		ax.set_xlabel('$r$', fontsize='large')
	for ax in axr[:,0]:
		ax.set_ylabel('$|e|$', fontsize='large')

	for ax,r in zip(axr.flatten(),soft):
		(list_f,label_f) = cut_rs(fld,r)
		zn = [x.sort_nodes()[0][1] for x in list_f]
		for i,f in enumerate(list_f):
			lbl=label_f[i] + ', $\\Omega_p = %.2e$'%f.evals[zn[i]].real
			ax.loglog(f.r,abs(f.edict[f.evals[zn[i]]]),label=lbl)
		ax.legend(loc='best')
		ax.set_title('$r_s = %.2e$' % r)
		ax.set_ylim(1e-3,1)

	subplots_adjust(wspace=0)

	figr, axm = subplots(2,3,sharex='col',sharey='row',figsize=(20,15))
	for ax in axm[1,:]:
		ax.set_xlabel('$r$', fontsize='large')
	for ax in axm[:,0]:
		ax.set_ylabel('$|e|$', fontsize='large')
	for ax,m in zip(axm.flatten(),moq):
		(list_f,label_f) = cut_moq(fld,m)
		zn = [x.sort_nodes()[0][1] for x in list_f]
		for i,f in enumerate(list_f):
			lbl=label_f[i] + ', $\\Omega_p = %.2e$'%f.evals[zn[i]].real
			ax.loglog(f.r,abs(f.edict[f.evals[zn[i]]]),label=lbl)
		ax.legend(loc='best')
		ax.set_title('$M/Q = %.2e$' % m)
		ax.set_ylim(1e-3,1)
	subplots_adjust(wspace=0)


def create_power_law_plots_2(fld,fld1,moq,soft,tstr=None):
	cut_moq = lambda fld,m: ([fld[(r,m)] for r in soft], ['$r_s = %.2e$'%r for r in soft])
	cut_rs = lambda fld,r: ([fld[(r,m)] for m in moq], ['$M/Q_{outer} = %.2e$'%m for m in moq])

	figr, axr = subplots(2,3,sharex='col',sharey='row',figsize=(20,15))
	for ax in axr[1,:]:
		ax.set_xlabel('$r$', fontsize='large')
	for ax in axr[:,0]:
		ax.set_ylabel('$|e|$', fontsize='large')

	for ax,r in zip(axr.flatten(),soft):
		(list_f,label_f) = cut_rs(fld,r)
		zn = [x.sort_nodes()[0][1] for x in list_f]
		for i,f in enumerate(list_f):
			lbl=label_f[i] #+ ', $\\Omega_p = %.2e$'%f.evals[zn[i]].real
			ax.loglog(f.r,abs(f.edict[f.evals[zn[i]]]),label=lbl)

		ax.set_color_cycle(None)
		(list_f,label_f) = cut_rs(fld1,r)
		zn = [x.sort_nodes()[0][1] for x in list_f]
		for i,f in enumerate(list_f):
			lbl=label_f[i] #+ ', $\\Omega_p = %.2e$'%f.evals[zn[i]].real
			ax.loglog(f.r,abs(f.edict[f.evals[zn[i]]]),'--')

		ax.legend(loc='best')
		if tstr != None:
			ax.set_title(tstr + ', '+'$r_s = %.2e$' % r)
		ax.set_ylim(1e-3,1)

	subplots_adjust(wspace=0)

	figr, axm = subplots(2,3,sharex='col',sharey='row',figsize=(20,15))
	for ax in axm[1,:]:
		ax.set_xlabel('$r$', fontsize='large')
	for ax in axm[:,0]:
		ax.set_ylabel('$|e|$', fontsize='large')
	for ax,m in zip(axm.flatten(),moq):
		(list_f,label_f) = cut_moq(fld,m)
		zn = [x.sort_nodes()[0][1] for x in list_f]
		for i,f in enumerate(list_f):
			lbl=label_f[i] #+ ', $\\Omega_p = %.2e$'%f.evals[zn[i]].real
			ax.loglog(f.r,abs(f.edict[f.evals[zn[i]]]),label=lbl)

		ax.set_color_cycle(None)
		(list_f,label_f) = cut_moq(fld1,m)
		zn = [x.sort_nodes()[0][1] for x in list_f]
		for i,f in enumerate(list_f):
			lbl=label_f[i] #+ ', $\\Omega_p = %.2e$'%f.evals[zn[i]].real
			ax.loglog(f.r,abs(f.edict[f.evals[zn[i]]]),'--')

		ax.legend(loc='best')
		if tstr != None:
			ax.set_title(tstr + ', '+'$M/Q = %.2e$' % m)
		ax.set_ylim(1e-3,1)
	subplots_adjust(wspace=0)


def sigma_from_machq_power_law(moq,h,ro,mu,delta):
	return h*h*moq*ro**(1+mu-delta) / pi

def sigma_from_machq_taper(moq,h,ro,rt,mu,delta):
	return sigma_from_machq_power(moq,h,ro,mu,delta)*(rt**3+rt**(mu))


def visc_runs(alpha,params,save_results = False,direct='saved_class/visc_diss/inner_taper/adiabatic/'):

# 	alpha =  array([ 0.    ,  0.0001,  0.001 ,  0.01  ,  0.1   ])
# params = {'Nplanets': 0, \
#  'alpha_b': 0, \
#  'alpha_s': 0.10000000000000001, \
#  'beta': 0, \
#  'flare_ind': 0.25, \
#  'gam': 2.3, \
#  'h0': 0.05, \
#  'mdisk': 0.0001, \
#  'np': 8, \
#  'nr': 400, \
#  'ri': 0.1, \
#  'ro': 10, \
#  'rs': 0.1, \
#  'sig_ind': -1.5, \
#  'tol': 1e-08}

	bcs = ['PRESBCIN','ZEROBCIN']
	rinner = [.01,.1]
	nrad = [600,400]

	fld = [[None for a in alpha] for j in range(4)]


	add_remove_option(bcs[0],bcs[1],True)


# [.01, 10]

	params['nr'] = nrad[0]
	params['ri'] = rinner[0]
	for i,a in enumerate(alpha):
		params['alpha_s'] = a
		if save_results:
			fname = direct + 'bcP_ri.01_alpha%d.dat'%i
			save_Field(run_code(params),fname)
			fld[0][i] = load_Field(fname)
		else:
			fld[0][i] = run_code(params)


	add_remove_option(bcs[1],bcs[0],True)
	for i,a in enumerate(alpha):
		params['alpha_s'] = a
		if save_results:
			fname = direct + 'bcZ_ri.01_alpha%d.dat'%i
			save_Field(run_code(params),fname)
			fld[1][i] = load_Field(fname)
		else:
			fld[1][i] = run_code(params)

# [.1, 10]
	params['nr'] = nrad[1]
	params['ri'] = rinner[1]
	add_remove_option(bcs[0],bcs[1],True)

	for i,a in enumerate(alpha):
		params['alpha_s'] = a
		if save_results:
			fname = direct + 'bcP_ri.1_alpha%d.dat'%i
			save_Field(run_code(params),fname)
			fld[2][i] = load_Field(fname)
		else:
			fld[2][i] = run_code(params)

	add_remove_option(bcs[1],bcs[0],True)
	for i,a in enumerate(alpha):
		params['alpha_s'] = a
		if save_results:
			fname = direct + 'bcZ_ri.1_alpha%d.dat'%i
			save_Field(run_code(params),fname)
			fld[3][i] = load_Field(fname)
		else:
			fld[3][i] = run_code(params)



	return fld


def gamma_var(alpha,gammas,params,save_results = False,direct='saved_class/visc_diss/inner_taper/adiabatic/'):

# gammas =  [1.1, 1.4, 1.7, 2.3]

	fld=[[None for a in alpha] for g in gammas]

	for i,g in enumerate(gammas):
		call(['mkdir',direct+'g%.1f'%g])
		for j,a in enumerate(alpha):
			params['gam'] = g
			params['alpha_s'] = a
			if save_results:
				fname = direct +'g%.1f/bcP_ri0.1_alpha%d.dat'%(g,j)
				save_Field(run_code(params),fname)
				fld[i][j] = load_Field(fname)
			else:
				fld[i][j] = run_code(params)

	return fld
def visc_diss_summary(fld,alpha,ind=-3):

	zn0_ev = [ fld[0][i].evals[ind] for i in range(len(alpha)) ]
	zn0_e = [abs(fld[0][i].edict[ev]) for i,ev in enumerate(zn0_ev)]

	zn1_ev = [ fld[1][i].evals[ind] for i in range(len(alpha)) ]
	zn1_e = [abs(fld[1][i].edict[ev]) for i,ev in enumerate(zn1_ev)]

	zn2_ev = [ fld[2][i].evals[ind] for i in range(len(alpha)) ]
	zn2_e = [abs(fld[2][i].edict[ev]) for i,ev in enumerate(zn2_ev)]

	zn3_ev = [ fld[3][i].evals[ind] for i in range(len(alpha)) ]
	zn3_e = [abs(fld[3][i].edict[ev]) for i,ev in enumerate(zn3_ev)]

	print zn2_ev

	in_pres_fld = fld[0]
	in_z_fld = fld[1]
	pres_fld = fld[2]
	z_fld = fld[3]


	fig,((ax0,ax2),(ax1,ax3)) = subplots(2,2,sharex='col')

	ax1.set_xlabel('$r$',fontsize='large')
	ax3.set_xlabel('$r$',fontsize='large')

	ax0.set_ylabel('$|e|$',fontsize='large')
	ax2.set_ylabel('$|e|$',fontsize='large')

	ax0.set_title('$\\delta P(r_i) = 0$, $r_i = 0.01$',fontsize='large')
	ax2.set_title('$\\delta P(r_i) = 0$, $r_i = 0.1$',fontsize='large')
	ax1.set_title('$e(r_i) = 0$, $r_i = 0.01$',fontsize='large')
	ax3.set_title('$e(r_i) = 0$, $r_i = 0.1$',fontsize='large')


	for i,a in enumerate(alpha):
		ax0.loglog(fld[0][i].r,zn0_e[i],label='$\\alpha = %.2e$, $\\gamma = %.2e$'%(a,imag(zn0_ev[i])))
		ax1.loglog(fld[1][i].r,zn1_e[i],label='$\\alpha = %.2e$, $\\gamma = %.2e$'%(a,imag(zn1_ev[i])))
		ax2.loglog(fld[2][i].r,zn2_e[i],label='$\\alpha = %.2e$, $\\gamma = %.2e$'%(a,imag(zn2_ev[i])))
		ax3.loglog(fld[3][i].r,zn3_e[i],label='$\\alpha = %.2e$, $\\gamma = %.2e$'%(a,imag(zn3_ev[i])))

	ax0.legend(loc='best')
	ax1.legend(loc='best')
	ax2.legend(loc='best')
	ax3.legend(loc='best')


	figure()

	loglog(alpha[1:],abs(imag(zn0_ev[1:])),label='$\\delta P(r_i) = 0$, $r_i = 0.01$')
	loglog(alpha[1:],abs(imag(zn1_ev[1:])),label='$e(r_i) = 0$, $r_i = 0.01$')
	loglog(alpha[1:],abs(imag(zn2_ev[1:])),label='$\\delta P(r_i) = 0$, $r_i = 0.1$')
	loglog(alpha[1:],abs(imag(zn3_ev[1:])),label='$e(r_i) = 0$, $r_i = 0.1$')

	xlabel('$\\alpha$',fontsize='large')
	ylabel('$|\\gamma|$',fontsize='large')
	legend(loc='best')
