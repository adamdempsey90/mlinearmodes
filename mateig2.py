from scipy.integrate import cumtrapz
from scipy.special import ellipe,ellipk
from subprocess import call
from copy import deepcopy
class Mode():
	def __init__(self,ev,emode,r,dmat,d2mat,omega,sigma):
		self.ev = ev
		self.r = r
		self.e = emode
		self.dre = dot(dmat,self.e)
		self.d2re = dot(d2mat,self.e)
		self.ex = emode.real
		self.ey = emode.imag
		self.nodes = count_nonzero(sign(self.ex[1:])-sign(self.ex[:-1]))
		self.sig  = -dot(dmat,sigma*self.e)
		self.u = 1j*r*omega*self.e
		self.v = .5*1j*self.u
		self.vort = .75 * omega * (self.e - 2 *self.dre)
		self.vortens =  self.vort/sigma - self.sig * omega / (2*sigma**2)


	def plot(self,logscale=True):

		fig, ((ax_ex,ax_sig),(ax_ey,ax_vort),(ax_e,ax_vortens)) = subplots(3,2,sharex='col')


		if logscale:
			ax_ex.semilogx(self.r,self.ex)
			ax_ey.semilogx(self.r,self.ey)
			ax_e.semilogx(self.r,abs(self.e))

			ax_sig.semilogx(self.r,self.sig)
			ax_vort.semilogx(self.r,self.vort)
			ax_vortens.semilogx(self.r,self.vortens)
		else:
			ax_ex.plot(self.r,self.ex)
			ax_ey.plot(self.r,self.ey)
			ax_e.plot(self.r,abs(self.e))

			ax_sig.plot(self.r,self.sig.real,self.r,self.sig.imag)
			ax_vort.plot(self.r,self.vort.real,self.r,self.vort.imag)
			ax_vortens.plot(self.r,self.vortens,self.r,self.vortens.imag)

		ax_e.set_xlabel('$r$',fontsize='large')
		ax_vortens.set_xlabel('$r$',fontsize='large')

		ax_ex.set_ylabel('$e_x$',fontsize='large')
		ax_ey.set_ylabel('$e_y$',fontsize='large')
		ax_e.set_ylabel('$|e|$',fontsize='large')

		ax_sig.set_ylabel('$\\sigma$',fontsize='large')
		ax_vort.set_ylabel('$\\omega_z$',fontsize='large')
		ax_vortens.set_ylabel('$\\xi_z$',fontsize='large')

		ax_ex.set_title('$\\Omega_p = %.2e + %.2ei$,\tN = %d' % (self.ev.real,self.ev.imag,self.nodes))



class Field():
	def __init__(self,params):
		dat=loadtxt('globals.dat')
		emat=loadtxt('eigen.dat')
		self.params = deepcopy(params)
		self.matrix=loadtxt('matrix.dat')
		self.matrix = self.matrix[:,::2] + 1j*self.matrix[:,1::2]
		self.dmat = loadtxt('D.dat')
		self.d2mat = loadtxt('D2.dat')

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

		self.Q = self.omega * sqrt(self.temp)/(pi*self.sigma)
		self.mdisk = trapz(self.r*self.sigma,x=self.r)*2*pi

		evals = emat[:,0] + 1j*emat[:,1]
		emat = emat[:,2:]
		evecs = emat[:,::2] + 1j*emat[:,1::2]

		inds = argsort(evals)



		self.evals = evals[inds]
		self.evecs = evecs[inds,:]

		self.modes = [Mode(self.evals[i],self.evecs[i,:],self.r,self.dmat,self.d2mat,self.omega,self.sigma) for i in range(len(self.evals))]
		self.sorted_inds = argsort([x.nodes for x in self.modes])
		self.node_list = [(self.modes[i].nodes, i) for i in self.sorted_inds]
