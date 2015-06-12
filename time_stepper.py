
from scipy.integrate import ode

def compute_matrices(mat,ei):
	e,x = eig(mat)

	xinv = inv(x)
	yi = dot(xinv,ei)

	return [e,x,xinv,yi]


def evolve(tend,matrices):
	e = matrices[0]
	x = matrices[1]
	xinv = matrices[2]
	yi = matrices[3]

	yf = zeros(yi.shape,dtype='complex')

	for i,y in enumerate(yi):
		yf[i] = y * exp(e[i]*tend)

	ef = dot(x,yf)

	return ef


def time_stepper(t,mat,ei,plot_results=False,logr=True,r=None):

	matrices = compute_matrices(mat,ei)

	ef = zeros(ei.shape + (len(t),),dtype='complex')
	ef[:,0] = ei
	for i,a in enumerate(t[1:],start=1):
		ef[:,i] = evolve(t[i],matrices)

	if plot_results:
		fig,(axr,axi,axe) = subplots(3,1,sharex='col')

		axe.set_xlabel('$r$',fontsize='large')
		axi.set_ylabel('$e_y$',fontsize='large')
		axr.set_ylabel('$e_x$',fontsize='large')
		axe.set_ylabel('$|e|$',fontsize='large')
		axr.set_title('t = %.2f' % t[0])
		axr.set_ylim((ef.real.min(),ef.real.max()))
		axi.set_ylim((ef.imag.min(),ef.imag.max()))
		axe.set_ylim((abs(ef).min(),abs(ef).max()))
		if logr:
			liner,=axr.semilogx(r,ef[:,0].real)
			linei,=axi.semilogx(r,ef[:,0].imag)
			linee,=axe.semilogx(r,abs(ef[:,0]))
		else:
			liner,=axr.plot(r,ef[:,0].real)
			linei,=axi.plot(r,ef[:,0].imag)
			linee,=axe.plot(r,abs(ef[:,0]))

		for i,ti in enumerate(t[1:],start=1):
			liner.set_ydata(ef[:,i].real)
			linei.set_ydata(ef[:,i].imag)
			linee.set_ydata(abs(ef[:,i]))
			axr.set_title('t = %.2f' % ti)
			fig.canvas.draw()

	return ef

def load_matrix():
	mat = loadtxt('matrix.dat')
	mat = mat[:,::2] + 1j*mat[:,1::2]
	mat *= -1j;
#	mat[0,:] = 0; mat[-1,:] = 0
	return mat

def integ_func(t,y,mat):
	return dot(mat,y)
#
# def advance_time(t_end, ei,dt,mat):
# 	nt = int(t_end/dt)
# 	t = linspace(0,t_end,nt)
# 	dt = diff(t)[0]
# 	print 'Using %d steps' % nt
# 	r = ode(integ_func).set_integrator('zvode', method='bdf')
# 	r.set_initial_value(ei, 0).set_f_params(mat)
# 	ef = zeros(ei.shape + (nt,),dtype='complex')
# 	ef[:,0] = ei
# 	for i,t1 in enumerate(t[1:]):
# 		while r.successful() and r.t < t1:
# 			r.integrate(r.t+dt)
# 			ef[:,i+1] = r.y
# 	return ef
#

def stepper(t,ei,mat,plot_results=False,r=None,logr=False):

	mat[0,:] = 0; mat[-1,:]=0;

	lbmat = eye(ei.shape[0])
	rbmat = eye(ei.shape[0])

	lbmat[0,0]= -1; lbmat[0,1] = 1;
 	lbmat[-1,-1] = 1; lbmat[-1,-2] = -1;
#
 	rbmat = lbmat
#	rbmat[-1,-1] = 0; rbmat[0,0]= 0;

	dt = diff(t)[0]

	ef = zeros(ei.shape + (len(t),),dtype='complex')
	ef[:,0] = ei
	err = zeros((len(t),),dtype='complex')
	err[0] = 0


	for i,t1 in enumerate(t[1:],start=1):
		print 'Advancing to time %f' % t1
		ef[:,i] = advance_time(dt, lbmat, rbmat,mat,ef[:,i-1])
		err[i] = compute_err(ef[:,i-1],ef[:,i])

	if plot_results:
		fig,(axr,axi)=subplots(2,1,sharex='col');
		if logr:
			liner,=axr.semilogx(r,ei.real,'-k')
			linei,=axi.semilogx(r,ei.imag,'-k')
		else:
			liner,=axr.plot(r,ei.real,'-k')
			linei,=axi.plot(r,ei.imag,'-k')
		axi.set_xlabel('$r$',fontsize='large')
		axr.set_ylabel('$Re(e)$',fontsize='large')
		axi.set_ylabel('$Im(e)$',fontsize='large')
		axr.set_ylim((ef.real.min(),ef.real.max()))
		axi.set_ylim((ef.imag.min(),ef.imag.max()))
		for i,t1 in enumerate(t[1:],start=1):
			liner.set_ydata(ef[:,i].real)
			linei.set_ydata(ef[:,i].imag)
			axr.set_title('$t = %.1f$' % t1)
			fig.canvas.draw()


	return ef,err
def advance_time(dt, lbmat,rbmat,mat,ei):
	lhs_matrix = lbmat - dt*mat;
	rhs = dot(rbmat,ei)
	ef =solve(lhs_matrix,rhs)
	return ef

def compute_err(ei,ef):
	return norm(ei-ef)
