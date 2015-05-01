
from scipy.integrate import ode

def load_matrix():
	mat = loadtxt('matrix.dat')
	mat = mat[:,::2] + 1j*mat[:,1::2]
	mat *= -1j;
	mat[0,:] = 0; mat[-1,:] = 0
	return mat

def integ_func(t,y,mat):
	return dot(mat,y)

def advance_time(t_end, ei,dt,mat):
	nt = int(t_end/dt)
	t = linspace(0,t_end,nt)
	dt = diff(t)[0]
	print 'Using %d steps' % nt
	r = ode(integ_func).set_integrator('zvode', method='bdf')
	r.set_initial_value(ei, 0).set_f_params(mat)
	ef = zeros(ei.shape + (nt,),dtype='complex')
	ef[:,0] = ei
	for i,t1 in enumerate(t[1:]):
		while r.successful() and r.t < t1:
			r.integrate(r.t+dt)
			ef[:,i+1] = r.y
	return ef
	
	
def stepper(t,ei,mat,plot_results=False,r=None):
	bmat = eye(ei.shape[0])
	bmat[0,0]= -1; bmat[0,1] = 1;
	bmat[-1,-1] = 1; bmat[-1,-2] = -1;
	
	dt = diff(t)[0]
	
	ef = zeros(ei.shape + (len(t),),dtype='complex')
	ef[:,0] = ei
	err = zeros((len(t),),dtype='complex')
	err[0] = 0
	
		
	for i,t1 in enumerate(t[1:],start=1):
		print 'Advancing to time %f' % t1
		ef[:,i] = advance_time(dt, bmat, mat,ef[:,i-1])	
		err[i] = compute_err(ef[:,i-1],ef[:,i])
		
	if plot_results:
		fig,(axr,axi)=subplots(2,1,sharex='col');
		liner,=axr.semilogx(r,ei.real,'-k')
		linei,=axi.semilogx(r,ei.imag,'-k')
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
def advance_time(dt, bmat,mat,ei):
	lhs_matrix = bmat - dt*mat;
	rhs = dot(bmat,ei)
	ef =solve(lhs_matrix,rhs)
	return ef

def compute_err(ei,ef):
	return norm(ei-ef)