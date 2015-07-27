from scipy.integrate import cumtrapz
from scipy.special import ellipe,ellipk
from subprocess import call
from copy import deepcopy
from matplotlib.pyplot import *
from numpy import *
import pickle
import h5py

class Matrices():
    def __init__(self,dat):
        self.matrix =  dat['Matrix'][:]
        self.bcmatrix = dat['BC'][:]
        self.D = dat['D'][:]
        self.D2 = dat['D2'][:]

class Parameters():
    def __init__(self,dat,fname='results.hdf5'):
    	self.m= int(dat['m'])
    	self.nr= int(dat['nr'])
        self.nf= int(dat['nf'])
    	self.ri= float(dat['ri'])
    	self.ro= float(dat['ro'])
    	self.mdisk= float(dat['mdisk'])
    	self.rs= float(dat['rs'])
    	self.h0= float(dat['h0'])
    	self.sig_ind= float(dat['sig_ind'])
    	self.flare_ind= float(dat['flare_ind'])
    	self.alpha_s= float(dat['alpha_s'])
    	self.alpha_b= float(dat['alpha_b'])
    	self.np= int(dat['np'])
    	self.gam= float(dat['gam'])
    	self.beta= float(dat['beta'])
    	self.tol= float(dat['tol'])
    	self.Nplanets= int(dat['Nplanets'])
    	self.outputname= fname.split('.')[0]


class Globals():
    def __init__(self,dat):
        self.lr = dat['lr'][:]
        self.r = dat['r'][:]
        self.omega = dat['Omega'][:]
        self.c2 = dat['c2'][:]
        self.sigma = dat['Sigma'][:]
        self.hor = dat['hor'][:]
        self.pres = dat['P'][:]
        self.temp = dat['T'][:]
        self.kappa2 = dat['kappa2'][:]
        self.dldc2 = dat['dc2'][:]
        self.dlds = dat['ds'][:]
        self.dldpres = dat['dp'][:]
        self.dldtemp = dat['dt'][:]
        self.d2lds = dat['d2s'][:]
        self.d2ldpres = dat['d2p'][:]
        self.d2ldtemp = dat['d2t'][:]
        self.kappa = sqrt(self.kappa2)
        self.wp = (self.omega**2 - self.kappa2)/(self.omega+self.kappa)
        self.Q = self.kappa*sqrt(self.c2)/(pi*self.sigma)
        self.Q0 = self.omega*sqrt(self.c2)/(pi*self.sigma)
        self.mdisk = trapz(self.r*self.sigma,x=self.r)*2*pi

    def get2d(self,rlim=None,Nphi=200):
        if rlim != None:
            r = self.r[self.r<rlim]

        phi = linspace(0,2*pi,Nphi)
    	rr,pp = meshgrid(r,phi)
    	xx = rr*cos(pp)
    	yy = rr*sin(pp)
        return xx,yy

    def loglog(self,q):
        self.plot(q,True,True)
    def semilogx(self,q):
        self.plot(q,True,False)
    def semilogy(self,q):
        self.plot(q,False,True)
    def plot(self,q,logr=False,logy=False):
        arglist = vars(self).keys()
    	if q not in arglist:
    		print 'Bad argument'
    		return

    	dat = getattr(self,q)
        figure();
        if logr:
            if logy:
                loglog(self.r,dat)
            else:
                semilogx(self.r,dat)
        else:
            if logy:
                semilogy(self.r,dat)
            else:
                plot(self.r,dat)

        xlabel('$r$',fontsize=20)
        title(q,fontsize=20)



class Mode():
    def __init__(self,ev,evec,params,glbls):
            self.ev = ev
            self.u = evec[::params.nf]
            self.v = evec[1::params.nf]
            self.s = evec[2::params.nf]
            if params.nf==4:
                self.p = evec[3::params.nf]
            else:
                self.p = self.s

            self.m = params.m
            self.r = glbls.r
            self.e = self.u/(1j*self.m*self.r*(glbls.omega-self.ev))
            ilr_rhs = glbls.omega - glbls.kappa/self.m - self.ev.real
            olr_rhs = glbls.omega + glbls.kappa/self.m - self.ev.real
            cor_rhs = glbls.omega - self.ev.real
            self.ilr = self.r[sign(ilr_rhs)[1:] - sign(ilr_rhs)[:-1] != 0]
            self.olr = self.r[sign(olr_rhs)[1:] - sign(olr_rhs)[:-1] != 0]
            self.cor = self.r[sign(cor_rhs)[1:] - sign(cor_rhs)[:-1] != 0]
            overlap = sign(self.u.real[1:])-sign(self.u.real[:-1])
            self.nodes = len(overlap[overlap != 0])
            self.dbar = glbls.sigma
            self.vpbar = self.r*glbls.omega
            self.freq = self.m*(self.ev.real - glbls.omega)/glbls.kappa
            self.freq *= self.freq
            self.Qbarr = glbls.Q**2*(1-self.freq)


    def semilogx(self,points=True):
        self.plot(logxy=(True,False),points=points)
    def semilogy(self,points=True):
        self.plot(logxy=(False,True),points=points)
    def loglog(self,points=True):
        self.plot(logxy=(True,True),points=points)

    def plot(self,logxy=(False,False),points=True):
        fig,axes = subplots(2,2,sharex='row',figsize=(15,10))

        rline = '-k'
        iline = '--k'
        if points:
            rline = '-.k'
            iline = '--ok'

        if logxy == (False,False):
            axes[0,0].plot(self.r,self.u.real,rline,self.r,self.u.imag,iline)
            axes[0,1].plot(self.r,self.v.real,rline,self.r,self.v.imag,iline)
            axes[1,0].plot(self.r,self.s.real,rline,self.r,self.s.imag,iline)
#            axes[1,1].plot(self.r,self.p.real,'.k',self.r,self.p.imag,'ok')
            axes[1,1].plot(self.r,self.freq.real,'b',self.r,self.Qbarr,'m')
        elif logxy == (True,False):
            axes[0,0].semilogx(self.r,self.u.real,rline,self.r,self.u.imag,iline)
            axes[0,1].semilogx(self.r,self.v.real,rline,self.r,self.v.imag,iline)
            axes[1,0].semilogx(self.r,self.s.real,rline,self.r,self.s.imag,iline)
#            axes[1,1].semilogx(self.r,self.p.real,'.k',self.r,self.p.imag,'ok')
            axes[1,1].semilogx(self.r,self.freq.real,'b',self.r,self.Qbarr,'m')

        elif logxy == (False,True):
            axes[0,0].semilogy(self.r,self.u.real,rline,self.r,self.u.imag,iline)
            axes[0,1].semilogy(self.r,self.v.real,rline,self.r,self.v.imag,iline)
            axes[1,0].semilogy(self.r,self.s.real,rline,self.r,self.s.imag,iline)
#            axes[1,1].semilogy(self.r,self.p.real,'.k',self.r,self.p.imag,'ok')
            axes[1,1].semilogy(self.r,self.freq.real,'b',self.r,self.Qbarr,'m')
        elif logxy == (True,True):
            axes[0,0].loglog(self.r,self.u.real,rline,self.r,self.u.imag,iline)
            axes[0,1].loglog(self.r,self.v.real,rline,self.r,self.v.imag,iline)
            axes[1,0].loglog(self.r,self.s.real,rline,self.r,self.s.imag,iline)
#            axes[1,1].loglog(self.r,self.p.real,'.k',self.r,self.p.imag,'ok')
            axes[1,1].loglog(self.r,self.freq.real,'b',self.r,self.Qbarr,'m')

        axes[1,1].axhline(1,color='r',linestyle='--')
        axes[1,1].axhline(-1,color='r')
    #    axes[1,1].set_ylim(-1000,1000)
        axes[1,1].set_yscale('symlog')
        if len(self.ilr) != 0:
            for xin in self.ilr:
                for ax in axes.flatten():
                    ax.axvline(x=xin,color='r')
        if len(self.olr) != 0:
            for xin in self.olr:
                for ax in axes.flatten():
                    ax.axvline(x=xin,color='r',linestyle='--')
        if len(self.cor) != 0:
            for xin in self.cor:
                for ax in axes.flatten():
                    ax.axvline(x=xin,color='r',linewidth=2)

        axes[0,0].set_ylabel('u',fontsize=20)
        axes[0,1].set_ylabel('v',fontsize=20)
        axes[1,0].set_ylabel('$\\sigma$',fontsize=20)
        axes[1,1].set_ylabel('p',fontsize=20)
        axes[1,0].set_xlabel('r',fontsize=20)
        axes[1,1].set_xlabel('r',fontsize=20)

        axes[0,0].set_title('$\\Omega_p = %.2e + %.2ei$'%(self.ev.real,self.ev.imag),fontsize=20)

    def semilogxreal(self,Nphi=200,rmax=None,plotbar=True):
        self.plot(logx=True,logy=False,Nphi=200,rmax=None,plotbar=True)
    def semilogyreal(self,Nphi=200,rmax=None,plotbar=True):
        self.plot(logx=False,logy=True,Nphi=200,rmax=None,plotbar=True)
    def loglogreal(self,Nphi=200,rmax=None,plotbar=True):
        self.plot(logx=True,logy=True,Nphi=200,rmax=None,plotbar=True)
    def plotreal(self,logx=False,logy=False,Nphi=200,rmax=None,plotbar=True):
        phi = linspace(-pi,pi,Nphi)
        fig,axes = subplots(2,2,figsize=(8,6),dpi=80)

        if rmax != None:
            r = self.r[self.r < rmax]
            rf = self.r[self.r < rmax]
            sval = self.s[self.r<rmax]
        else:
            r = self.r
            rf = self.r
            sval = self.s
        if logx:
            r = log10(r)


        rr,pp = meshgrid(r,phi)


        x = rr*cos(pp)
        y = rr*sin(pp)

        rho = zeros(x.shape)
        vr = zeros(x.shape)
        vp = zeros(x.shape)
        p = zeros(x.shape)

        for i,r in enumerate(r):

            rho[:,i] =real(self.s[i]*exp(1j*self.m*phi))
            vr[:,i] = real(self.u[i]*exp(1j*self.m*phi))
            vp[:,i] = real(self.v[i]*exp(1j*self.m*phi))
            if plotbar:
                rho[:,i] += self.dbar[i]
                vp[:,i] += self.vpbar[i]

        if logy:
            rho = log10(rho)
            tstr = '$\\log_10 \\Sigma$'
        else:
            tstr = '$\\Sigma$'

        axes[0,0].pcolormesh(x,y,rho)
        axes[0,1].pcolormesh(x,y,vp)
        axes[1,0].pcolormesh(x,y,vr)
        axes[1,1]
        axes[0,0].set_title('$\\Sigma$')
        axes[0,1].set_title('$v_\\phi$')
        axes[1,0].set_title('$v_r$')


        if logx:
            if logy:
                axes[1,1].loglog(rf,sval.real,'-k',rf,sval.imag,'--k')
            else:
                axes[1,1].semilogx(rf,sval.real,'-k',rf,sval.imag,'--k')
        else:
            if logy:
                axes[1,1].semilogy(rf,sval.real,'-k',rf,sval.imag,'--k')
            else:
                axes[1,1].plot(rf,sval.real,'-k',rf,sval.imag,'--k')


        if len(self.ilr) != 0:
            for xin in self.ilr:
                axes[1,1].axvline(x=xin,color='r')
        if len(self.olr) != 0:
            for xin in self.olr:
                axes[1,1].axvline(x=xin,color='r',linestyle='--')
        if len(self.cor) != 0:
            for xin in self.cor:
                axes[1,1].axvline(x=xin,color='r',linewidth=2)

        axes[1,1].set_ylabel('$\\sigma$',fontsize=14)
        axes[1,1].set_xlabel('r',fontsize=14)
        axes[1,1].set_title('$\\Omega_p = %.2e + %.2ei$'%(self.ev.real,self.ev.imag),fontsize=14)


class Field():
    def __init__(self,fname='results.hdf5'):
        with h5py.File(fname,'r') as f:
            self.params = Parameters(f['Mateig/Parameters/Parameters'],fname)
            self.globals = Globals(f['Mateig/Globals'])
            self.matrices = Matrices(f['Mateig/Matrices'])
            self.evals = f['Mateig/Results/Evals'][:]
            self.evecs = f['Mateig/Results/Evecs'][:]

        self.inds = argsort(self.evals)
        self.sevals = self.evals[self.inds]
        self.modes = [None]*self.params.nr*self.params.nf
        for i,j in enumerate(self.inds):
            self.modes[i] = Mode(self.evals[j],self.evecs[j,:],self.params,self.globals)
        self.modes = array(self.modes)
        self.nodes = array([x.nodes for x in self.modes])
        self.modes = self.modes[argsort(self.nodes)]
    def argand(self,linrange=(1e-8,1e-8)):
        figure();
        plot(self.evals.real,self.evals.imag,'.')
        xscale('symlog',linthreshx=linrange[0])
        yscale('symlog',linthreshy=linrange[1])
        xlabel('Re($\\Omega_p$)',fontsize=20)
        ylabel('Im($\\Omega$)',fontsize=20)
