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


    def semilogx(self):
        self.plot(logxy=(True,False))
    def semilogy(self):
        self.plot(logxy=(False,True))
    def loglog(self):
        self.plot(logxy=(True,True))

    def plot(self,logxy=(False,False)):
        fig,axes = subplots(2,2,sharex='row',figsize=(15,10))

        if logxy == (False,False):
            axes[0,0].plot(self.r,self.u.real,'-b',self.r,self.u.imag,'-r')
            axes[0,1].plot(self.r,self.v.real,'-b',self.r,self.v.imag,'-r')
            axes[1,0].plot(self.r,self.s.real,'-b',self.r,self.s.imag,'-r')
            axes[1,1].plot(self.r,self.p.real,'-b',self.r,self.p.imag,'-r')

        elif logxy == (True,False):
            axes[0,0].semilogx(self.r,self.u.real,'-b',self.r,self.u.imag,'-r')
            axes[0,1].semilogx(self.r,self.v.real,'-b',self.r,self.v.imag,'-r')
            axes[1,0].semilogx(self.r,self.s.real,'-b',self.r,self.s.imag,'-r')
            axes[1,1].semilogx(self.r,self.p.real,'-b',self.r,self.p.imag,'-r')

        elif logxy == (False,True):
            axes[0,0].semilogy(self.r,self.u.real,'-b',self.r,self.u.imag,'-r')
            axes[0,1].semilogy(self.r,self.v.real,'-b',self.r,self.v.imag,'-r')
            axes[1,0].semilogy(self.r,self.s.real,'-b',self.r,self.s.imag,'-r')
            axes[1,1].semilogy(self.r,self.p.real,'-b',self.r,self.p.imag,'-r')
        elif logxy == (True,True):
            axes[0,0].loglog(self.r,self.u.real,'-b',self.r,self.u.imag,'-r')
            axes[0,1].loglog(self.r,self.v.real,'-b',self.r,self.v.imag,'-r')
            axes[1,0].loglog(self.r,self.s.real,'-b',self.r,self.s.imag,'-r')
            axes[1,1].loglog(self.r,self.p.real,'-b',self.r,self.p.imag,'-r')



        axes[0,0].set_ylabel('u',fontsize=20)
        axes[0,1].set_ylabel('v',fontsize=20)
        axes[1,0].set_ylabel('$\\sigma$',fontsize=20)
        axes[1,1].set_ylabel('p',fontsize=20)
        axes[1,0].set_xlabel('r',fontsize=20)
        axes[1,1].set_xlabel('r',fontsize=20)

        axes[0,0].set_title('$\\Omega_p = %.2e + %.2ei$'%(self.ev.real,self.ev.imag),fontsize=20)

class Field():
    def __init__(self,fname='results.hdf5'):
        with h5py.File(fname,'r') as f:
            self.params = Parameters(f['Mateig/Parameters/Parameters'],fname)
            self.globals = Globals(f['Mateig/Globals'])
            self.matrices = Matrices(f['Mateig/Matrices'])
            self.evals = f['Mateig/Results/Evals'][:]
            self.evecs = f['Mateig/Results/Evecs'][:]

        self.inds = argsort(self.evals)
        self.modes = [None]*self.params.nr*self.params.nf
        for i,j in enumerate(self.inds):
            self.modes[i] = Mode(self.evals[j],self.evecs[j,:],self.params,self.globals)
        self.modes = array(self.modes)
