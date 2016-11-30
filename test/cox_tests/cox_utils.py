import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
import ocelot as oclt
from cox_configs import *

def make_beam(bm):
	parts0 = np.zeros((6,bm['Np']))
	parts0[0] = bm['Rx']*np.random.randn(bm['Np'])
	parts0[1] = bm['Ox']*np.random.randn(bm['Np'])
	parts0[2] = bm['Ry']*np.random.randn(bm['Np'])
	parts0[3] = bm['Oy']*np.random.randn(bm['Np'])
	parts0[4] = bm['Lz']*np.random.randn(bm['Np'])
	parts0[5] = bm['dE']*np.random.randn(bm['Np'])
	p_array = oclt.ParticleArray()
	p_array.particles = parts0.T.flatten()
	p_array.E = bm['E']
	p_array.s = 0.0
	if 'Q' in bm:
		p_array.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array.q_array = np.ones(bm['Np'])
	return p_array

def make_beam_contin(bm):
	p_array = oclt.ParticleArray()

	E1 = bm['E'][0]
	E2 = bm['E'][1]
	dE = (E2-E1)/(E2+E1)
	p_array.E = 0.5*(E2+E1)

	parts0 = np.zeros((6,bm['Np']))
	parts0[0] = bm['Rx']*np.random.randn(bm['Np'])
	parts0[1] = bm['Ox']*np.random.randn(bm['Np'])
	parts0[2] = bm['Ry']*np.random.randn(bm['Np'])
	parts0[3] = bm['Oy']*np.random.randn(bm['Np'])
	parts0[4] = bm['Lz']*np.random.randn(bm['Np'])
	parts0[5] = dE*(2*np.random.rand(bm['Np'])-1)

	p_array.particles = parts0.T.flatten()
	p_array.s = 0.0
	if 'Q' in bm:
		p_array.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array.q_array = np.ones(bm['Np'])
	return p_array

def make_line(Drifts = Drifts,QuadLengths=QuadLengths,QuadGradients=QuadGradients,DipLengths=DipLengths,\
  DipAngles=DipAngles,UndulConfigs=UndulConfigs,cell_keys=cell_keys):
	DriftObjs = {}
	QuadObjs = {}
	DipObjs = {}
	LattObjs = {}
	Latt = []

	for key in Drifts.keys(): LattObjs[key] = oclt.Drift(l=Drifts[key])
	for key in QuadLengths.keys(): LattObjs[key] = oclt.Quadrupole(l=QuadLengths[key],k1=QuadGradients[key])
	for key in DipLengths.keys(): LattObjs[key] = oclt.RBend(l=DipLengths[key],angle=DipAngles[key])
	for key in  ['IMG1','IMG2','IMG4','IMG5',]: LattObjs[key] = oclt.Marker()
	LattObjs['UNDL'] = oclt.Undulator(Kx=UndulConfigs['Strength'],nperiods=UndulConfigs['NumPeriods'],lperiod=UndulConfigs['Period'])
	return LattObjs

def damp_particles(p_array, Rx,Ry):
	Np = p_array.size
	s0 = p_array.s
	Rx0, Ry0 = Rx(s0), Ry(s0)
	indx = np.nonzero((np.abs(p_array.x())<Rx0)*(np.abs(p_array.y())<Ry0))[0]
	p_array.particles = p_array.particles.reshape((p_array.size(),6))[indx,:].flatten()
	return p_array, p_array.size()

def make_shot(p_array_init,QAP_T,cell,meth=1):
	D1,QAP1,D2,QAP2,D3,QAP3,D4 = cell
	method = oclt.MethodTM()
	if meth>1:method.global_method = oclt.SecondTM

	p_array = oclt.deepcopy(p_array_init)
	QAP1_Tx,QAP1_Tz,QAP2_Tx,QAP2_Tz,QAP3_Tx,QAP3_Tz = QAP_T
	#print QAP1_Tx,QAP1_Tz,QAP2_Tx,QAP2_Tz,QAP3_Tx,QAP3_Tz

	lat = oclt.MagneticLattice(cell,method=method,start=D1,stop=D1)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	p_array.particles[::6] -= QAP1_Tx*1e-3
	p_array.particles[2::6]-= QAP1_Tz*1e-3
	lat = oclt.MagneticLattice(cell,method=method,start=QAP1,stop=QAP1)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	p_array.particles[::6] += QAP1_Tx*1e-3
	p_array.particles[2::6]+= QAP1_Tz*1e-3

	lat = oclt.MagneticLattice(cell,method=method,start=D2,stop=D2)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	p_array.particles[::6] -= QAP2_Tx*1e-3
	p_array.particles[2::6]-= QAP2_Tz*1e-3
	lat = oclt.MagneticLattice(cell,method=method,start=QAP2,stop=QAP2)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	p_array.particles[::6] += QAP2_Tx*1e-3
	p_array.particles[2::6]+= QAP2_Tz*1e-3

	lat = oclt.MagneticLattice(cell,method=method,start=D3,stop=D3)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	p_array.particles[::6] -= QAP3_Tx*1e-3
	p_array.particles[2::6]-= QAP3_Tz*1e-3
	lat = oclt.MagneticLattice(cell,method=method,start=QAP3,stop=QAP3)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	p_array.particles[::6] += QAP3_Tx*1e-3
	p_array.particles[2::6]+= QAP3_Tz*1e-3

	lat = oclt.MagneticLattice(cell,method=method,start=D4)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	return p_array, np.array([p_array.x().mean()*1e3,p_array.y().mean()*1e3])

def beam_plot(p_array_init,QAPshifts,cell,ax=None,meth=1,**kw):
	if ax==None: fig,ax = plt.subplots(1,1,**kw)
	p_array, cents = make_shot(p_array_init,QAPshifts.flatten(),cell,meth=meth)
	ax.hist2d(p_array.x()*1e3,p_array.y()*1e3,bins=400,range=[[-17,17],[-8,8]],cmap=plt.cm.viridis)
	print 'beam center is at:', cents
	if ax==None: plt.show()

def make_shot_aligned(p_array_init,cell,meth=1):
	method = oclt.MethodTM()
	if meth>1:method.global_method = oclt.SecondTM
	p_array = oclt.deepcopy(p_array_init)
	lat = oclt.MagneticLattice(cell,method=method)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	return p_array, np.array([p_array.x().mean()*1e3,p_array.y().mean()*1e3])

def beam_plot_aligned(p_array_init,cell,ax=None,meth=1,**kw):
	if ax==None: fig,ax = plt.subplots(1,1)
	p_array, cents = make_shot_aligned(p_array_init,cell,meth=meth)
	ax.hist2d(p_array.x()*1e3,p_array.y()*1e3,bins=400,cmap=plt.cm.viridis,**kw)
#	ax.hist2d(p_array.x()*1e3,p_array.y()*1e3,bins=400,range=[[-17,17],[-8,8]],cmap=plt.cm.viridis,**kw)
	if ax==None: plt.show()

