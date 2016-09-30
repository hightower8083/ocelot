import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
import ocelot as oclt

def make_beam(bm):
	parts0 = np.zeros((6,bm['Np']))
	parts0[0] = bm['Rx']*np.random.randn(bm['Np'])
	parts0[1] = bm['Ox']*np.random.randn(bm['Np'])
	parts0[2] = bm['Ry']*np.random.randn(bm['Np'])
	parts0[3] = bm['Oy']*np.random.randn(bm['Np'])
	parts0[4] = bm['Lz']*np.random.randn(bm['Np'])
	parts0[5] = bm['dE']*np.random.rand(bm['Np'])
	p_array_init = oclt.ParticleArray()
	p_array_init.particles = parts0.T.flatten()
	p_array_init.E = bm['E']
	p_array_init.s = 0.0
	if 'Q' in bm:
		p_array_init.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array_init.q_array = np.ones(bm['Np'])
	return p_array_init

def make_shot(p_array_init, QAP_T, corr_qap):
	D1  = oclt.Drift(l=0.04685)
	D2  = oclt.Drift(l=0.10295)
	D3  = oclt.Drift(l=0.09895)
	D4 = oclt.Drift(l=0.52085)
	QAP1 = oclt.Quadrupole(l=0.047,  k1= 170.*corr_qap[0])
	QAP2 = oclt.Quadrupole(l=0.0511, k1=-170.*corr_qap[1])
	QAP3 = oclt.Quadrupole(l=0.0323, k1= 170.*corr_qap[2])
	cell = (D1, QAP1, D2, QAP2, D3, QAP3,D4)

	method = oclt.MethodTM()
	#method.global_method = oclt.SecondTM

	p_array = oclt.deepcopy(p_array_init)
	QAP1_Tz,QAP1_Tx,QAP2_Tz,QAP2_Tx,QAP3_Tz,QAP3_Tx = QAP_T

	lat = oclt.MagneticLattice(cell,method=method,start=D1,stop=D1)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	p_array.particles[::6] += QAP1_Tz*1e-3
	p_array.particles[2::6]+= QAP1_Tx*1e-3
	lat = oclt.MagneticLattice(cell,method=method,start=QAP1,stop=QAP1)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	p_array.particles[::6] -= QAP1_Tz*1e-3
	p_array.particles[2::6]-= QAP1_Tx*1e-3

	lat = oclt.MagneticLattice(cell,method=method,start=D2,stop=D2)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	p_array.particles[::6] += QAP2_Tz*1e-3
	p_array.particles[2::6]+= QAP2_Tx*1e-3
	lat = oclt.MagneticLattice(cell,method=method,start=QAP2,stop=QAP2)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	p_array.particles[::6] -= QAP2_Tz*1e-3
	p_array.particles[2::6]-= QAP2_Tx*1e-3

	lat = oclt.MagneticLattice(cell,method=method,start=D3,stop=D3)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	p_array.particles[::6] += QAP3_Tz*1e-3
	p_array.particles[2::6]+= QAP3_Tx*1e-3
	lat = oclt.MagneticLattice(cell,method=method,start=QAP3,stop=QAP3)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)
	p_array.particles[::6] -= QAP3_Tz*1e-3
	p_array.particles[2::6]-= QAP3_Tx*1e-3

	lat = oclt.MagneticLattice(cell,method=method,start=D4,stop=D4)
	navi = oclt.Navigator(lat)
	tw, p_array = oclt.track(lat,p_array,navi,print_progress=False)

	return p_array, np.array([p_array.x().mean()*1e3,p_array.y().mean()*1e3])

def test_bba(QAP_T,QAP_T_fixed,positions,corrections,p_array_init):
	err = 0.0
	for i in range(len(positions)):
		p_array, cents = make_shot(p_array_init, np.r_[QAP_T_fixed,QAP_T],corrections[i])
		err += (cents[0]-positions[i][0])**2+(cents[1]-positions[i][1])**2
	return err

def beam_plot(p_array_init,QAPshifts,corr_qap=1,ax=None,**kw):
	if ax==None: fig,ax = plt.subplots(1,1,**kw)
	p_array, cents = make_shot(p_array_init,QAPshifts.T.flatten(),[corr_qap,]*3)
	ax.hist2d(p_array.x()*1e3,p_array.y()*1e3,bins=400,range=[[-12,12],[-12,12]],cmap=plt.cm.viridis);
	print 'beam center is at:', cents
	if ax==None: plt.show()
