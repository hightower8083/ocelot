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
	parts0[5] = bm['dE']*np.random.randn(bm['Np'])
	p_array_init = oclt.ParticleArray()
	p_array_init.particles = parts0.T.flatten()
	p_array_init.E = bm['E']
	p_array_init.s = 0.0
	if 'Q' in bm:
		p_array_init.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array_init.q_array = np.ones(bm['Np'])
	return p_array_init

def make_cell(corr_qap=[1,1,1], grad_nom=[0.,0.,0.],\
	  drifts=[0.04685,0.10295,0.09895,0.52085],lengths=[0.047,0.0511,0.0323]):
	D1 = oclt.Drift(l=drifts[0])
	D2 = oclt.Drift(l=drifts[1])
	D3 = oclt.Drift(l=drifts[2])
	D4 = oclt.Drift(l=drifts[3])
	QAP1 = oclt.Quadrupole(l=lengths[0],  k1= grad_nom[0]*corr_qap[0])
	QAP2 = oclt.Quadrupole(l=lengths[1], k1= grad_nom[1]*corr_qap[1])
	QAP3 = oclt.Quadrupole(l=lengths[2], k1= grad_nom[2]*corr_qap[2])
	cell = [D1, QAP1, D2, QAP2, D3, QAP3,D4]
	return oclt.deepcopy(cell)

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

def test_bba(QAP_T,QAP_T_fixed,positions,corrections,p_array_init,grad_nom=[0.,0.,0.],meth=1):
	err = 0.0
	for i in range(len(positions)):
		cell = make_cell(corrections[i],grad_nom)
		p_array, cents = make_shot(p_array_init, np.r_[QAP_T_fixed,QAP_T],cell,meth=meth)
		err += (cents[0]-positions[i][0])**2+(cents[1]-positions[i][1])**2
	return err

def test_bba1(QAP_T13,QAP_T_fixed2,positions,corrections,p_array_init,grad_nom=[0.,0.,0.],meth=1):
	err = 0.0
	QAP_T13 = QAP_T13.reshape(2,2)
	for i in range(len(positions)):
		cell = make_cell(corrections[i],grad_nom)

		arr = np.concatenate((QAP_T13[0][None,:],QAP_T_fixed2[None,:],QAP_T13[1][None,:]),axis=0)
		p_array, cents = make_shot(p_array_init, arr.flatten(),cell,meth=meth)
		err += (cents[0]-positions[i][0])**2+(cents[1]-positions[i][1])**2
	return err

def test_bba2(QAP_T1,QAP_T_fixed23,positions,corrections,p_array_init,grad_nom=[0.,0.,0.],meth=1):
	err = 0.0
	for i in range(len(positions)):
		cell = make_cell(corrections[i],grad_nom)
		arr = np.concatenate((QAP_T1[None,:],QAP_T_fixed23),axis=0)
		p_array, cents = make_shot(p_array_init, arr.flatten(),cell,meth=meth)
		err += (cents[0]-positions[i][0])**2+(cents[1]-positions[i][1])**2
	return err

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

