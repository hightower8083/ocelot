import matplotlib.pyplot as plt
import numpy as np
import sys
import ocelot as oclt
from cox_configs import *

method = oclt.MethodTM()
method.global_method = oclt.SecondTM

def sliced_spectrum(Emin,Emax,dg = 0.001,):
	e_m = Emin
	e_p = Emin
	nrg_sliced = [e_m,]
	while e_p<=Emax:
		e_p = e_m*(1+dg)
		nrg_sliced.append(e_p)
		e_m = e_p
	return np.array(nrg_sliced)

def make_beam_sliced(bm,dg=0.002):
	p_arrays = []
	E1 = bm['E'][0]
	E2 = bm['E'][1]
	nrg_sliced =  sliced_spectrum(E1,E2,dg=dg)
	nrg_cntrs = 0.5*(nrg_sliced[1:]+nrg_sliced[:-1])
	Nbeams = nrg_sliced.shape[0]-1
	part_per_slice = np.round(bm['Np']*(nrg_sliced[1:]-nrg_sliced[:-1])/(nrg_sliced[-1]-nrg_sliced[0])).astype('i')
	for i in range(Nbeams):
		beam = oclt.deepcopy(bm)
		beam['Np'] = part_per_slice[i]
		beam['Q'] /= Nbeams
		beam['E'] = (nrg_sliced[i],nrg_sliced[i+1])
		p_arrays.append(make_beam_contin(beam))
	return p_arrays

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
  DipAngles=DipAngles,UndulConfigs=UndulConfigs,BeamEnergy=BeamEnergy_ref, BeamEnergy_ref=BeamEnergy_ref):
	DriftObjs = {}
	QuadObjs = {}
	DipObjs = {}
	LattObjs = {}
	Latt = []
	Ecorr = BeamEnergy_ref/BeamEnergy

	for key in Drifts.keys(): LattObjs[key] = oclt.Drift(l=Drifts[key])
	for key in QuadLengths.keys(): LattObjs[key] = oclt.Quadrupole(l=QuadLengths[key],k1=QuadGradients[key]*Ecorr)
	for key in DipLengths.keys(): LattObjs[key] = oclt.RBend(l=DipLengths[key],angle=DipAngles[key]*Ecorr)
	for key in  ['IMG1','IMG2','IMG4','IMG5',]: LattObjs[key] = oclt.Marker()
	LattObjs['UNDL'] = oclt.Undulator(Kx=UndulConfigs['Strength'],nperiods=UndulConfigs['NumPeriods'],lperiod=UndulConfigs['Period'])
	return LattObjs

def aligh_slices(p_arrays,QuadGradients=QuadGradients,DipAngles=DipAngles,UndulConfigs=UndulConfigs, BeamEnergy=BeamEnergy_ref,\
  BeamEnergy_ref=BeamEnergy_ref, stop_key=None, method=method):
	r56 = []
	sys.stdout.write('\nAligning '+str(len(p_arrays))+' slices');sys.stdout.flush()
	for i in range(len(p_arrays)):
		latt_elmts = make_line(QuadGradients=QuadGradients,DipAngles=DipAngles,UndulConfigs=UndulConfigs,\
		  BeamEnergy=p_arrays[i].E, BeamEnergy_ref=BeamEnergy_ref)
		if stop_key==None:
			lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys],method=method)
		else:
			lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys],method=method,stop=latt_elmts[stop_key])
		oclt.cpbd.optics.lattice_transfer_map_R(lat,p_arrays[i].E)
		r56.append(lat.R[4,5])

	ds = [0,]
	s_loc = 0.0
	for i in range(1,len(p_arrays)):
		de = 1-p_arrays[i].E/p_arrays[i-1].E
		s_loc += de*r56[i]
		ds.append(s_loc)
	for i in range(len(p_arrays)):
		p_arrays[i].s += ds[i]
	return p_arrays

def make_shot(p_arrays,QuadGradients=QuadGradients,DipAngles=DipAngles,UndulConfigs=UndulConfigs, BeamEnergy=BeamEnergy_ref,\
  BeamEnergy_ref=BeamEnergy_ref, stop_key=None, method=method, Nit = 1,output = None, damping = None):

	latt_elmts = make_line(UndulConfigs=UndulConfigs)
	if stop_key==None:
		lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys],method=method)
	else:
		lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys],method=method,stop=latt_elmts[stop_key])
	dz = lat.totalLen/Nit
	if type(p_arrays)!=list: p_arrays = [p_arrays,]

	outputs = {}
	if damping!=None:
		sys.stdout.write('Particle losses activated \n');sys.stdout.flush()
		outputs['staying'] = np.zeros(Nit)
		to_pop = []
	if output!=None:
		for key in output: outputs[key] = np.zeros(Nit)
		outputs['s'] = dz*np.arange(Nit)

	for i in range(len(p_arrays)):
		latt_elmts = make_line(QuadGradients=QuadGradients,DipAngles=DipAngles,UndulConfigs=UndulConfigs,\
        BeamEnergy=p_arrays[i].E, BeamEnergy_ref=BeamEnergy_ref)
		lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys],method=method,stop=latt_elmts['QEM4-UNDL'])
		navi = oclt.Navigator(lat)
		sss = '\r'+'Transporing slice '+str(i+1)+' of '+str(len(p_arrays))+':'
		#sys.stdout.write('\r'+'Transporing slice '+str(i+1)+' of '+str(len(p_arrays))+':')
		#sys.stdout.flush()
		for j in range(Nit):
			sys.stdout.write(sss+'step '+str(j+1)+' of '+str(Nit))
			oclt.tracking_step(lat, p_arrays[i], dz,navi)
			if damping!=None:
				p_arrays[i], np_loc = damp_particles(p_arrays[i],*damping)
				outputs['staying'][j] += np_loc
				if np_loc<1:
					to_pop.append(i)
					break
			if output!=None:
				for key in output: outputs[key][j] += beam_diags(p_arrays[i],key)

	if damping!=None:
		poped = 0
		for i in to_pop:
			p_arrays.pop(i-poped)
			poped += 1
		if poped>0: 
			sys.stdout.write('\n'+str(poped)+' slices are lost');sys.stdout.flush()

	return p_arrays, outputs

def beam_diags(p_array, key):
	if key == 'x':  return( p_array.x().sum())
	if key == 'y':  return( p_array.y().sum())
	if key == 'x2': return((p_array.x()**2).sum())
	if key == 'y2': return((p_array.y()**2).sum())
	if key == 'n':  return( p_array.size())

def damp_particles(p_array, Rx,Ry):
	Np = p_array.size
	s0 = p_array.s
	Rx0, Ry0 = Rx(s0), Ry(s0)
	indx = np.nonzero((np.abs(p_array.x())<Rx0)*(np.abs(p_array.y())<Ry0))[0]
	p_array.particles = p_array.particles.reshape((p_array.size(),6))[indx,:].flatten()
	return p_array, p_array.size()

