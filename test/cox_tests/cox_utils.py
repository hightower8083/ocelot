"""
This module contains the functions needed for OCELOT to model 
the COXINEL experiment. In particular it deals with the broad-spectrum 
beam transport.

by I.A. Andriyash
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import ocelot as oclt
from cox_configs import *
from scipy.constants import m_e,c,e

mc2_GeV = m_e*c**2/e*1e-9

method = oclt.MethodTM()
method.global_method = oclt.SecondTM

def sliced_spectrum(Emin,Emax,dg = 0.001,):
	"""
	Gets the array silced with fixed relative width

	Parameters
	----------
	Emin : float
	  Minimal value
	Emax : float
	  Maximal value
	dg : float
	  relative width of the slices

	Returns
	-------
	nrg_sliced : array
	  boundaries of the slices

	"""

	e_m = Emin
	e_p = Emin
	nrg_sliced = [e_m,]
	while e_p<=Emax:
		e_p = e_m*(1+dg)
		nrg_sliced.append(e_p)
		e_m = e_p
	return np.array(nrg_sliced)

def make_beam_sliced(bm,dg=0.002,div_chirp=None):
	"""
	Makes a contineous spectrum beam as a list of slices

	Parameters
	----------
	bm : dictionary
	  beam input parameters as defined for make_beam_contin function
	dg : float
	  relative width of the slices

	Returns
	-------
	p_arrays : list of oclt.ParticleArray objects
	  list of the beams

	See Also
	--------
	  make_beam and make_beam_contin

	"""

	p_arrays = []
	E1 = bm['E'][0]
	E2 = bm['E'][1]
	nrg_sliced =  sliced_spectrum(E1,E2,dg=dg)
	nrg_cntrs = 0.5*(nrg_sliced[1:]+nrg_sliced[:-1])
	Nbeams = nrg_sliced.shape[0]-1
	part_per_slice = np.round(bm['Np'] \
	  *(nrg_sliced[1:]-nrg_sliced[:-1]) \
	  /(nrg_sliced[-1]-nrg_sliced[0])).astype('i')
	for i in range(Nbeams):
		beam = oclt.deepcopy(bm)
		beam['Np'] = part_per_slice[i]
		beam['Q'] /= Nbeams
		beam['E'] = (nrg_sliced[i],nrg_sliced[i+1])
		p_arrays.append(make_beam_contin(beam,div_chirp=div_chirp))
	return p_arrays

def make_beam(bm):
	"""
	Makes a beam with Gaussian phase distributions

	Parameters
	----------
	bm : dictionary
	  beam input parameters (in SI units): 
	    Np -- full number of particles
	    Lz, Rx, Ry -- spatial dimensions of the beam
	    dE, Ox, Oy -- normalized momenta dispersions
	    E -- beam mean energy (in GeV)
	    Q -- full beam charge (in C)

	Returns
	-------
	p_array : an oclt.ParticleArray object

	"""

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

def make_beam_contin(bm, div_chirp=None):
	"""
	Makes a beam with Gaussian phase distributions, 
	exept for energy spread which is flat-top

	Parameters
	----------
	bm : dictionary
	  beam input parameters as in make_beam exept for: 
	    bm['E'] = (E1,E2) -- range of electron energies (in GeV)

	Returns
	-------
	p_array : an oclt.ParticleArray object

	See Also
	--------
	  make_beam

	"""
	p_array = oclt.ParticleArray()

	E1 = bm['E'][0]
	E2 = bm['E'][1]
	dE = (E2-E1)/(E2+E1)
	p_array.E = 0.5*(E2+E1)

	parts0 = np.zeros((6,bm['Np']))
	parts0[0] = bm['Rx']*np.random.randn(bm['Np'])
	parts0[2] = bm['Ry']*np.random.randn(bm['Np'])
	parts0[4] = bm['Lz']*np.random.randn(bm['Np'])
	parts0[5] = dE*(2*np.random.rand(bm['Np'])-1)

	if div_chirp!=None:
		pz0 = np.sqrt( (div_chirp/mc2_GeV)**2-1. )
		pz = (p_array.E/mc2_GeV)*(1+parts0[5])

		px =  bm['Ox']*pz0*np.random.randn(bm['Np'])
		py =  bm['Oy']*pz0*np.random.randn(bm['Np'])
		
		parts0[1] = px/pz
		parts0[3] = py/pz
	else:
		parts0[1] = bm['Ox']*np.random.randn(bm['Np'])
		parts0[3] = bm['Oy']*np.random.randn(bm['Np'])

	p_array.particles = parts0.T.flatten()
	p_array.s = 0.0
	if 'Q' in bm:
		p_array.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array.q_array = np.ones(bm['Np'])
	return p_array

def make_line(Drifts = Drifts,QuadLengths=QuadLengths,\
  QuadGradients=QuadGradients,DipLengths=DipLengths,\
  DipAngles=DipAngles,UndulConfigs=UndulConfigs,\
	BeamEnergy=BeamEnergy_ref, BeamEnergy_ref=BeamEnergy_ref):

	"""
	Makes a COXINEL-type transport line with drifts, dipoles, 
	quadrupoles, undulator and few screens

	Parameters
	----------
	Drifts, QuadLengths, QuadGradients, DipLengths,\
	 DipAngles,UndulConfigs : dictionaries with elements as 
	 defined in cox_configs
	BeamEnergy: float
	  energy of the beam (in GeV)
	BeamEnergy_ref: float
	  reference energy of the beam for which lattice gradients 
	  and angles are defined (in GeV)

	Returns
	-------
	LattObjs : dictionary with ocelot.cpbd.elements

	See Also
	--------
	  cox_configs

	"""

	DriftObjs = {}
	QuadObjs = {}
	DipObjs = {}
	LattObjs = {None:None,}
	Latt = []
	Ecorr = BeamEnergy_ref/BeamEnergy

	for key in Drifts.keys(): 
		LattObjs[key] = oclt.Drift(l=Drifts[key])

	for key in QuadLengths.keys():
		LattObjs[key] = oclt.Quadrupole( \
		  l=QuadLengths[key], k1=QuadGradients[key]*Ecorr )
		
	for key in DipLengths.keys():
		if key == 'DIP1':
			e1 = 0.0;
			e2 = DipAngles[key]*Ecorr;
		elif key == 'DIP2':
			e1 = DipAngles[key]*Ecorr;
			e2 =  0.0;
		elif key == 'DIP3':
			e1 =  0.0;
			e2 = -DipAngles[key]*Ecorr;
		elif key == 'DIP4':
			e1 = -DipAngles[key]*Ecorr;
			e2 = 0.0;

		LattObjs[key] = oclt.RBend( \
		  l=DipLengths[key], angle=DipAngles[key]*Ecorr)

	for key in  ['IMG1','IMG2','IMG4','IMG5',]: 
		LattObjs[key] = oclt.Marker()

	LattObjs['UNDL1'] = oclt.Undulator(Kx=UndulConfigs['Strength'],\
	  nperiods=UndulConfigs['NumPeriods'],lperiod=UndulConfigs['Period'])
	LattObjs['UNDL2'] = oclt.Undulator(Kx=UndulConfigs['Strength'],\
	  nperiods=UndulConfigs['NumPeriods'],lperiod=UndulConfigs['Period'])
	return LattObjs

def make_shot(p_arrays,QuadGradients=QuadGradients,DipAngles=DipAngles,\
  UndulConfigs=UndulConfigs, BeamEnergy=BeamEnergy_ref, \
  BeamEnergy_ref=BeamEnergy_ref, stop_key=None, method=method, \
  Nit = 1,output = None, damping = None):

	"""
	Performs the transport simulation through a COXINEL-type transport line.

	Parameters
	----------
	p_arrays: list or single oclt.ParticleArray object
	Drifts, QuadLengths, QuadGradients, DipLengths,\
	 DipAngles,UndulConfigs,BeamEnergy, BeamEnergy_ref : line parameters as 
	 defined in make_line
	stop_key : str
	  key-name of lattice element at which to stop the transport
	method : oclt.methodTM
	  transport method (works with first or second orders)
	Nit : int
	  number of steps to perform (for output and damping resolutions)
	output : dictionary or None
	  outputs to perform as defined in beam_diags
	damping: list, tuple or None
	  If not None the particle losses in the limitied width pipe is 
	  performed as defined in damp_particles

	Returns
	-------
	p_arrays : list of oclt.ParticleArray objects
	outputs : dictionary with arrays

	See Also
	--------
	  make_line, beam_diags, damp_particles

	"""

	latt_elmts = make_line(UndulConfigs=UndulConfigs)
	if stop_key==None:
		lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys], \
		  method=method)
	else:
		lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys],\
		  method=method,stop=latt_elmts[stop_key])
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
		latt_elmts = make_line(QuadGradients=QuadGradients,DipAngles=DipAngles, \
		  UndulConfigs=UndulConfigs, BeamEnergy=p_arrays[i].E, \
		  BeamEnergy_ref=BeamEnergy_ref)
		lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys], \
		  method=method,stop=latt_elmts[stop_key])
		navi = oclt.Navigator(lat)
		sss = '\r'+'Transporing slice '+str(i+1)+' of '+str(len(p_arrays))+':'
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
			sys.stdout.write('\n'+str(poped)+' slices are lost')
			sys.stdout.flush()

	return p_arrays, outputs

def aligh_slices(p_arrays,QuadGradients=QuadGradients, \
  DipAngles=DipAngles,UndulConfigs=UndulConfigs, BeamEnergy=BeamEnergy_ref,\
  BeamEnergy_ref=BeamEnergy_ref, stop_key=None, method=method):

	"""
	Alignes the slices of the spectrally sliced beam at the end of 
	the transport simulation.

	Parameters
	----------
	p_arrays: list of oclt.ParticleArray objects
	Drifts, QuadLengths, QuadGradients, DipLengths,\
	 DipAngles,UndulConfigs,BeamEnergy, BeamEnergy_ref : line parameters as 
	 defined in make_line and make_shot
	stop_key : str
	  key-name of lattice element at which to stop the transport
	method : oclt.methodTM
	  transport method (works with first or second orders)

	Returns
	-------
	p_arrays : list or single oclt.ParticleArray object

	See Also
	--------
	  make_line, make_shot

	"""

	r56 = []
	sys.stdout.write('\nAligning '+str(len(p_arrays))+' slices')
	sys.stdout.flush()
	for i in range(len(p_arrays)):
		latt_elmts = make_line(QuadGradients=QuadGradients, \
		  DipAngles=DipAngles,UndulConfigs=UndulConfigs,\
		  BeamEnergy=p_arrays[i].E, BeamEnergy_ref=BeamEnergy_ref)
		if stop_key==None:
			lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys], \
			  method=method)
		else:
			lat = oclt.MagneticLattice([latt_elmts[key] for key in cell_keys], \
			  method=method,stop=latt_elmts[stop_key])
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

def beam_diags(p_array, key):
	"""
	Returns the sums of particles coordinates/momenta defined by key-names

	Parameters
	----------
	p_array: oclt.ParticleArray object
	key : str
	  key-name of the coordinate to sum and return; can be:
	    'x', 'y' for sum(x) and sum(y)
	    'x2', 'y2' for sum(x**2) and sum(y**2) *
	    'px', 'py' for sum(px) and sum(py)
	    'px2', 'py2' for sum(px**2) and sum(py**2) *
	    'n' for full number of particles
	Returns
	-------
	value : float

	See Also
	--------
	  make_shot

	"""

	if key == 'x':  val = p_array.x().sum()
	if key == 'y':  val = p_array.y().sum()
	if key == 'z':  val = -p_array.tau().sum()
	if key == 'x2': val = (p_array.x()**2).sum()
	if key == 'y2': val = (p_array.y()**2).sum()
	if key == 'z2': val = (p_array.tau()**2).sum()
	if key == 'px':  val = p_array.px().sum()
	if key == 'py':  val = p_array.py().sum()
	if key == 'px2': val = (p_array.px()**2).sum()
	if key == 'py2': val = (p_array.py()**2).sum()
	if key == 'n':  val =  p_array.size()
	return val

def damp_particles(p_array, Rx,Ry):
	"""
	Dumps the particles which overpass the limited pipe width

	Parameters
	----------
	p_array: oclt.ParticleArray object
	Rx,Ry : python functions
	  X and Y pipe widths as the functions of distance: Rx(s) and Ry(s)
	Returns
	-------
	p_array: oclt.ParticleArray object
	Np : Number of remaining particles

	See Also
	--------
	  make_shot

	"""

	s0 = p_array.s
	Rx0, Ry0 = Rx(s0), Ry(s0)
	indx = np.nonzero((np.abs(p_array.x())<Rx0)*(np.abs(p_array.y())<Ry0))[0]
	p_array.particles = \
	  p_array.particles.reshape((p_array.size(),6))[indx,:].flatten()
	Np = p_array.size()
	return p_array, Np

def ocelot_to_chimera(p_arrays,beam,lam0,keep_orig=True):
	"""
	Exports the list of ParticleArrays from OCELOT to CHIMERA

	Parameters
	----------
	p_arrays: list of oclt.ParticleArray objects
	  beam represented in the form of slices list
	beam : chimera.Species object
	  CHIMERA species obejct to populate with particles
	lam0: float
	  normalization length for CHIMERA in meters

	Returns
	-------
	beam : chimera.Species object
		CHIMERA species obejct populated with particles
	"""

	Np = np.sum([p_array.size() for p_array in p_arrays])
	g0 = p_array.E/mc2_GeV

	xx = np.hstack(([-p_array.tau()/lam0 for p_array in p_arrays]))
	yy = np.hstack(([p_array.y()/lam0 for p_array in p_arrays]))
	zz = np.hstack(([p_array.x()/lam0 for p_array in p_arrays]))

	gg = np.hstack(([(p_array.p()+1)*g0 for p_array in p_arrays]))
	oy = np.hstack(([p_array.py() for p_array in p_arrays]))
	oz = np.hstack(([p_array.px() for p_array in p_arrays]))

	qq = np.hstack(([p_array.q_array*1e12 for p_array in p_arrays]))

	px = np.sqrt( (gg**2-1.)/(1+oy**2+oz**2) )
	py = px*oy
	pz = px*oz
	gg, oy, oz = None, None, None

	beam.Data['coords'] = np.zeros((3,Np))
	beam.Data['momenta'] = np.zeros((3,Np))
	beam.Data['weights'] = np.zeros(Np)

	beam.Data['coords'][0] = xx
	beam.Data['coords'][1] = yy
	beam.Data['coords'][2] = zz

	beam.Data['momenta'][0] = px
	beam.Data['momenta'][1] = py
	beam.Data['momenta'][2] = pz

	beam.Data['weights'][:] = -qq/(beam.weight2pC*lam0*1e6)
	beam.Data['coords'][0] -= (beam.Data['coords'][0]*beam.Data['weights']\
	  ).sum()/(beam.Data['weights']).sum()
	beam.Data['coords_halfstep'] = beam.Data['coords'].copy()
	if keep_orig is False: 
		del p_arrays
	return  beam

def print_latt_chars(lat, BeamEnergy):
	"""
	Prints some elements of the transport matrices
	Parameters
	----------
	lat: MagneticLattice object
	BeamEnergy : float
	  Electron energy in GEVs
   """
	oclt.cpbd.optics.lattice_transfer_map_RT(lat,BeamEnergy)
	params = [lat.R[0,0],lat.R[2,2],lat.R[4,5]*1000,\
	  lat.T[1,1,5], lat.T[3,3,5],lat.T[0,1,5],lat.T[2,3,5]]
	print latt_par_string.format(*params)
