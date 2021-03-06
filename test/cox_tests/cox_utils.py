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
from scipy.constants import m_e,c,e
from scipy.optimize import fmin

np.set_printoptions(precision=6, suppress=True, linewidth=120)

mc2_GeV = m_e*c**2/e*1e-9
method = oclt.MethodTM()
method.global_method = oclt.SecondTM

latt_par_string = """
###############################
R_11 = {0:.3g}, R_33 = {1:.3g},
R_56 = {2:.3g} mm,
R_226 = {3:.3g}, R_446 = {4:.3g},
R_126 = {5:.3g}, R_346 = {6:.3g}
################################
"""

def sliced_spectrum(Emin,Emax,dg = 0.001,):
	"""
	Gets the array sliced with fixed relative width

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

	p_array = oclt.ParticleArray()

	if 'Features' not in bm:
		bm['Features'] = {}
	if 'RandomSeed' not in bm['Features']:
		bm['Features']['RandomSeed'] = int(np.random.rand()*1e7)
	rnd = np.random.RandomState(seed = bm['Features']['RandomSeed'])

	if 'X0' not in bm: bm['X0'] = 0.
	if 'Y0' not in bm: bm['Y0'] = 0.
	if 'Z0' not in bm: bm['Z0'] = 0.

	parts0 = np.zeros((6,bm['Np']))

	if 'FlatX' in bm['Features']:
		parts0[0] = bm['X0'] + bm['Rx']*(2*rnd.rand(bm['Np'])-1)
	else:
		parts0[0] = bm['X0'] + bm['Rx']*rnd.randn(bm['Np'])

	if 'FlatY' in bm['Features']:
		parts0[2] = bm['Y0'] + bm['Ry']*(2*rnd.rand(bm['Np'])-1)
	else:
		parts0[2] = bm['Y0'] + bm['Ry']*rnd.randn(bm['Np'])

	if 'FlatZ' in bm['Features']:
		parts0[4] = bm['Z0'] + bm['Lz']*(2*rnd.rand(bm['Np'])-1)
	else:
		parts0[4] = bm['Z0'] + bm['Lz']*rnd.randn(bm['Np'])

	parts0[1] = bm['Ox']*rnd.randn(bm['Np'])
	parts0[3] = bm['Oy']*rnd.randn(bm['Np'])
	parts0[5] = bm['dE']*rnd.randn(bm['Np'])

	p_array.particles = parts0.T.flatten()
	p_array.E = bm['E']
	p_array.s = 0.0

	p_array.x_c = 0.0
	p_array.y_c = 0.0
	p_array.z_c = 0.0

	if 'Q' in bm:
		p_array.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array.q_array = np.ones(bm['Np'])
	return p_array

def make_beam_sliced(bm,dg=0.002,div_chirp=None):
	"""
	Makes a continuous spectrum beam as a list of slices

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
	if 'Features' not in bm:
		bm['Features'] = {}
	if 'RandomSeed' not in bm['Features']:
		bm['Features']['RandomSeed'] = int(np.random.rand()*1e7)
	rnd = np.random.RandomState(seed=bm['Features']['RandomSeed'])

	E1 = bm['E'][0]
	E2 = bm['E'][1]
	nrg_sliced =  sliced_spectrum(E1,E2,dg=dg)
	nrg_cntrs = 0.5*(nrg_sliced[1:]+nrg_sliced[:-1])
	Nbeams = nrg_sliced.shape[0]-1
	part_per_slice = np.round(bm['Np'] \
	  *(nrg_sliced[1:]-nrg_sliced[:-1]) \
	  /(nrg_sliced[-1]-nrg_sliced[0])).astype('i')
	rnd_seeds = (rnd.rand(Nbeams)*1e7).astype('l')

	for i in range(Nbeams):
		beam = oclt.deepcopy(bm)
		beam['Np'] = part_per_slice[i]
		beam['Q'] /= int(bm['Np']/beam['Np'])
		beam['E'] = (nrg_sliced[i],nrg_sliced[i+1])
		beam['Features']['RandomSeed'] = rnd_seeds[i]
		p_arrays.append(make_beam_contin(beam,div_chirp=div_chirp))
	return p_arrays

def make_beam_contin(bm, div_chirp=None):
	"""
	Makes a beam with Gaussian phase distributions,
	except for energy spread which is flat-top

	Parameters
	----------
	bm : dictionary
	  beam input parameters as in make_beam except for:
	    bm['E'] = (E1,E2) -- range of electron energies (in GeV)

	Returns
	-------
	p_array : an oclt.ParticleArray object

	See Also
	--------
	  make_beam

	"""
	p_array = oclt.ParticleArray()

	if 'Features' not in bm:
		bm['Features'] = {}
	if 'RandomSeed' not in bm['Features']:
		bm['Features']['RandomSeed'] = int(rnd.rand()*1e7)
	rnd = np.random.RandomState(seed=bm['Features']['RandomSeed'])

	E1 = bm['E'][0]
	E2 = bm['E'][1]
	dE = (E2-E1)/(E2+E1)
	p_array.E = 0.5*(E2+E1)

	if 'X0' not in bm: bm['X0'] = 0.
	if 'Y0' not in bm: bm['Y0'] = 0.
	if 'Z0' not in bm: bm['Z0'] = 0.

	parts0 = np.zeros((6,bm['Np']))

	if 'FlatX' in bm['Features']:
		parts0[0] = bm['X0'] + bm['Rx']*(2*rnd.rand(bm['Np'])-1)
	else:
		parts0[0] = bm['X0'] + bm['Rx']*rnd.randn(bm['Np'])

	if 'FlatY' in bm['Features']:
		parts0[2] = bm['Y0'] + bm['Ry']*(2*rnd.rand(bm['Np'])-1)
	else:
		parts0[2] = bm['Y0'] + bm['Ry']*rnd.randn(bm['Np'])

	if 'FlatZ' in bm['Features']:
		parts0[4] = bm['Z0'] + bm['Lz']*(2*rnd.rand(bm['Np'])-1)
	else:
		parts0[4] = bm['Z0'] + bm['Lz']*rnd.randn(bm['Np'])

	parts0[5] = dE*(2*rnd.rand(bm['Np'])-1)

	if div_chirp!=None:
		pz0 = np.sqrt( (div_chirp/mc2_GeV)**2-1. )
		pz = (p_array.E/mc2_GeV)*(1+parts0[5])

		px =  bm['Ox']*pz0*rnd.randn(bm['Np'])
		py =  bm['Oy']*pz0*rnd.randn(bm['Np'])

		parts0[1] = bm['Ox']*rnd.randn(bm['Np'])*(pz0/pz)**1.4
		parts0[3] = bm['Oy']*rnd.randn(bm['Np'])*(pz0/pz)**1.4

	else:
		parts0[1] = bm['Ox']*rnd.randn(bm['Np'])
		parts0[3] = bm['Oy']*rnd.randn(bm['Np'])

	p_array.particles = parts0.T.flatten()
	p_array.s = 0.0
	if 'Q' in bm:
		p_array.q_array = (bm['Q']/bm['Np'])*np.ones(bm['Np'])
	else:
		p_array.q_array = np.ones(bm['Np'])

	p_array.x_c = 0.0
	p_array.y_c = 0.0
	p_array.z_c = 0.0

	return p_array

def make_line(lattice_elements, cell_keys, BeamEnergy, BeamEnergy_ref):
	"""
	Makes a COXINEL-type transport line with drifts, dipoles,
	quadrupoles, undulator and few screens

	Parameters
	----------
	lattice_elements: dictionary which contains COXINEL elements 
	  as defined in cox_configs (shoud include Drifts, QuadLengths, 
	  QuadGradients, DipLengths, DipAngles,UndulConfigs)
	cell_keys: list with the keys to define the lattice sequence
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

	for key in lattice_elements['Drifts'].keys():
		LattObjs[key] = oclt.Drift(l=lattice_elements['Drifts'][key])

	for key in lattice_elements['QuadLengths'].keys():
		LattObjs[key] = oclt.Quadrupole( \
		  l=lattice_elements['QuadLengths'][key], \
		  k1 = lattice_elements['QuadGradients'][key]*Ecorr )

	for key in lattice_elements['DipLengths'].keys():
		if key == 'DIP1':
			e1 = 0.0;
			e2 = lattice_elements['DipAngles'][key]*Ecorr;
		elif key == 'DIP2':
			e1 = lattice_elements['DipAngles'][key]*Ecorr;
			e2 =  0.0;
		elif key == 'DIP3':
			e1 =  0.0;
			e2 = -lattice_elements['DipAngles'][key]*Ecorr;
		elif key == 'DIP4':
			e1 = -lattice_elements['DipAngles'][key]*Ecorr;
			e2 = 0.0;

		LattObjs[key] = oclt.RBend( \
		  l = lattice_elements['DipLengths'][key], 
		  angle = lattice_elements['DipAngles'][key]*Ecorr )

	for key in  ['IMG1','IMG2','IMG4','IMG5',]:
		LattObjs[key] = oclt.Marker()

	LattObjs['UNDL1'] = oclt.Undulator(\
	  Kx = lattice_elements['UndulConfigs']['Strength'],\
	  nperiods = lattice_elements['UndulConfigs']['NumPeriods'],\
	  lperiod = lattice_elements['UndulConfigs']['Period'])
	LattObjs['UNDL2'] = oclt.Undulator(\
	  Kx = lattice_elements['UndulConfigs']['Strength'],\
	  nperiods = lattice_elements['UndulConfigs']['NumPeriods'],\
	  lperiod = lattice_elements['UndulConfigs']['Period'])

	LattObjs = [LattObjs[key] for key in cell_keys]

	return LattObjs

def make_shot(p_arrays, lattice_elements, cell_keys, \
  BeamEnergy_ref,start_key=None, stop_key=None, \
  method=method, Nit = 1,output = None,verbose=False):

	"""
	Performs the transport simulation through a COXINEL-type transport line.

	Parameters
	----------
	p_arrays: list or single oclt.ParticleArray object
	lattice_elements: line parameters as defined in make_line
	stop_key : str
	  key-name of lattice element at which to stop the transport (exclusive)
	method : oclt.methodTM
	  transport method (works with first or second orders)
	Nit : int
	  number of steps to perform (for output and damping resolutions)
	output : dictionary or None
	  outputs to perform as defined in beam_diags

	Returns
	-------
	p_arrays : list of oclt.ParticleArray objects
	outputs : dictionary with arrays

	See Also
	--------
	  make_line, beam_diags, damp_particles

	"""

	latt_objs = make_line(lattice_elements, cell_keys, \
	  BeamEnergy=BeamEnergy_ref, BeamEnergy_ref=BeamEnergy_ref,)

	if stop_key!=None and start_key!=None:
		stop_indx = cell_keys.index(stop_key)-1
		start_indx = cell_keys.index(start_key)
		lat = oclt.MagneticLattice(latt_objs, method=method, \
		  start=latt_objs[start_indx], stop=latt_objs[stop_indx])
	elif stop_key!=None and start_key==None:
		stop_indx = cell_keys.index(stop_key)-1
		lat = oclt.MagneticLattice(latt_objs, method=method, \
		  stop=latt_objs[stop_indx])
	elif stop_key==None and start_key!=None:
		start_indx = cell_keys.index(start_key)
		lat = oclt.MagneticLattice(latt_objs, method=method,\
		  start=latt_objs[start_indx])
	else:
		lat = oclt.MagneticLattice(latt_objs, method=method)

	dz = lat.totalLen/Nit
	if type(p_arrays)!=list: p_arrays = [p_arrays,]

	outputs = {}
	if output!=None:
		for key in output: outputs[key] = np.zeros(Nit+1)

	for i in range(len(p_arrays)):
		latt_objs = make_line(lattice_elements, cell_keys, \
		  BeamEnergy=p_arrays[i].E, BeamEnergy_ref=BeamEnergy_ref)

		if stop_key!=None and start_key!=None:
			stop_indx = cell_keys.index(stop_key)-1
			start_indx = cell_keys.index(start_key)
			lat = oclt.MagneticLattice(latt_objs, method=method,\
			  start=latt_objs[start_indx],stop=latt_objs[stop_indx])
		elif stop_key!=None and start_key==None:
			stop_indx = cell_keys.index(stop_key)-1
			lat = oclt.MagneticLattice(latt_objs, method=method,\
			  stop=latt_objs[stop_indx])
		elif stop_key==None and start_key!=None:
			start_indx = cell_keys.index(start_key)
			lat = oclt.MagneticLattice(latt_objs,method=method,\
			  start=latt_objs[start_indx])
		else:
			lat = oclt.MagneticLattice(latt_objs, method=method)

		navi = oclt.Navigator(lat)
		sss = '\r'+'Transporing slice '+str(i+1)+' of '+str(len(p_arrays))+':'
		for j in range(Nit+1):
			if verbose: sys.stdout.write(sss+'step '+str(j+1)+' of '+str(Nit))
			oclt.tracking_step(lat, p_arrays[i], dz,navi)
			if output!=None:
				for key in output:
					outputs[key][j] += beam_diags(p_arrays[i],key)

	return p_arrays, outputs

def aligh_slices(p_arrays, lattice_elements, cell_keys,\
  BeamEnergy_ref, start_key=None, stop_key=None,\
  method=method,verbose=False):

	"""
	Aligns the slices of the spectrally sliced beam at the end of 
	the transport simulation.

	Parameters
	----------
	p_arrays: list of oclt.ParticleArray objects
	lattice_elements: line parameters as defined in make_line and make_shot
	stop_key : str
	  key-name of lattice element at which to stop the transport (exclusive)
	method : oclt.methodTM
	  transport method (works with first or second orders)

	Returns
	-------
	p_arrays : list or single oclt.ParticleArray object

	See Also
	--------
	  make_line, make_shot

	"""

	if verbose:
		sys.stdout.write('\nAligning '+str(len(p_arrays))+' slices')
		sys.stdout.flush()

	r16 = []
	r36 = []
	r56 = []
	method_ord1 = oclt.MethodTM()

	for i in range(len(p_arrays)):
		latt_objs = make_line(lattice_elements, cell_keys, \
		  BeamEnergy=p_arrays[i].E, BeamEnergy_ref=BeamEnergy_ref)

		if stop_key!=None and start_key!=None:
			stop_indx = cell_keys.index(stop_key)-1
			start_indx = cell_keys.index(start_key)
			lat = oclt.MagneticLattice(latt_objs, method=method_ord1,\
			  start=latt_objs[start_indx],stop=latt_objs[stop_indx])
		elif stop_key!=None and start_key==None:
			stop_indx = cell_keys.index(stop_key)-1
			lat = oclt.MagneticLattice(latt_objs, method=method_ord1,\
			  stop=latt_objs[stop_indx])
		elif stop_key==None and start_key!=None:
			start_indx = cell_keys.index(start_key)
			lat = oclt.MagneticLattice(latt_objs, method=method_ord1,\
			  start=latt_objs[start_indx])
		else:
			lat = oclt.MagneticLattice(latt_objs, method=method_ord1)

		oclt.cpbd.optics.lattice_transfer_map_R(lat,p_arrays[i].E)
		r16.append(lat.R[0,5])
		r36.append(lat.R[2,5])
		r56.append(lat.R[4,5])

	dx = [0,]
	dy = [0,]
	ds = [0,]
	x_loc = 0.0
	y_loc = 0.0
	s_loc = 0.0

	for i in range(1,len(p_arrays)):
		de = 1-p_arrays[i].E/p_arrays[i-1].E
		x_loc += de*r16[i]
		y_loc += de*r36[i]
		s_loc += de*r56[i]
		dx.append(x_loc)
		dy.append(y_loc)
		ds.append(s_loc)

	ee = np.array([p_array.E for p_array in p_arrays])
#	ee_centr = 0.5*(ee.max()+ee.min())
	ee_centr = BeamEnergy_ref
	indx_centr = (ee_centr-ee>0).sum()

	dx = np.array(dx)
	dy = np.array(dy)
	ds = np.array(ds)
	dx -= dx[indx_centr]
	dy -= dy[indx_centr]
	ds -= ds[indx_centr]
	for i in range(len(p_arrays)):
		p_arrays[i].x_c = -dx[i]
		p_arrays[i].y_c = -dy[i]
		p_arrays[i].z_c = ds[i]

	return p_arrays

def insert_slit(p_arrays,cntr = 0.0, width=np.inf, comp='x'):
	to_pop = []
	for i in range(len(p_arrays)):
		if comp=='x':
			x = p_arrays[i].x() + p_arrays[i].x_c
			coord = (x-cntr)/width*2
		elif comp=='y':
			y = p_arrays[i].y() + p_arrays[i].y_c
			coord = (y-cntr)/width*2
		elif comp=='rectangular':
			x = p_arrays[i].x() + p_arrays[i].x_c
			y = p_arrays[i].y() + p_arrays[i].y_c
			coord = (x-cntr[0])**x + (y-cntr[1])**2
		elif comp=='ellipse':
			x = p_arrays[i].x() + p_arrays[i].x_c
			y = p_arrays[i].y() + p_arrays[i].y_c
			coord = np.sqrt((x-cntr[0])**2/width[0]**2 \
			  + (y-cntr[1])**2/width[1]**2)

		Num_loc = coord.shape[0]
		indx = np.nonzero(np.abs(coord)<=1)[0]

		if indx.shape[0]==0:
			to_pop.append(i)
		else:
			p_arrays[i].particles = p_arrays[i].particles \
			  .reshape((Num_loc,6))[indx,:].flatten()
			p_arrays[i].q_array = p_arrays[i].q_array[indx]

	for i in to_pop[::-1]: 
		p_arrays.pop(i)
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
	if key == 's':  val = p_array.s*p_array.size()

	if key == 'x2': val = (p_array.x()**2).sum()
	if key == 'y2': val = (p_array.y()**2).sum()
	if key == 'z2': val = (p_array.tau()**2).sum()

	if key == 'xp':  val = p_array.px().sum()
	if key == 'yp':  val = p_array.py().sum()
	if key == 'zp':  val = p_array.p().sum()

	if key == 'xp2': val = (p_array.px()**2).sum()
	if key == 'yp2': val = (p_array.py()**2).sum()
	if key == 'zp2': val = (p_array.p()**2).sum()

	if key == 'px':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		val = (p_array.px()*gg).sum()
	if key == 'py':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		val = (p_array.py()*gg).sum()
	if key == 'pz':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		val = gg.sum()

	if key == 'px2':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		pz = np.sqrt(gg**2-1/(1+p_array.px()**2+p_array.py()**2))
		val = ((p_array.px()*pz)**2).sum()
	if key == 'py2':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		pz = np.sqrt(gg**2-1/(1+p_array.px()**2+p_array.py()**2))
		val = ((p_array.py()*pz)**2).sum()
	if key == 'pz2':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		pz = np.sqrt(gg**2-1/(1+p_array.px()**2+p_array.py()**2))
		val = (pz**2).sum()

	if key == 'xxp': val = (p_array.x()*p_array.px()).sum()
	if key == 'yyp': val = (p_array.y()*p_array.py()).sum()

	if key == 'xpx':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		pz = np.sqrt(gg**2-1/(1+p_array.px()**2+p_array.py()**2))
		val = (p_array.x()*p_array.px()*pz).sum()
	if key == 'ypy':
		gg = (1+p_array.p())*p_array.E/mc2_GeV
		pz = np.sqrt(gg**2-1/(1+p_array.px()**2+p_array.py()**2))
		val = (p_array.y()*p_array.py()*pz).sum()

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
	  p_array.particles.reshape((np.int(p_array.size()),6))[indx,:].flatten()
	p_array.q_array = p_array.q_array[indx]
	Np = np.int(p_array.size())
	return p_array, Np

def print_latt_chars(latt_objs, BeamEnergy):
	"""
	Prints some elements of the transport matrices
	Parameters
	----------
	lat: MagneticLattice object
	BeamEnergy : float
	  Electron energy in GEVs
   """
	lat = oclt.MagneticLattice(latt_objs,method=method)
	oclt.cpbd.optics.lattice_transfer_map_RT(lat,BeamEnergy)
	params = [lat.R[0,0],lat.R[2,2],lat.R[4,5]*1000,\
	  lat.T[1,1,5], lat.T[3,3,5],lat.T[0,1,5],lat.T[2,3,5]]
	print(latt_par_string.format(*params))

R_par_str =\
"""
 ( {:g}, {:g}, {:g}, {:g}, {:g}, {:g} )
 ( {:g}, {:g}, {:g}, {:g}, {:g}, {:g} )
 ( {:g}, {:g}, {:g}, {:g}, {:g}, {:g} )
 ( {:g}, {:g}, {:g}, {:g}, {:g}, {:g} )
 ( {:g}, {:g}, {:g}, {:g}, {:g}, {:g} )
 ( {:g}, {:g}, {:g}, {:g}, {:g}, {:g} )
"""

def print_R(lat, BeamEnergy):
	"""
	Prints some elements of the transport matrices
	Parameters
	----------
	lat: MagneticLattice object
	BeamEnergy : float
	  Electron energy in GEVs
   """
	oclt.cpbd.optics.lattice_transfer_map_R(lat,BeamEnergy)

	print(R_par_str.format(*lat.R.flatten()))

def supermatch_matrix(latt_objs,nrg,varz,MAG=True, printout=True):
	lat = oclt.MagneticLattice(latt_objs,method=method)

	def error_func(x):
		for i in range(len(varz)):
			if varz[i].__class__ == oclt.Quadrupole:
				varz[i].k1 = x[i]
				varz[i].transfer_map = lat.method.create_tm(varz[i])

		oclt.cpbd.optics.lattice_transfer_map_RT(lat, nrg)

        
		err = np.abs( lat.R[0,1])**2 + np.abs(lat.R[2,3] )**2 \
		    + np.abs(lat.T[1,1,5])**2 + np.abs(lat.T[3,3,5])**2 \
		    + (np.abs(lat.T[0,1,5])-np.abs(lat.T[2,3,5]))**2

		if MAG is True:
			err += (np.abs(lat.R[0,0])-np.abs(lat.R[2,2]))**2
		elif MAG is not True and MAG is not None:
			err += (np.abs(lat.R[0,0])-MAG[0])**2 + (np.abs(lat.R[2,2])-MAG[1])**2

		return err

	x = [0.0] * len(varz)

	for i in range(len(varz)):
		if varz[i].__class__ == oclt.Quadrupole:
			x[i] = varz[i].k1

	if printout: print("initial value: x = ", x)
	fmin(error_func, x, xtol=1e-9, maxiter=20000, maxfun=20000,disp=printout)
    
	for i in range(len(varz)):
		if varz[i].__class__ == oclt.Quadrupole:
			x[i] = varz[i].k1
	if printout: print("final value: x = ", x)
	return x
