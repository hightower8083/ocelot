import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt2d

import ocelot as oclt
from ocelot.test.cox_tests.cox_configs import *
from ocelot.test.cox_tests import cox_utils as cox

from ipywidgets import interactive

method = oclt.MethodTM()
method.global_method = oclt.SecondTM

def get_envs(QAP1=1., QAP2=1., QAP3=1.):
	energies = np.r_[0.02:0.25:300j]
	mags = []
    
	for BeamEnergy in energies:
		indx_stop = np.nonzero(np.array(cell_keys)=='IMG1')[0][0]
		QuadGradients1 = oclt.deepcopy(QuadGradients)
		QuadGradients1['QAP1'] *= QAP1
		QuadGradients1['QAP2'] *= QAP2
		QuadGradients1['QAP3'] *= QAP3
		latt_elmts = cox.make_line( \
		  BeamEnergy=BeamEnergy,QuadGradients=QuadGradients1 )
		lat = oclt.MagneticLattice( \
		  [latt_elmts[key] for key in cell_keys[:indx_stop+1]], \
		  method=method,stop=latt_elmts['IMG1'])
		oclt.cpbd.optics.lattice_transfer_map_R(lat,BeamEnergy)
		mags.append([lat.R[0,1],lat.R[2,3],lat.R[0,0],lat.R[2,2]])
    
	mags = np.array(mags)
	sx = np.abs(mags[:,0]*3e-3 + mags[:,2]*1e-6)*1e3
	sy = np.abs(mags[:,1]*3e-3 + mags[:,3]*1e-6)*1e3

	plt.plot(energies*1e3,sx,lw=1.5)
	plt.plot(energies*1e3,sy,lw=1.5)
	plt.plot(energies*1e3,np.sqrt(sx*sy),'--',c='k',lw=1.5)
	plt.legend(('$\sigma_x$','$\sigma_z$','$\sqrt{\sigma_x\sigma_z}$'),loc=1)
	plt.ylabel('Beam sizes (mm)')
	plt.xlabel('Electron energy (MeV)')
	plt.ylim(0,15);
    
def plot_beam(v,beam):
	p_arrays_init = cox.make_beam_sliced(beam, \
	  div_chirp=0.5*(beam['E'][0]+beam['E'][1]))

	fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
	p_arrays = oclt.deepcopy(p_arrays_init)

	QuadGradients1 = oclt.deepcopy(QuadGradients)
	QuadGradients1['QAP1'] *= v.children[0].value
	QuadGradients1['QAP2'] *= v.children[1].value
	QuadGradients1['QAP3'] *= v.children[2].value

    
	p_arrays,outputs = cox.make_shot(p_arrays, \
	  QuadGradients=QuadGradients1,stop_key='IMG1',)
	p_arrays = cox.aligh_slices(p_arrays, \
	  QuadGradients=QuadGradients1,stop_key='IMG1')
    
	xx = np.hstack(([p_array.x()*1e3 for p_array in p_arrays]))
	yy = np.hstack(([p_array.y()*1e3 for p_array in p_arrays]))
	ee = np.hstack(([(p_array.p()+1)*p_array.E*1e3 for p_array in p_arrays])) 
	rr = np.sqrt(xx**2+yy**2)

	h,x,y = np.histogram2d(xx,yy,bins=200,range=[[-0.5,0.5],[-0.5,0.5]])
	pl = ax1.imshow(medfilt2d(h,5).T,cmap=plt.cm.spectral,aspect='auto', \
	  extent=(x.min(),x.max(),y.min(),y.max()),origin='lower')

	h,x,y = np.histogram2d(ee,rr,bins=120, \
	  range=[[50,200],[0,3]])
	pl = ax2.imshow(medfilt2d(np.abs(h),5).T,cmap=plt.cm.spectral,\
	  aspect='auto',extent=(x.min(),x.max(),y.min(),y.max()),origin='lower')

	ax1.set_xlabel('X (mm)')
	ax1.set_ylabel('Y (mm)')
	ax2.set_xlabel('Electron energy (MeV)')
	ax2.set_ylabel('Radius (mm)')
	print('\nFraction of charge in 10 um cylinder is {0:g}%'.format( \
	  (rr<10e-3).sum()*100./rr.shape[0] )  )


qap_widgt = interactive(get_envs, \
  QAP1=(0.87,1.6,0.02), \
  QAP2=(0.88,1.6,0.02), \
  QAP3=(0.96,1.7,0.02),)
