import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt2d

import ocelot as oclt
from ocelot.test.cox_tests.cox_configs import *
from ocelot.test.cox_tests import cox_utils as cox

from ipywidgets import interactive
from ipywidgets import interact

method = oclt.MethodTM()
method.global_method = oclt.SecondTM

def get_envs(QAP1=1., QAP2=1., QAP3=1., \
             SRC=Drifts['LPWA-QAP1']*1e3,TARG=Drifts['QAP3-IMG1']*1e2, ):

	emin,emax = 0.04,0.25
	energies = np.r_[emin:emax:300j]
	mags = []

	for BeamEnergy in energies:
		indx_stop = np.nonzero(np.array(cell_keys)=='IMG1')[0][0]
		QuadGradients1 = oclt.deepcopy(QuadGradients)
		QuadGradients1['QAP1'] *= QAP1
		QuadGradients1['QAP2'] *= QAP2
		QuadGradients1['QAP3'] *= QAP3
		Drifts1 = oclt.deepcopy(Drifts)
		Drifts1['LPWA-QAP1'] = SRC*1e-3
		Drifts1['QAP3-IMG1'] = TARG*1e-2

		latt_elmts = cox.make_line( \
		  BeamEnergy=BeamEnergy,QuadGradients=QuadGradients1,Drifts=Drifts1 )
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
	plt.plot(energies*1e3,np.sqrt(sx*sy),c='k',lw=1.5)
	plt.legend(('$\sigma_x$','$\sigma_z$','$\sqrt{\sigma_x\sigma_z}$'),loc=1)
	plt.ylabel('Beam sizes (mm)')
	plt.xlabel('Electron energy (MeV)')
	plt.ylim(0,4.);

def plot_beam(v,beam,plot_XY=False,plot_spect=True,**imshowargs):
	p_arrays_init = cox.make_beam_sliced(beam, \
	  div_chirp=0.5*(beam['E'][0]+beam['E'][1]))

	if plot_spect:
		fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
	else:
		fig, ax1 = plt.subplots(1,1,figsize=(5,5))

	if plot_XY:
		fig2,(ax3,ax4) = plt.subplots(1,2,figsize=(10,4))

	p_arrays = oclt.deepcopy(p_arrays_init)

	QuadGradients1 = oclt.deepcopy(QuadGradients)
	QuadGradients1['QAP1'] *= v.children[0].value
	QuadGradients1['QAP2'] *= v.children[1].value
	QuadGradients1['QAP3'] *= v.children[2].value
	Drifts1 = oclt.deepcopy(Drifts)
	Drifts1['LPWA-QAP1'] = v.children[3].value*1e-3
	Drifts1['QAP3-IMG1'] = v.children[4].value*1e-2

	p_arrays,outputs = cox.make_shot(p_arrays, \
	  QuadGradients=QuadGradients1,Drifts=Drifts1,stop_key='IMG1',)
	p_arrays = cox.aligh_slices(p_arrays, \
	  QuadGradients=QuadGradients1,Drifts=Drifts1,stop_key='IMG1')

	xx = np.hstack(([p_array.x()*1e3 for p_array in p_arrays]))
	yy = np.hstack(([p_array.y()*1e3 for p_array in p_arrays]))
	zz = np.hstack(([-p_array.tau()*1e3 for p_array in p_arrays]))
	px = np.hstack(([p_array.px()*1e3 for p_array in p_arrays]))
	py = np.hstack(([p_array.py()*1e3 for p_array in p_arrays]))
	ee = np.hstack(([(p_array.p()+1)*p_array.E*1e3 for p_array in p_arrays]))
	rr = np.sqrt(xx**2+yy**2)

	h,x,y = np.histogram2d(xx,yy,bins=800,range=[[-0.1,0.1],[-0.1,0.1]])
	h = medfilt2d(np.abs(h),5)
	pl = ax1.imshow(h.T,aspect='auto', \
	  extent=(x.min(),x.max(),y.min(),y.max()),origin='lower',**imshowargs)

	ax1.set_xlabel('X (mm)')
	ax1.set_ylabel('Y (mm)')

	hx = h[:,h.shape[1]/2-1:h.shape[1]/2+2].mean(-1)
	hy = h[h.shape[0]/2-1:h.shape[0]/2+2,:].mean(0)

	hx -= hx.min()
	hat_ind = np.nonzero(hx>0.5*hx.max())[0]
	sx_fwhm = x[hat_ind[-1]] - x[hat_ind[0]]
	hy -= hy.min()
	hat_ind = np.nonzero(hy>0.5*hy.max())[0]
	sy_fwhm = y[hat_ind[-1]] - y[hat_ind[0]]

	if plot_spect:
		h,x,y = np.histogram2d(ee,rr,bins=400, \
		  range=[[20,250],[0,0.15]])
		h = medfilt2d(np.abs(h),5)
		pl = ax2.imshow(h.T,\
		  aspect='auto',extent=(x.min(),x.max(),y.min(),y.max()), \
		  origin='lower',**imshowargs)
		ax2.set_xlabel('Electron energy (MeV)')
		ax2.set_ylabel('Radius (mm)')

	if plot_XY:
		h,x,y = np.histogram2d(ee,xx,bins=400, \
		  range=[[20,250],[-3,3]])
		pl = ax3.imshow(medfilt2d(np.abs(h),5).T,\
		  aspect='auto',extent=(x.min(),x.max(),y.min(),y.max()), \
		  origin='lower',**imshowargs)

		h,x,y = np.histogram2d(ee,yy,bins=400, \
		  range=[[20,250],[-3,3]])
		pl = ax4.imshow(medfilt2d(np.abs(h),5).T,\
		  aspect='auto',extent=(x.min(),x.max(),y.min(),y.max()), \
		  origin='lower',**imshowargs)

		ax3.set_ylabel('X (mm)')
		ax4.set_ylabel('Y (mm)')
		for ax in (ax3,ax4,): ax.set_xlabel('Electron energy (MeV)')

	part_select_mask = (rr<100e-3)
	part_select_ind = np.nonzero(part_select_mask)[0]
	num_parts_full = rr.shape[0]
	num_parts_select = part_select_mask.sum()
	select_frac = 100.*num_parts_select/num_parts_full
	sx_rms = xx[part_select_ind].std()
	sy_rms = yy[part_select_ind].std()
	sz = zz[part_select_ind].std()
	ox = px[part_select_ind].std()
	oy = py[part_select_ind].std()
	e_cent = ee[part_select_ind].mean()
	se = ee[part_select_ind].std()/e_cent

	print('\nFraction in r<100 um {0:g}%'.format( \
	  select_frac )  )
	print( \
	  "Sx={0:g} ({1:g}), Sy={2:g} ({3:g}), Sz={4:g} [um] RMS (FWHM)\n" \
	  .format(sx_rms*1e3,sx_fwhm*1e3,sy_rms*1e3,sy_fwhm*1e3,sz*1e3) \
	  +"Ox={0:g}, Oy={1:g}" \
	  .format(ox,oy) \
	  + " [mrad]\ne={0:g} [MeV], de={1:g}%"\
	  .format(e_cent,se*100) )

	out_dat = np.vstack((xx,yy,zz,px,py,ee  ))
	return out_dat

qap_widgt = interactive(get_envs, \
  QAP1=(0.87,1.6,0.005), \
  QAP2=(0.88,1.6,0.005), \
  QAP3=(0.96,1.7,0.005), \
  SRC=(0,100,2), \
  TARG=(0,200,2.), \
)

qap_widgt.children[3].description = 'SRC (mm)'
qap_widgt.children[4].description = 'TARG (cm)'
