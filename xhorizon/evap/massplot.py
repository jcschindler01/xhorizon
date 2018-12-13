
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import matplotlib.patches


def massplot(reglist, chainparams, params):
	"""
	"""
	#### params

	## ideal total mass
	M = 0.5 * params['Rmax']

	## horizon radii
	Rh = chainparams['Rh']

	## initialize masses
	m = 0. * Rh
	for i in range(len(reglist)):
		fp = reglist[i].metfunc.fparams
		if 'R' in fp.keys():
			m[i] = 0.5*fp['R']

	## junction locations
	r0, u0, v0 = chainparams['ps_r0'], chainparams['ps_u0'], chainparams['ps_v0']
	r1, u1, v1 = chainparams['fs_r0'], chainparams['fs_u0'], chainparams['fs_v0']


def massplot2():
	"""
	"""
	## setup figure #################
	## tex
	plt.rcParams['text.usetex'] = True
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.serif'] = []
	plt.rcParams['text.latex.preamble'] = []
	## rc
	plt.rcParams['font.size'] = 7
	plt.rcParams['axes.linewidth'] = .4
	plt.rcParams['xtick.major.size'] = 2
	plt.rcParams['ytick.major.size'] = 2
	## figure size
	fig_width  = 2.5  ## inches
	fig_aspect = 1.68   ## height/width
	fig_height = fig_width * fig_aspect
	## define figure and axes
	plt.figure(99, figsize=(fig_width,fig_height))
	## axes
	h, hh = .08, .33
	ax1 = plt.axes([.08,1.-h-hh,.89,hh]) # [left,bottom,width,height]
	ax2 = plt.axes([.08,h      ,.89,hh]) # [left,bottom,width,height]
	ax = [ax1,ax2]
	## format axes
	for axx in [ax1,ax2]:
		plt.sca(axx)
		## axes
		plt.xticks([0,1])
		plt.yticks([0,1])
		plt.ylabel(r"$m/M$", labelpad=-4)
		plt.xlim(-.3,1.3)
		plt.ylim(-.1,1.1)
		plt.grid(1, alpha=.15, zorder=0)
	## labels
	ax1.set_xlabel(r"$v/\tau_{f}$", labelpad=-6)
	ax2.set_xlabel(r"$u/\tau_{e}$", labelpad=-6)
	## labels
	ax1.set_title("$(d)$ Formation Dynamics")
	ax2.set_title("$(e)$ Evaporation Dynamics")
	################################


	####################### fake input #######################
	## mass
	m  = np.array([ 0.        ,  0.04921912,  0.07382869,  0.09843825,  0.09      ,  0.        ])

	## params
	params =    params =   {'A': 2.0, 'dv_evap': 0.5, 'Rmin': 0.18, 'B': 0.4, 'veta': 1.0, 'voff': -0.4, 'uoff': 0.0, 
							'Naccrete': 3, 'Rmax': 0.2, 'l': 0.1, 'ueta': 0.0, 'fparams0': {}, 'fparams1': {'l': 0.01}, 
						   }

	## chainparams
	chainparams =   {'Rh': np.array([ 0.        ,  0.09740062,  0.14697381,  0.19636592,  0.17944098,  0.        ]),
					 'fs_matchmode': ['rv', 'rv', 'rv', 'rv', 'rv', 'rv'],
					 'fs_r0': np.array([ 8.03671267,  7.98307469,  7.86133028,  0.29636592,  0.27944098,
					        7.73488086]),
					 'fs_u0': np.array([-16.47342534, -16.92254034, -16.78550319,   0.43013935,
					         0.90914107, -14.46976172]),
					 'fs_v0': np.array([-0.4 , -0.1 ,  0.1 ,  0.75,  1.25,  1.  ]),
					 'i0': 0,
					 'ps_matchmode': ['rv', 'rv', 'rv', 'rv', 'rv', 'rv'],
					 'ps_r0': np.array([ 8.03671267,  8.03671267,  7.98307469,  7.86133028,  0.29636592,
					        0.27944098]),
					 'ps_u0': np.array([-16.47342534, -17.23115108, -17.23361662, -16.90990623,
					         0.3189892 ,   0.44111805]),
					 'ps_v0': np.array([-0.4 , -0.3 , -0.1 ,  0.25,  0.75,  1.  ])
					 }	
#######################################################################################

	"""
	Parameters.
	"""

	## ideal total mass
	M = 0.5 * params['Rmax']
	
	## ideal lifetimes
	Tacc  = params['B']
	Tevap = params['A']

	## horizon radii
	Rh = chainparams['Rh']

	## initialize masses
	m = 1.*m    ##change for real case

	## icrit=nacc and mcrit
	Nacc = params['Naccrete']
	icrit = np.argmax(m)
	mcrit = m[icrit]

	## junction locations
	r0, u0, v0 = chainparams['ps_r0'], chainparams['ps_u0'], chainparams['ps_v0']
	r1, u1, v1 = chainparams['fs_r0'], chainparams['fs_u0'], chainparams['fs_v0']

	## inf
	npoints = 5001
	inf = 100.

	"""
	Accretion.
	v is split into Nacc+1 intervals, each of length dv
	    vv = past starting point of each interval
	    mm = mass in each interval
	"""

	## init
	vacc  = 1.*np.zeros(Nacc+1)
	macc  = 1.*np.zeros(Nacc+1)

	## fill vacc values
	vacc[0] = -1.*inf
	dvacc = (v1-v0)[1:Nacc]
	vacc[1:] = np.cumsum(np.concatenate([np.array([0.]),dvacc]))

	## fill macc values
	macc = 1.*m[:Nacc+1]

	## reference vacc arrays
	vv = Tacc * np.linspace(-1, 2, npoints)

	## mass step function
	mmacc = 0.*vv
	for i in range(Nacc+1):
		mmacc[vv>vacc[i]] = macc[i]

	## ideal mass function
	x = 1.*vv/Tacc
	macc_ideal = 0.*vv
	macc_ideal[x>0.] = M*(0.5 + 0.5*vv[x>0.]/Tacc)
	macc_ideal[x>1.] = M + 0.*vv[x>1.]

	"""
	Evap.
	u is split into Nevap+1 intervals, each of length du
	    uu = past ending point of each interval
	    mm = mass in each interval
	"""

	## Nevap
	Nevap = len(m[icrit:])-1

	## init
	uevap  = 1.*np.zeros(Nevap+1)
	mevap  = 1.*np.zeros(Nevap+1)

	## fill uevap values
	uevap[0] = -1.*inf
	duevap = (u1-u0)[-Nevap:-1]
	uevap[1:] = np.cumsum(np.concatenate([ np.array([0.]), duevap]))

	## fill mevap values
	mevap = 1.*m[-Nevap-1:]

	## reference uevap arrays
	uu = Tevap * np.linspace(-1, 2, npoints)

	## mass step function
	mmevap = 0.*uu
	for i in range(Nevap+1):
		mmevap[uu>uevap[i]] = mevap[i]

	## ideal mass function
	y = 1.*uu/Tevap
	mevap_ideal = 1.*M + 0.*uu
	mevap_ideal[y>0.] = M * (1. - (uu[y>0.]/Tevap))**(1./3.)
	mevap_ideal[y>1.] = 0.*uu[y>1.]



	"""
	Plot.
	"""



	## style
	fill_sty = dict(alpha=.15, edgecolor='none', zorder=100)
	border_sty = dict(alpha=.3, lw=1, zorder=150)
	ideal_sty = dict(alpha=1., lw=1.5, ls='--')

	## ax1
	plt.sca(ax1)

	## plot accretion
	plt.fill_between(vv/Tacc, 0., mmacc/M, color='b', **fill_sty)
	plt.plot(vv/Tacc, mmacc/M, color='b', **border_sty)
	plt.plot(vv/Tacc, macc_ideal/M, color='b', **ideal_sty)

	## legend
	blue1, blue2 = (0.,0.,1.,.15), (0.,0.,1.,.4)
	a1, = plt.plot([], [], label=r"$m_v\!(v)$ (ideal)", color='b', linestyle="dashed", dashes=(3,3), lw=1.5)
	a2 = matplotlib.patches.Patch(label=r"$m_v\!(v)$ (numerical)", color=blue1, ec=blue2, lw=.7)

	## make legend
	leg1 = plt.legend(handles=[a1,a2], loc='lower right', numpoints=10, fontsize=6)
	leg1.zorder = 900
	leg1.get_frame().set_linewidth(.4)

	## ax2
	plt.sca(ax2)

	## plot evaportion
	plt.fill_between(uu/Tevap, 0., mmevap/M, color='r', **fill_sty)
	plt.plot(uu/Tevap, mmevap/M, color='r', **border_sty)
	plt.plot(uu/Tevap, mevap_ideal/M, color='r', **ideal_sty)

	## legend
	red1, red2 = (1.,0.,0.,.15), (1.,0.,0.,.4)
	a1, = plt.plot([], [], label=r"$m_u\!(u)$ ideal", color='r', linestyle="dashed", dashes=(3,3), lw=1.5)
	a2 = matplotlib.patches.Patch(label=r"$m_u\!(u)$ numerical", color=red1, ec=red2, lw=.7)

	## make legend
	leg1 = plt.legend(handles=[a1,a2], loc='lower left', numpoints=10, fontsize=6)
	leg1.zorder = 900
	leg1.get_frame().set_linewidth(.4)


	## save
	if True:
		plt.savefig("temp.png", dpi=800)


	## show
	if False:
		plt.show()




############################################################





if __name__=='__main__':
	if False:
		import xhorizon as xh
		xh.evap.demo()
	else:
		massplot2()





