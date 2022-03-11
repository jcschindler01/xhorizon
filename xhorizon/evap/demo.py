
def demo():
	import numpy as np
	import matplotlib.pyplot as plt
	import xhorizon as xh
	from xhorizon.evap.massplot import massplot

	"""
	"""

	## setup figure #################
	## rc
	plt.rcParams['font.size'] = 7
	plt.rcParams['axes.linewidth'] = .4
	plt.rcParams['xtick.major.size'] = 2
	plt.rcParams['ytick.major.size'] = 2
	plt.rcParams['lines.solid_capstyle'] = 'butt'
	plt.rcParams['lines.dash_capstyle'] = 'butt'
	## figure size
	fig_width  = 2.5  ## inches
	fig_aspect = 1.68   ## height/width
	fig_height = fig_width * fig_aspect
	## axis size
	width  = 0.86   ## x fraction
	aspect = 1.84   ## absolute height/width
	height = aspect * width / fig_aspect    ## y fraction
	## axis loc
	left, bottom = .12, .05     ## x fraction, y fraction
	## define figure and axes
	plt.figure(1,figsize=(fig_width,fig_height))
	ax1 = plt.axes([left,bottom,width,height], aspect=1.)
	ax = [ax1]
	## format axes
	for axx in ax:
		## set axis
		plt.sca(axx)
		## labels
		plt.xlabel('$V-U$', labelpad=-7)
		plt.ylabel('$V+U$', labelpad=-4)
		## ticks
		plt.xticks([1,2])
		plt.yticks([-1,0,1])
		## lims
		sq = 1.14 + .4
		x0 = .93 -.4
		plt.xlim(x0,x0+sq)
		plt.ylim(-0.5*aspect*sq,0.5*aspect*sq)
	################################
	

	## save
	path = "temp-figs/demo"
	sfp = dict(dpi=800)
	temp_only = False

	## draw
	draw = True

	## label
	label = '$(a)$'

	# func type
	ftype = 2

	## input
	l  = 1e-4
	le = 1e-4

	## input
	R  = 2e-3

	## if background
	L = 10.

	## evap
	Nevap = 12
	Tevap = 10.

	## accrete
	Nacc = 5
	Tacc = .5

	## seed
	seed = 0

	## vv
	voff = -Tacc - 0.
	veta = 1.

	## uu
	uoff = 0.
	ueta = 0.


	## params
	params = dict()
	## funcs
	#### schwarz
	if ftype==0:
		params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.schwarzschild, fparams1=dict()))
	#### hayward
	if ftype==1:
		params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.hayward, fparams1=dict(l=1.*l)))
	#### AdS
	if ftype==2:
		params.update(dict(functype0=xh.mf.AdS, fparams0=dict(L=1.*L), functype1=xh.mf.Hay_AdS, fparams1=dict(l=1.*l, L=1.*L)))
	## evap
	params.update(dict(Rmax=1.*R, le=1.*le, Tevap=1.*Tevap, Nevap=1*Nevap))
	## accrete
	params.update(dict(Tacc=1.*Tacc, Naccrete=1*Nacc))
	## offset
	params.update(dict(voff=1.*voff, veta=1.*veta, uoff=1.*uoff, ueta=1.*ueta))


##########################################
	## change lim if seed
	if seed not in [0,-1]:
		## format axes
		for axx in ax:
			## set axis
			plt.sca(axx)
			## ticks
			plt.xticks([-1,0,1,2,3])
			plt.yticks([-1,0,1,2,3])
			## lims
			sq = 3.4
			x0, y0 = -1.2, 1.2
			plt.xlim(x0,x0+sq)
			plt.ylim(y0-0.5*aspect*sq,y0+0.5*aspect*sq)
##########################################


	## go
	reglist, chainparams = xh.evap.create_evap(params, seed=seed)

	## draw
	if draw:
		print("plot")
		pp = dict(l=1.*l, R=1.*R)
		xh.evap.drawreg(reglist, chainparams, fparams=pp)


	## label
	if True:
		plt.annotate(s=label, xy=(.95,.97), xycoords='axes fraction', ha='right', va='top', size=8)


	## auto param label
	if False:
		plabel = [r"$l_{ev}=%3s$"%(le), r"", r"$l=%3s$"%(l), r"$2M=%3s$"%(R)]
		plabel += [r"", r"$\tau_{acc}=%3s$"%(Tacc), r"$\tau_{ev}=%3s$"%(Tevap)]
		plabel = "\n".join(plabel)
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)

	## manual param label
	if False:
		plabel = [r"$\tau_{ev}=%s$"%(Tevap), r"$\tau_{acc}=%s$"%(Tacc)]
		plabel += [r"", r"$M=10^{-3}$", r"$l_{ev}=l=10^{-4}$"]
		plabel = "\n".join(plabel)
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)

	## mass plot
	mpgo = True
	if mpgo:
		massplot(chainparams, params)

	## save
	if True:
		xh.evap.evapsave(path=path, params=params, chainparams=chainparams, seed=seed, sfp=sfp, temp_only=temp_only, massplot=mpgo)


	## show
	if False:
		plt.show()


	## return
	return reglist, chainparams




if __name__=="__main__":
	demo()


