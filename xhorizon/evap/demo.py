
def demo():
	import numpy as np
	import matplotlib.pyplot as plt
	import xhorizon as xh

	"""
	"""

	## setup figure #################
	## rc
	plt.rcParams['font.size'] = 7
	plt.rcParams['axes.linewidth'] = .4
	plt.rcParams['xtick.major.size'] = 2
	plt.rcParams['ytick.major.size'] = 2
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
		sq = 1.14
		x0 = .93
		plt.xlim(x0,x0+sq)
		plt.ylim(-0.5*aspect*sq,0.5*aspect*sq)
	################################
	

	## save
	path = "temp-figs/demo"
	sfp = dict(dpi=800)
	temp_only = True

	## label
	label = ''

	# func type
	ftype = 1

	## input
	l = .05
	R = 1.

	## evap
	R0 = .8
	tauB = 1.
	dv = .75
	le = .005

	## accrete
	Nacc = 2
	tauA = 1.

	## seed
	seed=0

	## vv
	voff = 0.
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
	## evap
	params.update(dict(Rmin=1.*R0, Rmax=1.*R, dv_evap=1.*dv, l=1.*le, A=1.*tauB))
	## accrete
	params.update(dict(B=1.*tauB, Naccrete=1*Nacc))
	## offset
	params.update(dict(voff=1.*voff, veta=1.*veta, uoff=1.*uoff, ueta=1.*ueta))

	## go
	reglist, chainparams = xh.evap.create_evap(params, seed=seed)

	## draw
	print("plot")
	pp = dict(l=1.*l, R=1.*R)
	xh.evap.drawreg(reglist, chainparams, fparams=pp)


	if seed not in [0,-1]:
		plt.xlim(-3,3)
		plt.ylim(-1.5,4.5)

	## label
	if True:
		plt.annotate(s=label, xy=(.95,.97), xycoords='axes fraction', ha='right', va='top', size=8)


	## param label
	if True:
		plabel = "\n".join([r"$l_{ev}=%s$"%(le), r"", r"$l=%s$"%(l), r"$2M=%s$"%(R)])
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)



	## save
	if True:
		xh.evap.evapsave(path=path, params=params, chainparams=chainparams, seed=seed, sfp=sfp, temp_only=temp_only)


	## show
	if False:
		plt.show()


	## return
	return reglist, chainparams




if __name__=="__main__":
	demo()


