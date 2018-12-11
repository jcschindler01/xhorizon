
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
	fig_width  = 2.  ## inches
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
	

	## params
	params = dict()
	## funcs
	#
	params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.hayward, fparams1=dict(l=.05)))
	#params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.schwarzschild, fparams1=dict()))
	## evap
	params.update(dict(Rmin=1., Rmax=1., dv_evap=.5, l=.05, A=.2))
	## accrete
	params.update(dict(B=.5, Naccrete=1))
	## offset
	params.update(dict(voff=0., veta=1., uoff=0., ueta=0.))

	## seed
	seed = 0

	## go
	reglist, chainparams = xh.evap.create_evap(params, seed=seed)

	## draw
	print("plot")
	xh.evap.drawreg(reglist, chainparams)


	if seed not in [0,-1]:
		plt.xlim(-3,3)
		plt.ylim(-1.5,4.5)

	## save
	if True:
		path = "temp-figs/demo"
		sfp = dict(dpi=800)
		temp_only = False
		xh.evap.evapsave(path=path, params=params, chainparams=chainparams, seed=seed, sfp=sfp, temp_only=temp_only)

	## show
	if False:
		plt.show()

	## return
	return reglist, chainparams




if __name__=="__main__":
	demo()


