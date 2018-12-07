
def demo():
	import numpy as np
	import matplotlib.pyplot as plt
	import xhorizon as xh

	"""
	"""

	## params
	params = dict()
	## funcs
	#
	#params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.hayward, fparams1=dict(l=.05)))
	params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.schwarzschild, fparams1=dict()))
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
	print "plot"
	xh.evap.drawreg(reglist, chainparams)

	## format axes
	## ticks
	plt.xticks([1,2])
	plt.yticks([-1,0,1])
	## lims
	plt.xlim(.8,2.2)
	plt.ylim(-1.2,1.2)
	## labels
	plt.xlabel('$V-U$')
	plt.ylabel('$V+U$')

	if seed not in [0,-1]:
		plt.xlim(-3,3)
		plt.ylim(-1.5,4.5)

	## save
	if True:
		path = "temp-figs/demo"
		sfp = dict(dpi=400)
		temp_only = False
		xh.evap.evapsave(path=path, params=params, chainparams=chainparams, seed=seed, sfp=sfp, temp_only=temp_only)

	## show
	if False:
		plt.show()

	## return
	return reglist, chainparams




if __name__=="__main__":
	demo()


