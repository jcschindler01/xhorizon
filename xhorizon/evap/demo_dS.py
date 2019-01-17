
def demo():
	import numpy as np
	import matplotlib.pyplot as plt
	import xhorizon as xh
	from massplot import massplot

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
		sq = 6.
		x0 = -2.5
		y0 = 1.
		plt.xlim(x0,x0+sq)
		plt.ylim(y0-0.5*aspect*sq,y0+0.5*aspect*sq)
	################################
	

	## save
	path = "temp-figs/demo"
	sfp = dict(dpi=800)
	temp_only = True

	## draw
	draw = False

	## label
	label = '$(b)$'

	## input	
	L = 10.
	l  = .01
	le = .01

	## input
	R  = .2


	## accrete
	Tacc = .75


	## vv
	voff = 0. #-Tacc - 0.

	## params
	seed = 0
	mpgo=True
	params = dict()
	chainparams = dict()

##################################
	## go
	# func0 = xh.mf.dS(L=1.*L)
	# func1 = xh.mf.Hay_dS(l=1.*l, R=1.*R, L=1.*L)
	# func2 = xh.mf.Hay_dS_cutoff(l=1.*l, R=1.*R, L=1.*L)
	# func3 = xh.mf.dS(L=1.*L)
	# funclist = [func0, func1, func2, func3]



	# reg0 = xh.reg.MAXreg(func0, boundary=True, rlines=True)
	# reg1 = xh.reg.EFreg( func1, boundary=True, rlines=True)
	# reg2 = xh.reg.MAXreg(func2, boundary=True, rlines=True)
	# reg3 = xh.reg.MAXreg(func3, boundary=True, rlines=True)
	# reglist = [reg0, reg1, reg2, reg3]

	# for reg in reglist[2:3]:
	# 	reg.rplot()

	func2 = xh.mf.Hay_dS_cutoff(l=1.*l, R=1.*R, L=1.*L)

	print type(func2.rj[-1])
	print func2.f(func2.rj[-1])

######################################

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
		plabel = [r"$l_{ev}=%3s$"%(le), r"", r"$l=%3s$"%(l), r"$2M=%3s$"%(R), r"$L=%3s$"%(L)]
		#plabel += [r"", r"$\tau_{acc}=%3s$"%(Tacc), r"$\tau_{ev}=%3s$"%(Tevap)]
		plabel = "\n".join(plabel)
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)

	## manual param label
	if False:
		plabel = [r"$\tau_{ev}=%s$"%(Tevap), r"$\tau_{acc}=%s$"%(Tacc)]
		plabel += [r"", r"$M=10^{-3}$", r"$l_{ev}=l=10^{-4}$"]
		plabel = "\n".join(plabel)
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)

	## save
	if True:
		xh.evap.evapsave(path=path, params=params, chainparams=chainparams, seed=seed, sfp=sfp, temp_only=temp_only, massplot=mpgo)


	## show
	if False:
		plt.show()


	## return
	return None




if __name__=="__main__":
	demo()


