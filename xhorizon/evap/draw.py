"""
This module contains helpers to draw with evap code.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

import xhorizon as xh
from xhorizon.shell_junction import interpolators as interp


############## plot reglist ####################

def drawreg(reglist):
	"""
	Plot all regions in reglist.
	"""
	## lines
	rline_zero(reglist, sty=dict(), npoints=5001)
	rline_inf( reglist, sty=dict(), npoints=5001)
	rline_inf_col( reglist, sty=dict(), npoints=5001)
	rlines(reglist)
	#xh.evap.uv_lines(reglist, uv='u', sty=dict())
	#xh.evap.s0_lines(reglist, sty=dict(lw=3))
	## figure
	plt.figure(figsize=(6,6))
	plt.axes([.1, .1, .8, .8])
	plt.xlim(-3,3)
	plt.ylim(-3,3)
	plt.gca().set_aspect('equal')
	## plot
	for reg in reglist:
		reg.rplot()
	## fill
	#xh.evap.fill_by_R_2(reglist, col='r', amax=.2)

#############################################################


################# tools for drawing ########################################

def rline_zero(reglist, npoints=5001, sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			if b.j==0:
				## style
				style = dict(lw=0.9, c='0.5', zorder=2000)
				style.update(sty)
				## curve
				cv = xh.cm.rstarlines_special_2([0.], b.uvbounds, c=0., sty=style, inf=100., npoints=1001, eps=1e-12)
				b.add_curves_uv(cv)


def rline_inf(reglist, npoints=5001, inf=1000., sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			## outermost blocks only
			if b.j==len(b.master.metfunc.rj)-2:
				## style
				style = dict(lw=0.9, c='0.5', zorder=2010)
				style.update(sty)
				## rstar value
				rstar = b.master.metfunc.rstar_ref[-1]
				## curve
				#b.add_curves_uv(xh.cm.block_boundary_2(b, sty=style, inf=100., npoints=npoints)[-1:])
				b.add_curves_uv(xh.cm.uvlines([-inf,inf], uv='uv', uvbounds=b.uvbounds, sty=style, c=0., inf=inf, npoints=npoints))


def rline_inf_col(reglist, npoints=5001, sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			## outermost blocks only
			if b.j==len(b.master.metfunc.rj)-2:
				## style
				style = dict(lw=.9, c='0.5', zorder=2000)
				style.update(sty)
				## rstar value
				rstar = b.master.metfunc.rstar_ref[-1]
				## mass
				m = 0.
				if 'R' in b.master.metfunc.fparams.keys():
					m = 0.5 * b.master.metfunc.fparams['R']
				## color
				cmin = .3
				colval = cmin + (1.-cmin)*m
				style.update(dict(c=plt.cm.Blues(colval)))
				## curve
				b.add_curves_uv(xh.cm.block_boundary_2(b, sty=style, inf=100., npoints=npoints)[-1:])


def rlines(reglist, rmin=0.05, rmax=25., dr=.2, sty={}, cm=plt.cm.hsv_r, npoints=2001, inf=25.):
	"""
	Add colorscaled lines of constant radius to region.
	Useful to check for matching.
	"""
	rvals = np.arange(rmin,rmax,dr)
	colvals = (rvals - np.min(rvals)) / (np.max(rvals) - np.min(rvals))
	cols = cm(colvals)
	for reg in reglist:
		for b in reg.blocks:
			for i in range(len(rvals)):
				if (b.rj[0]<rvals[i]) and (rvals[i]<b.rj[1]):
					style = dict(c=cols[i], ls='-', markersize=1, lw=.5, zorder=1500)
					style.update(sty)
					rstarvals = 1.*b.master.metfunc.F(rvals[i:i+1])
					b.add_curves_uv(xh.cm.rstarlines_special_2(rstarvals, b.uvbounds, c=b.master.rparams['c'], sty=style, inf=2.*inf, npoints=npoints))
	return reglist


def s0_lines(reglist, sty={}, npoints=1001, inf=50., eps=1e-24):
	"""
	"""
	for reg in reglist:
		for b in reg.blocks:
			x = np.array([(2.*reg.rparams['s0'] - eps)])
			style1 = dict(c='m', alpha=.5, lw=3, ls=':', zorder=5000)
			style2 = dict(c='c', alpha=.5, lw=3, ls=':', zorder=5000)
			style1.update(sty)
			style2.update(sty)
			b.add_curves_uv(xh.cm.uvlines( x, uv='uv', uvbounds=b.uvbounds, sty=style1, c=0., inf=inf, npoints=npoints))
			b.add_curves_uv(xh.cm.uvlines(-x, uv='uv', uvbounds=b.uvbounds, sty=style2, c=0., inf=inf, npoints=npoints))
	return reglist


def uv_lines(reglist, uv='uv', sty={}, npoints=1001, inf=50., eps=1e-24):
	"""
	"""
	for reg in reglist:
		for b in reg.blocks:
			smin, smax, ds = -5., 25., 1.
			vals = np.arange(smin, smax, ds)
			cm = plt.cm.gist_rainbow
			cv = np.linspace(0,1,len(vals))
			for i in range(len(vals)):
				style1 = dict(c=cm(cv[i]), lw=.8, ls='-', zorder=6000)
				style1.update(sty)
				b.add_curves_uv(xh.cm.uvlines([vals[i]], uv=uv, uvbounds=b.uvbounds, sty=style1, c=0., inf=inf, npoints=npoints))
	return reglist

##########################################################################################




######################### tools for filling #############################################

def fillcols_by_R(reglist):
	"""
	Get fill color values based on radius.
	"""
	## get rvals
	Rvals = np.zeros(len(reglist))
	for i, reg in enumerate(reglist):
		if 'R' in reg.metfunc.fparams.keys():
			x = reg.metfunc.fparams['R']
			Rvals[i] = x
	## get colvals
	norm = np.max(Rvals)
	colvals = 0.2 + 0.7 * Rvals / norm
	#print "Rvals = %s"%(Rvals)
	#print "colvals = %s"%(colvals)
	## return
	return colvals

def alphas_by_R_2(reglist, amax=.9):
	"""
	Get fill color values based on radius.
	"""
	## rvals
	Rvals = np.zeros(len(reglist))
	for i, reg in enumerate(reglist):
		if 'R' in reg.metfunc.fparams.keys():
			x = reg.metfunc.fparams['R']
			Rvals[i] = 1.*x
	## log
	Rvals = np.log(1.+Rvals)
	## alphas
	amin, amax = 0., amax
	avals = amin + (amax-amin)*(Rvals/np.max(Rvals))
	## return
	return 1.*avals


def fill_by_R(reglist, cm=plt.cm.prism):
	colvals = fillcols_by_R(reglist)
	for i in range(len(reglist)):
		reg = reglist[i]
		col = cm(colvals[i])
		for b in reg.blocks:
			sty = dict(fc=col, alpha=.2)
			b.fill(sty=sty)

def fill_by_R_2(reglist, col='r', amax=.9):
	##
	alphas = alphas_by_R_2(reglist, amax=amax)
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in reg.blocks:
			sty = dict(fc=col, alpha=alphas[i])
			b.fill(sty=sty)

############################################################################################


################## more helpers #############



##########################################

































if __name__=='__main__':
	import os
	os.system("python C:\\Joe\\gdrive\\physics\\2018\\form_evap\\form_evap\\figures\\004\\004b.py")

