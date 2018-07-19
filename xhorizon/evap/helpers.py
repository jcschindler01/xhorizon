"""
This module contains helpers to help with evap code.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
import xhorizon as xh

################ tools for choosing junction corner radius ###############

def get_rinf_uv0(reglist, v0=[], u0=[]):
	"""
	Given a list of regions, assume that `outermost' region has index [-1].
	Assume r->inf as u->-inf.
	Find the r0 value corresponding to (v=v0 and u=2.*s0) in each region.
	Take minimum r0 over reglist.
	Return a value slightly less than r0.
	"""
	r0 = np.nan * np.ones(len(reglist))
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in [reg.blocks[-1]]:
			u = np.array([-2.*reg.rparams['s0']])
			v = np.array([v0[i]])
			uv = np.array([u,v])
			r0[i] = b.tr_of_uv(uv)[1]
	r0 = .95*np.min(r0)
	return r0


def get_rhawk_u0(reglist, u0=[]):
	"""
	Given a list of regions, assume that `outermost' region has index [-1].
	Assume r->inf as u->-inf.
	Find the r0 value corresponding to (u=v0 and u=2.*s0) in each region.
	Take maximum r0 over reglist.
	Return a value slightly more than r0.
	"""
	r0 = np.nan * np.ones(len(reglist))
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in [reg.blocks[-1]]:
			u = np.array([u0[i]])
			v = np.array([-2.*reg.rparams['s0']])
			uv = np.array([u,v])
			r0[i] = b.tr_of_uv(uv)[1]
	r0 = 1.05*np.max(r0)
	print r0
	return r0


###############################################################################
				

############# tools to rescale all coordinates at the end ###################

def UVcompose1(reg, fU=None, fV=None):
	"""
	Take the final UV coords for region and compose with a function in each null direction.
	"""
	reg2 = copy.deepcopy(reg)
	if not fU==None:
		reg.U_of_udl = lambda x: fU(reg2.U_of_udl(x))
	if not fV==None:
		reg.V_of_vdl = lambda x: fU(reg2.V_of_vdl(x))
	return reg

def UVcompose(reglist, fU=None, fV=None):
	"""
	Take the final UV coords for region and compose with a function in each null direction.
	Takes list input.
	"""
	for reg in reglist:
		reg = UVcompose1(reg, fU=fU, fV=fV)
	return reglist

############################################################################


################# tools for drawing ########################################3

def boundarylines(reglist, npoints=5001, sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			style = dict(lw=0.9, c='0.5', zorder=2000)
			style.update(sty)
			b.add_curves_uv(xh.cm.block_boundary(b, sty=style, npoints=npoints))


def colorlines(reglist, rmin=0.05, rmax=5., dr=.2, sty={}, npoints=2001, inf=25.):
	"""
	Add colorscaled lines of constant radius to region.
	Useful to check for matching.
	"""
	rvals = np.arange(rmin,rmax,dr)
	colvals = (rvals - np.min(rvals)) / (np.max(rvals) - np.min(rvals))
	cols = plt.cm.hsv_r(colvals)
	for reg in reglist:
		for b in reg.blocks:
			for i in range(len(rvals)):
				style = dict(c=cols[i], lw=.5, zorder=1500)
				style.update(sty)
				b.add_curves_tr(xh.cm.rlines([rvals[i]], sty=style, npoints=npoints, inf=inf))
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
	print "Rvals = %s"%(Rvals)
	#print "colvals = %s"%(colvals)
	## return
	return colvals


def fill_by_R(reglist, cm=plt.cm.prism):
	colvals = fillcols_by_R(reglist)
	for i in range(len(reglist)):
		reg = reglist[i]
		col = cm(colvals[i])
		for b in reg.blocks:
			sty = dict(fc=col, alpha=.2)
			b.fill(sty=sty)


############################################################################################

