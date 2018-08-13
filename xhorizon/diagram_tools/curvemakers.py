
"""
This module contains functious for generating useful curves which can be added to a diagram.

DEFAULT ZORDER:
100 = FILL
200 = UV GRID
400 = BOUNDARY
600 = RLINES


"""

import numpy as np

from curve_class import curve
from block_masks import rstar_limits


def rlines(rvals, sty={}, inf=50., t0=0., npoints=1000):
	"""
	Produce lines of constant radius for r in rvals. Return list of curves.
	"""
	## default style
	style = dict(zorder=600)
	## go
	crvlist = []
	for rval in rvals:
		## create new curve
		crv = curve()
		## define t and r arrays
		t = t0 + np.linspace(-inf, inf, 2.*npoints+1)
		r = rval + 0.*t
		## add coords to curve
		crv.tr = np.array([t,r])
		## set style
		style.update(sty)
		crv.sty.update(style)
		## add crv to output list
		crvlist += [crv]
	## return
	return crvlist




def rstarlines(vals, c=None, sty={}, inf=50., npoints=1000):
	"""
	Produce lines of constant rstar for rstar in vals. Return list of curves.
	"""
	crvlist = []
	for val in vals:
		## create new curve
		crv = curve()
		## define t array
		t = np.linspace(-inf, inf, 2.*npoints+1)
		## define u and v arrays
		c = float(c)
		u = t - val + c
		v = t + val - c
		## add coords to curve
		crv.uv = np.array([u,v])
		## update curve style from default
		style = dict(zorder=600)
		style.update(sty)
		crv.sty.update(style)
		## add crv to output list
		crvlist += [crv]
	## return
	return crvlist




def uvlines(uvvals, uv='uv', sty={}, c=0., eps=1e-3, inf=100., npoints=1000):
	"""
	Produce lines of constant radius for r in rvals. Return list of curves.
	"""
	crvlist = []
	for val in uvvals:
		## create new curve
		ucrv = curve()
		vcrv = curve()
		## define linearly spaced array
		sa = np.linspace(-inf, inf, 2*npoints+1)
		sb = val - 2.*c + eps*np.array([-1.,0.,1.])
		sc = val + 2.*c + eps*np.array([-1.,0.,1.])
		s = np.sort(np.concatenate([sa,sb,sc]))
		## ucrv is line of const, vcrv is line of const v
		ucrv.uv = np.array([val+0.*s, 1.*s])
		vcrv.uv = np.array([1.*s, val+0.*s])
		## update curve style from default
		style = dict(zorder=200)
		style.update(sty)
		ucrv.sty.update(style)
		vcrv.sty.update(style)
		## add crvs to output list
		if 'u' in uv:
			crvlist += [ucrv]
		if 'v' in uv:
			crvlist += [vcrv]
	## return
	return crvlist


def block_boundary(blk, sty={}, inf=100., npoints=1000):
	"""
	"""
	## change default style
	style = dict(c='k', ls='-', lw=1.5, zorder=400)
	style.update(sty)
	## get min and max rstar values
	j = blk.j
	lims = rstar_limits(blk.master.metfunc)
	lims = lims[j:j+2]
	## adjust to avoid accidental masking
	idx = np.argsort(lims)
	lims[idx] = lims[idx] + 2. * np.array([1.,-1.]) * blk.master.metfunc.Fparams['eps']
	## make list of boundary curves
	crvlist = []
	## make curves
	for rstar in lims:
		if rstar==np.inf:
			rstar = inf
		if rstar==-np.inf:
			rstar = -inf
		crvlist += rstarlines([rstar], c=blk.master.rparams['c'], sty=style, inf=2.*inf, npoints=npoints)
	## return
	return crvlist


def block_boundary_2(blk, sty={}, inf=100., npoints=1000):
	"""
	"""
	## change default style
	style = dict(c='k', ls='-', lw=1.5, zorder=400)
	style.update(sty)
	## get min and max rstar values
	j = blk.j
	lims = rstar_limits(blk.master.metfunc)
	lims = lims[j:j+2]
	## adjust to avoid accidental masking
	idx = np.argsort(lims)
	lims[idx] = lims[idx] + 2. * np.array([1.,-1.]) * blk.master.metfunc.Fparams['eps']
	## make list of boundary curves
	crvlist = []
	## make curves
	for rstar in lims:
		if rstar==np.inf:
			rstar = inf
		if rstar==-np.inf:
			rstar = -inf
		crvlist += rstarlines_special_2([rstar], c=blk.master.rparams['c'], sty=style, inf=2.*inf, npoints=npoints)
	## return
	return crvlist




def default_rlines(func, sty={}, inf=50., npoints=1000):
	## get radius length scales of the problem
	rscales = (func.rj[1:-1])
	## no length cale
	if len(rscales) == 0:
		dr = np.array([0.5])
		rr = np.array([0.,20.])
	## one length scale
	if len(rscales) == 1:
		dr = rscales / 10.
		rr = np.array([0.5*dr[0], 20.*rscales[0]])
	## multiple length scales
	if len(rscales) >= 2:
		dr = rscales / 10.
		rr = np.zeros(len(rscales)+1)
		rr[1:] = 5.*rscales[:]
		rr[-1] = 20.*rscales[-1]
		rr[0] = 0.5*rscales[0]
	## create curves
	cols = ['midnightblue', 'sage']
	crvlist = []
	for i in range(len(rr)-1):
		imax = len(rr)-2
		col = cols[(imax-i)%len(cols)]
		rvals = np.arange(rr[i], rr[i+1], dr[i])
		style = dict(c=col, zorder=599-i)
		style.update(sty)
		crvlist += rlines(rvals, sty=style, inf=inf, npoints=npoints)
	## return
	return crvlist




def rstarlines_special_1(vals, c=None,  s0=None, ki=None, A=1e4, sty={}, inf=100., npoints1=1001, npoints2=1001, npoints3=1001):
	"""
	Produce lines of constant rstar for rstar in vals. Return list of curves.
	These are specially supplemented to give many data points immediately after s0 for zoom views.
	"""
	crvlist = []
	for val in vals:
		## initialize params
		c = float(c)
		## create new curve
		crv = curve()
		## initialize u and v lists
		ulist = []
		vlist = []
		## small values
		ss = np.linspace(-1.02*s0,1.02*s0,npoints1)
		for x in [2.*ss]:
			## u values
			ulist += [x]
			vlist += [x + 2.*(val-c)]
			## v values
			vlist += [-x]
			ulist += [-x - 2.*(val-c)]
		## large values
		ss = np.linspace(s0,inf,npoints3)
		for x in [2.*ss, -2.*ss]:
			## u values
			ulist += [x]
			vlist += [x + 2.*(val-c)]
			## v values
			vlist += [-x]
			ulist += [-x - 2.*(val-c)]
		## prep intermediate values
		e = np.sign(ki)
		ds = ki**(-1) * np.log(1.+e*ki*A)
		np2 = [ int(np.abs(e[i])) * npoints2 for i in range(len(ki)) ]
		ss = [ np.linspace(e[i]*s0, e[i]*s0+ds[i], np2[i]) for i in range(len(ki)) ]
		## intermediate values
		for x in [2.*ss[0], 2.*ss[1]]:
			## u values
			ulist += [x]
			vlist += [x + 2.*(val-c)]
			## v values
			vlist += [-x]
			ulist += [-x - 2.*(val-c)]
		## concatenate to form output
		u = np.concatenate(ulist)
		v = np.concatenate(vlist)
		## sort
		mask = np.argsort(u)
		u = u[mask]
		v = v[mask]
		## add coords to curve
		crv.uv = np.array([u,v])
		## update curve style from default
		style = dict(zorder=600)
		style.update(sty)
		crv.sty.update(style)
		## add crv to output list
		crvlist += [crv]
	## return
	return crvlist





def rstarlines_special_2(vals, c=None, sty={}, inf=50., npoints=1000):
	"""
	Produce lines of constant rstar for rstar in vals. Return list of curves.
	"""
	crvlist = []
	for val in vals:
		## create new curve
		crv = curve()
		## define t array
		t = np.linspace(-inf, inf, 2.*npoints+1)
		## define u and v arrays
		c = float(c)
		u = t - val + c
		v = t + val - c
		## add coords to curve
		crv.uv = np.array([u,v])
		## update curve style from default
		style = dict(zorder=600)
		style.update(sty)
		crv.sty.update(style)
		## add crv to output list
		crvlist += [crv]
	## return
	return crvlist





def colnorm1(svals,smin,smax):
	"""
	Color bar normalizer for input into colormap.
	Input an array of svals.
	Return an array of the absolute values normalized to (0,1) and truncated at smin, smax.
	"""
	## initialize output array to abs val of input
	ss = 1. * np.abs(svals)
	## cutoff at min
	ss[ss<smin] = smin + 0.*ss[ss<smin]
	## cutoff at max	
	ss[ss>smax] = smax + 0.*ss[ss>smax]
	## normalize
	ss = (ss - smin) / (smax - smin)
	## return
	return ss




def rstarlines_special_2():
	"""
	"""
	## return




