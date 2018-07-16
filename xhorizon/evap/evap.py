
"""
This module provides method for making forming and evaporation BH diagrams.
This module imports the entire xhorizon package. 
It is meant for a higher level usage than the other subpackages, none of the
guts of xhorizon rely on this.
"""

import numpy as np
import matplotlib.pyplot as plt

import xhorizon as xh


def accrete(reg1, v1=0., v2=0., r0=3., R2=1., L2=0.1, rparams2={}, rlines=False, boundary=True):
	"""
	Accrete a Hayward of outer radius R2 onto the region reg1 on a line of constant v by 
	slicing reg1 at v1 and the new region reg2 at v2.

	Inputs:
		reg1 = existing region to accrete onto
		v1   = value of v to slice reg1
		v2   = value of v to slice new region
		R2   = Hayward outer radius of new region

	Returns:
		[reg1, reg2] where....
		reg1 = original region modified
		reg2 = new Hayward region

	Intended use is something like reglist += accrete(args).
	"""
	## create new region
	func2 = xh.mf.hayward(R=R2,l=L2)
	reg2 = xh.reg.EFreg(func2, rparams=rparams2, boundary=boundary, rlines=rlines)
	## edit uvbounds
	for b in reg1.blocks:
		b.uvbounds.update(dict(vmax=v1))
	for b in reg2.blocks:
		b.uvbounds.update(dict(vmin=v2))
	## passive slice of reg1
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), v0=v1, r0=r0)
	## warn if bad u0 or v0 value
	pslice, reg1 = xh.junc.slicecheck(pslice, reg1)
	## set U0(r) and V0(r) for target coords
	U0 = lambda r: pslice.U_of_r_at_v0(r)
	V0 = lambda r: pslice.V_of_r_at_u0(r)
	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[2], vblocks=[0,1,2], v0=v2, r0=r0, U0=U0, V0=V0)
	## warn if bad u0 or v0 value
	aslice, reg2 = xh.junc.slicecheck(aslice, reg2)
	## update coordinate transformations
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	## return
	return [reg1, reg2]





def colorlines(reglist, rmin=0.05, rmax=3., dr=.2, sty={}, npoints=2001, inf=25.):
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

