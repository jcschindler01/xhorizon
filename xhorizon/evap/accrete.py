
"""
This module provides method for making forming and evaporation BH diagrams.
This module imports the entire xhorizon package. 
It is meant for a higher level usage than the other subpackages, none of the
guts of xhorizon rely on this.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

import xhorizon as xh
from helpers import *


"""
Accretion.
From highest to lowest level functions.
"""

##################################################################################################################


######### chain of accretion steps with initial cap ##################3

def accretion(m, v, l=0.1, rparams={}):
	"""
	Stars with minkowski space, then adds a series of accretions as specified.

	Inputs:
		v = 1d length n array of v0 values for accretion disks
		m = 1d length n array for series of masses for the hayward spaces

	Returns:
		reglist = list of regions
	"""
	## initial minkowski region
	reglist = [xh.reg.EFreg(xh.mf.minkowski(),rlines=False, boundary=False,rparams=rparams)]
	## accrete
	for i in range(len(v)):
		vx, Rx = 1.*v[i], 2.*m[i]
		reglist += xh.evap.accrete(reglist.pop(), v1=vx, v2=vx, R2=Rx, L2=l, rlines=False, boundary=False, rparams2=rparams)
	## return
	return reglist


############### individual accretion steps ########################

def accrete(reg1, v1=0., v2=0., R2=1., L2=0.1, rparams2={}, rlines=False, boundary=False):
	"""
	Accrete a Hayward of outer radius R2 onto the region reg1 on a line of constant v by 
	slicing reg1 at v1 and the new region reg2 at v2.

	Use a corner junction at radius r0=rinf (as close to r->inf as possible withouth causing
	slice matching problems).

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
	## choose r0 value for corner junction
	r0 = get_rinf_uv0([reg1,reg2],v0=[v1,v2])
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






##################################################################################################################



"""
Tests.
Run if __name__='__main__'.
"""



if __name__=='__main__':
	pass
	#test5()


