
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
Evaporation.
From highest to lowest level functions.
"""

###############################################################################################################3

######## chain of evap steps with cap #############

def evaporation(R, u, reg0, l=0.1, rparams={}):
	"""
	Starts with region reg0, evaporates away to minkowski space, with intermediate masses m at times u.

	Last m should be zero, if not zero it will be replaced by zero and ignored.

	Inputs:
		reg0 = initial hayward region
		u = 1d length n array of u0 values for accretion disks
		R= 2m = 1d length n array for series of masses for the hayward spaces

	Returns:
		reglist = list of regions
	"""
	## check for valid input
	go = True
	if not len(R)==len(u):
		print "\nERROR: Mismatched input array lengths.\n"
		go = False
	if not R[-1]==0.:
		print "\nERROR: Final radius value must be zero.\n"
		go = False
	if go==False:
		return None
	## reglist
	reglist = [reg0]
	## iteration values
	ivals = range(len(u))
	## evap steps
	for i in ivals[:-1]:
		ux, Rx = u[i], R[i]
		print "NEXT EVAP STEP, i=%r, u=%r, R=%r"%(i,ux,Rx)
		reglist += xh.evap.evap(reglist.pop(), u1=ux, u2=ux, R2=Rx, L2=l, rlines=False, boundary=False, rparams2=rparams)
	## final minkowski region
	for i in ivals[-1:]:
		ux = u[i]
		print "FINAL EVAP STEP, i=%r, u=%r, R=%r"%(i,ux,0.)
		reglist += xh.evap.evap_final(reglist.pop(), u1=6., u2=6., rlines=False, boundary=False, rparams2=rparams)
	## return
	return reglist


############ individual evap steps #################

def evap(reg1, r0=np.nan, u1=0., u2=0., R2=1., L2=0.1, rparams2={}, rlines=False, boundary=False):
	"""
	Given reg1 = a Hayward region.
	Create an evaporation event at corner radius r0.
	Create a new Hayward region reg2 with outer radius R2 and inner radius L2, using rparams2 and other options.
	Match to r0,u1 in reg1 and to r0,u2 in reg2.
	Must split reg1 into two new regions.
	Return [reg1a,reg1b,reg2] = list of new regions to be incorporated into diagram.
	Intended use looks like reglist += evap(args). 
	"""
	print "BEGIN EVAP"
	## create new region
	func2 = xh.mf.hayward(R=R2,l=L2)
	reg2 = xh.reg.EFreg(func2, rparams=rparams2, boundary=boundary, rlines=rlines)
	## print starting info
	print "REG1"
	print reg1.metfunc.info
	print reg1.rparams
	print "REG2"
	print reg2.metfunc.info
	print reg2.rparams
	## choose r0 value for corner junction
	r0 = get_rhawk_u0([reg1,reg2], u0=[u1,u2])
	print "r0 = %r"%r0
	## passive slice of reg1
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), u0=u1, r0=r0)
	print pslice
	v1 = 1.*pslice.v0
	print "v1 = %r"%v1
	## warn if bad u0 or v0 value
	#pslice, reg1 = xh.junc.slicecheck(pslice, reg1)
	## set U0(r) and V0(r) for target coords
	U0 = lambda r: pslice.U_of_r_at_v0(r)
	V0 = lambda r: pslice.V_of_r_at_u0(r)
	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[2], vblocks=[0,1,2], u0=u2, r0=r0, U0=U0, V0=V0)
	v2 = 1.*aslice.v0
	## warn if bad u0 or v0 value
	#aslice, reg2 = xh.junc.slicecheck(aslice, reg2)
	## update coordinate transformations
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	## create another copy of region 1
	reg1b = copy.deepcopy(reg1)
	reg1b.blocks = [reg1b.blocks[2]]
	## edit uvbounds
	## reg1
	for b in reg1.blocks:
		b.uvbounds.update(dict(vmax=v1))
	## reg1b
	for b in reg1b.blocks:
		b.uvbounds.update(dict(vmin=v1,umax=u1))
	## reg2
	for b in reg2.blocks:
		b.uvbounds.update(dict(vmin=v2))
	for b in [reg2.blocks[2]]:
		b.uvbounds.update(dict(umin=u2))
	## print announcement
	print "EVAP --- NEW REGION AND SHELL"
	print "R  = %r"%(R2)
	print "l  = %r"%(L2)
	print "r0 = %r"%(r0)
	print "u1, u2 = %r, %r"%(u1,u2)
	print "v1, v2 = %r, %r"%(v1,v2)
	print "\n"
	## return
	return [reg1b, reg1, reg2]


############ final evap step #########################

def evap_final(reg1, r0=np.nan, u1=0., u2=0., rparams2={}, rlines=False, boundary=False):
	"""
	Final evaporation step.
	Given reg1 = a Hayward region.
	Create an evaporation event at corner radius r0.
	Create a new MInkowski region reg2.
	Match to r0,u1 in reg1 and to r0,u2 in reg2.
	Must split reg1 into two new regions.
	Return [reg1a,reg1b,reg2] = list of new regions to be incorporated into diagram.
	Intended use looks like reglist += evap(args). 
	"""
	## create new region
	func2 = xh.mf.minkowski()
	reg2 = xh.reg.EFreg(func2, rparams=rparams2, boundary=boundary, rlines=rlines)
	## choose r0 value for corner junction
	r0 = get_rhawk_u0([reg1], u0=[u1])
	## passive slice of reg1
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), u0=u1, r0=r0)
	v1 = pslice.v0
	## warn if bad u0 or v0 value
	pslice, reg1 = xh.junc.slicecheck(pslice, reg1)
	## set U0(r) and V0(r) for target coords
	U0 = lambda r: pslice.U_of_r_at_v0(r)
	V0 = lambda r: pslice.V_of_r_at_u0(r)
	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[0], vblocks=[0], u0=u2, r0=r0, U0=U0, V0=V0)
	v2 = aslice.v0
	## warn if bad u0 or v0 value
	aslice, reg2 = xh.junc.slicecheck(aslice, reg2)
	## update coordinate transformations
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	## create another copy of region 1
	reg1b = copy.deepcopy(reg1)
	reg1b.blocks = [reg1b.blocks[2]]
	## edit uvbounds
	## reg1
	for b in reg1.blocks:
		b.uvbounds.update(dict(vmax=v1))
	## reg1b
	for b in reg1b.blocks:
		b.uvbounds.update(dict(vmin=v1,umax=u1))
	## reg2
	for b in reg2.blocks:
		b.uvbounds.update(dict(vmin=v2))
	for b in [reg2.blocks[-1]]:
		b.uvbounds.update(dict(umin=u2))
	## print announcement
	print "FINAL_EVAP --- NEW REGION AND SHELL"
	print "R  = %6.3f"%(0.)
	print "l  = %6.3f"%(0.)
	print "r0 = %6.3f"%(r0)
	print "u1, u2 = %6.3f, %6.3f"%(u1,u2)
	print "v1, v2 = %6.3f, %6.3f"%(v1,v2)
	print "\n"
	## return
	return [reg1b, reg1, reg2]


##################################################################################################################



"""
Tests.
Run if __name__='__main__'.
"""



if __name__=='__main__':
	pass
	#test5()



