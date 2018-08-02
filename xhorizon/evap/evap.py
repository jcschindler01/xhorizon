
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

def evaporation(reg0, R=[], du=[], dv=[], u0=0., v0=0., l=0.1, rparams={}):
	"""
	"""
	##
	print "\nBEGIN EVAPORATION"
	## print announcement
	print "\n"
	print "EVAPORATION"
	print "reg0 = %r"%(reg0)
	print "R = %r"%(R)
	print "du = %r"%(du)
	print "dv = %r"%(dv)
	print "u0 = %r"%(u0)
	print "v0 = %r"%(v0)
	print "l = %r"%(l)
	print "rparams = %r"%(rparams)
	## input checker
	if not evaporation_input_good(R, du, dv):
		print "INPUT ERROR"
		return None
	## iteration values
	ii = range(len(R))
	## reglist
	reglist = [reg0]
	## evap steps
	for i in ii:
		## print announcement
		print "\nNEXT EVAP STEP, i=%r, R=%r, du=%r, dv=%r"%(i, R[i], du[i], dv[i])
		## get reg to build from
		reg1 = reglist.pop()
		## reg1 params tell you where in the existing region shell starts
		## get reg1 params for first step
		if i==0:
			u1, v1 = u0, v0
		## get reg1 params for subsequent steps
		if i>0:
			print "\nreg1.blocks[-1].uvbounds = %r"%(reg1.blocks[-1].uvbounds)
			u1 = reg1.blocks[-1].uvbounds['umin'] + du[i-1]
			v1 = reg1.blocks[-1].uvbounds['vmin'] + dv[i-1]
		## round reg1 params
		uvlen = 6
		#u1, v1 = np.round(u1,uvlen), np.round(v1,uvlen)
		## reg2 params tell you where next region starts
		## get reg2 params
		u2, R2 = u1, R[i]  ## could choose u2=0 if desired
		## print
		print "\n"
		print "u1, v1     = %r, %r"%(u1, v1)
		print "u2, R2, l2 = %r, %r, %r"%(u2, R2, l)
		print "rparams2 = %r"%(rparams)
		## do evap
		print "\nreglist = [length=%r]"%(len(reglist))
		reglist += evap(reg1, u1=u1, v1=v1, u2=u2, R2=R2, L2=l, rparams2=rparams)
	## final minkowski region
	## print announcement
	print "\nFINAL EVAP STEP"
	## get reg to build from
	reg1 = reglist.pop()
	## get reg1 params
	print "\nreg1.blocks[-1].uvbounds = %r"%(reg1.blocks[-1].uvbounds)
	u1 = reg1.blocks[-1].uvbounds['umin'] + du[-1]
	v1 = reg1.blocks[-1].uvbounds['vmin'] + dv[-1]
	## round reg1 params
	uvlen = 6
	u1, v1 = np.round(u1,uvlen), np.round(v1,uvlen)
	## get reg2 params
	u2 = u1 
	## print
	print "\n"
	print "u1, v1     = %r, %r"%(u1, v1)
	print "u2 = %r"%(u2)
	print "rparams2 = %r"%(rparams)
	## do evap
	print "\nreglist = [length=%r]"%(len(reglist))
	reglist += evap(reg1, func2=xh.mf.minkowski(), u1=u1, v1=v1, u2=u2, R2=None, L2=None, rparams2=rparams)
	##
	print "\nEND EVAPORTION"
	## return
	return reglist


############ individual evap steps #################

def evap(reg1, func2=None, u1=0., v1=0., u2=0., R2=1., L2=0.1, rparams2={}, rlines=False, boundary=False):
	"""
	Given reg1 = a Hayward region.
	Create an evaporation event and adjoining Hayward region reg2.
	Provide func2 = minkowski() for final step.
	Shell at:
		u1,v1,r1 in reg1
		u2,v2,r2 in reg2
	Provide input:
		reg1 = past region
		u1,v1 = coords in reg1
		u2 = u coord in reg2
	Calculate:
		r1 = r1(u1,v1)
		r2 = r1
		v2 = v2(u2,r2)
	"""
	print "\nBEGIN EVAP"
	## create new hayward metfunc unless otherwise provided
	if func2==None:
		func2 = xh.mf.hayward(R=R2,l=L2)
	## create new region
	reg2 = xh.reg.EFreg(func2, rparams=rparams2, boundary=boundary, rlines=rlines)
	## print starting info
	print "\n"
	print "REG1"
	print reg1.metfunc.info
	print reg1.rparams
	print "REG2"
	print reg2.metfunc.info
	print reg2.rparams
	## get r1 value for corner junction
	r1 = reg1.blocks[-1].r_of_uv(np.array([[u1],[v1]]))[0]
	#r1 = np.round(r1, 6)
	r1a = r1.copy()
	print "\n"
	print "r1a = %r"%r1a
	print "r1  = %r"%r1
	## passive slice of reg1
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), u0=u1, r0=r1)
	## get reg1 uvr values
	u1, v1, r1 = pslice.u0, pslice.v0, pslice.r0
	print "r1a = %r"%r1a
	print "r1  = %r"%r1
	## print reg1 params
	print "\n"
	print pslice
	print "u1, v1, r1 = %r, %r, %r"%(u1, v1, r1)
	## set U0(r) and V0(r) for target coords
	U0 = lambda r: pslice.U_of_r_at_v0(r)
	V0 = lambda r: pslice.V_of_r_at_u0(r)
	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), u0=u2, r0=r1, U0=U0, V0=V0)
	## get reg2 uvr values
	u2, v2, r2 = aslice.u0, aslice.v0, aslice.r0
	print "r1a = %r"%r1a
	print "r2  = %r"%r2
	## print reg2 params
	print "\n"
	print aslice
	print "u2, v2, r2 = %r, %r, %r"%(u2, v2, r2)
	print "r1a = %r"%r1a
	print "r2  = %r"%r2
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
	print "\n"
	print "EVAP --- NEW REGION AND SHELL"
	print "R2  = %r"%(R2)
	print "l2  = %r"%(L2)
	print "r1, u1, v1 = %r, %r, %r"%(r1, u1, v1)
	print "r2, u2, v2 = %r, %r, %r"%(r2, u2, v2)
	print " reg1.blocks[-1].uvbounds = %s"%( reg1.blocks[-1].uvbounds)
	print "reg1b.blocks[-1].uvbounds = %s"%(reg1b.blocks[-1].uvbounds)
	print " reg2.blocks[-1].uvbounds = %s"%( reg2.blocks[-1].uvbounds)
	##
	print "\nEND EVAP\n"
	## return
	return [reg1b, reg1, reg2]



##################################################################################################################



"""
Tests.
Run if __name__='__main__'.
"""



if __name__=='__main__':
	pass
	#test4()



