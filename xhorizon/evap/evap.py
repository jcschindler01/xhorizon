
"""
This module provides method for making forming and evaporation BH diagrams.
This module imports the entire xhorizon package. 
It is meant for a higher level usage than the other subpackages, none of the
guts of xhorizon rely on this.
"""

import numpy as np
import matplotlib.pyplot as plt

import xhorizon as xh




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
	## create new region
	func2 = xh.mf.hayward(R=R2,l=L2)
	reg2 = xh.reg.EFreg(func2, rparams=rparams2, boundary=boundary, rlines=rlines)
	## choose r0 value for corner junction
	r0 = get_rinf_uv0([reg1,reg2],v0=[u1,u2])
	## edit uvbounds
	for b in reg1.blocks:
		b.uvbounds.update(dict(vmax=u1))
	for b in reg2.blocks:
		b.uvbounds.update(dict(vmin=u2))
	## passive slice of reg1
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), v0=u1, r0=r0)
	## warn if bad u0 or v0 value
	pslice, reg1 = xh.junc.slicecheck(pslice, reg1)
	## set U0(r) and V0(r) for target coords
	U0 = lambda r: pslice.U_of_r_at_v0(r)
	V0 = lambda r: pslice.V_of_r_at_u0(r)
	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[2], vblocks=[0,1,2], v0=u2, r0=r0, U0=U0, V0=V0)
	## warn if bad u0 or v0 value
	aslice, reg2 = xh.junc.slicecheck(aslice, reg2)
	## update coordinate transformations
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	## return
	return [reg1, reg2]








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
				

def boundarylines(reglist, npoints=5001, sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			style = dict(lw=0.9, c='0.5', zorder=2000)
			style.update(sty)
			b.add_curves_uv(xh.cm.block_boundary(b, sty=style, npoints=npoints))


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
	reglist = [xh.reg.EFreg(xh.mf.minkowski(),rlines=False, boundary=False)]
	## accrete
	for i in range(len(v)):
		vx, Rx = 1.*v[i], 2.*m[i]
		reglist += xh.evap.accrete(reglist.pop(), v1=vx, v2=vx, R2=Rx, L2=l, rlines=False, boundary=False, rparams2=rparams)
	## return
	return reglist















"""
Tests.
Run if __name__='__main__'.
"""

def test1():
	"""
	Test functionality of get_rinf_uv0.
	"""
	##
	print "\nTEST 1\n"
	## regions
	reg1 = xh.reg.EFreg(xh.mf.schwarzschild())
	reg2 = xh.reg.EFreg(xh.mf.hayward())
	reg3 = xh.reg.EFreg(xh.mf.minkowski())
	reglist = [reg1,reg2,reg3]
	## v values
	v0 = [-1.,1.,0.]
	## get
	rinf = get_rinf_uv0(reglist,v0=v0)
	## print
	print rinf
	##
	print "\nEND TEST 1\n"



def test2():
	"""
	Test functionality of accrete.
	"""
	##
	print "\nTEST 2\n"
	## create initial region
	reglist = [xh.reg.EFreg(xh.mf.minkowski(),rlines=False,boundary=False)]
	## create accreted regions
	reglist += xh.evap.accrete(reglist.pop(), v1=0., v2=0., R2=0.8)
	reglist += xh.evap.accrete(reglist.pop(), v1=.5, v2=.5, R2=.9)
	reglist += xh.evap.accrete(reglist.pop(), v1=1., v2=1., R2=1.)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	plt.figure()
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	plt.show()
	##
	print "\nEND TEST 2\n"



def test3():
	"""
	Test functionality of accretion.
	"""
	##
	print "\nTEST 3\n"
	## params
	v = np.linspace(0,1,5)
	m = .2 + .8*v
	## create initial region
	reglist = accretion(m,v)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	plt.figure()
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	plt.show()
	##
	print "\nEND TEST 3\n"



def test4():
	"""
	Test functionality of accrete.
	"""
	##
	print "\nTEST 4\n"
	## create initial region
	reglist = [xh.reg.EFreg(xh.mf.hayward(R=1.),rlines=False,boundary=False)]
	## create evaporated regions
	reglist += xh.evap.evap(reglist.pop(), u1=0., u2=0., R2=0.8)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	plt.figure()
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	plt.show()
	##
	print "\nEND TEST 4\n"


if __name__=='__main__':
	#test1()
	#test2()
	#test3()
	test4()


