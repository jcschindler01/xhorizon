


"""
Tests for evap module.
"""

import numpy as np
import matplotlib.pyplot as plt
import xhorizon as xh

from evap import *
from helpers import *
from accrete import *
from formevap import *


def main():
	#test1()
	#test2()
	#test3()
	#test4()
	test5()
	#test6()
	#test7()

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
	Test functionality of evap.
	"""
	##
	print "\nTEST 4\n"
	## create initial region
	reg0 = xh.reg.EFreg(xh.mf.hayward(R=1.,l=0.1),rlines=False,boundary=False)
	reglist = [reg0]
	## create evaporated regions
	reglist += xh.evap.evap(reglist.pop(), u1=0., v1=0., u2=0., R2=0.99)
	reglist += xh.evap.evap(reglist.pop(), u1=2., v1=2., u2=0., R2=0.98)
	## draw diagram?
	if True:
		## add lines
		xh.evap.colorlines(reglist)
		xh.evap.boundarylines(reglist)
		## draw
		xh.newfig(tex=False,sqaxis=3)
		plt.title("Test 4")
		for reg in reglist:
			reg.rplot()
		## fill
		fill_by_R(reglist)
		## show plot
		plt.savefig("temp-figs/test4.png", dpi=200)
		##
	print "\nEND TEST 4\n"


def test5():
	"""
	Test functionality of evaporation().
	"""
	##
	print "\nTEST 5\n"
	## create initial region
	reg0 = xh.reg.EFreg(xh.mf.hayward(R=1.,l=0.1),rlines=False,boundary=False)
	reglist = [reg0]
	## init params
	u0 = 5.
	v0 = -3.
	## subsequent params
	Nreg = 2
	tt = np.linspace(0.,1.,Nreg+2)[1:-1]
	R = 1. * (1.-tt)**(1./3.)
	du = 0.4 + 0.*tt
	dv = 0.2 + 0.*tt
	## create evaporated regions
	reglist += xh.evap.evaporation(reglist.pop(), R=R, du=du, dv=dv, u0=u0, v0=v0, l=0.1, rparams={})
	############
	## check ubound values
	uu = [reglist[0].blocks[-1].uvbounds['umin'], reglist[0].blocks[-1].uvbounds['umax']]
	for reg in reglist:
		a, b = reg.blocks[-1].uvbounds['umin'], reg.blocks[-1].uvbounds['umax']
		if np.isfinite(a) and np.isfinite(b):
			uu += [a,b]
	print "\n"
	print "uu = %r"%(uu)
	## draw regions?
	if True:
		## add lines
		xh.evap.colorlines(reglist)
		xh.evap.boundarylines(reglist)
		## draw
		xh.newfig(tex=False,sqaxis=3)
		plt.title('Test 5')
		for reg in reglist:
			reg.rplot()
		## fill
		fill_by_R(reglist, cm=plt.cm.prism)
		## show plot
		plt.savefig("temp-figs/test5.png", dpi=200)
		##
	print "\nEND TEST 5\n"


def test6():
	"""
	Test functionality of formevap().
	"""
	##
	print "\nTEST 6\n"
	## accrete params
	N = 4.
	R0, R = .8, 1.
	accrete_R = R0 + (R-R0)*np.linspace(0,1,N)
	v0, v = -1., -.5
	accrete_v = np.linspace(v0,v,N)
	## evap params
	evap_R = [.9 ,.0]
	evap_u = [6.,9.]
	## create evaporated regions
	reglist = xh.evap.formevap(accrete_R=accrete_R, accrete_v=accrete_v, evap_R=evap_R, evap_u=evap_u, l=0.01, rparams=dict(c=0., s0=10.))
	## resquish
	fU = None
	fV = None
	reglist = UVcompose(reglist, fU=fU, fV=fV)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	xh.newfig(tex=False,sqaxis=3)
	#plt.figure()
	plt.title('Test 6')
	plt.xlim(0.5,2.5)
	plt.ylim(-1,1)
	for reg in reglist:
		reg.rplot()
	## fill
	fill_by_R(reglist)
	## mink compare
	#xh.reg.EFreg(xh.mf.minkowski()).rplot()
	## show plot
	plt.savefig("temp-figs/test6.png", dpi=600)
	##
	print "\nEND TEST 6\n"


def test7():
	"""
	Test functionality of get_uvdl_of_UV().
	"""
	##
	print "\nTEST 7\n"
	## create evaporated regions
	reglist = xh.evap.formevap(rparams=dict(s0=10.))
	## resquish
	fU, fV = get_uvdl_of_UV(reglist[-1])
	## plot
	xx = np.linspace(-2,2,50001)
	plt.figure()
	plt.plot(xx, fU(xx), 'rx-')
	plt.plot(xx, fV(xx), 'bx-')
	plt.grid()
	plt.show()
	## edit regions
	UVcompose(reglist, fU=fU, fV=fV)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	xh.newfig(tex=False,sqaxis=3)
	plt.title('Test 7')
	for reg in reglist:
		reg.rplot()
	fill_by_R(reglist)
	## show plot
	plt.show()
	plt.savefig("temp-figs/test7.png", dpi=400)
	##
	print "\nEND TEST 7\n"





main()
