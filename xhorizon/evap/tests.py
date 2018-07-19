


"""
Tests for evap module.
"""

import numpy as np
import matplotlib.pyplot as plt
import xhorizon as xh

from evap import *
from helpers import *


def main():
	#test1()
	#test2()
	#test3()
	test4()
	test5()


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
	Test functionality of evap and evap_final.
	"""
	##
	print "\nTEST 4\n"
	## create initial region
	reg0 = xh.reg.EFreg(xh.mf.hayward(R=1.,l=0.1),rlines=False,boundary=False)
	# for b in reg0.blocks:
	# 	b.uvbounds.update(dict(vmin=-5.))
	# for b in [reg0.blocks[2]]:
	# 	b.uvbounds.update(dict(umin=-5.))
	reglist = [reg0]
	## create evaporated regions
	reglist += xh.evap.evap(reglist.pop(), u1=3., u2=3., R2=0.9)
	reglist += xh.evap.evap_final(reglist.pop(), u1=6., u2=6.)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## plot list
	#reglist = reglist[0:3]
	## draw
	plt.figure()
	plt.title("Test 4")
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	## fill
	fill_by_R(reglist)
	## show plot
	#plt.show()
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
	## params
	m = 0.5 * np.array([.9,0.])
	u = np.array([.3,.6])
	## create evaporated regions
	reglist += xh.evap.evaporation(m, u, reglist.pop(), l=0.1, rparams={})
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	plt.figure()
	plt.title('Test 5')
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	## fill
	fill_by_R(reglist)
	## show plot
	plt.show()
	##
	print "\nEND TEST 5\n"


main()
