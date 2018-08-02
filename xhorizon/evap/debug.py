
import numpy as np
import matplotlib.pyplot as plt

import xhorizon as xh



def go():
	## params
	R1 = 1.
	R2 = .5
	l = .1
	rparams = dict(s0=10.)
	u1 = 8.
	v1 = 0.
	u2 = 1.*u1
	## region 1
	func1 = xh.mf.hayward(R=R1,l=l)
	reg1 = xh.reg.EFreg(func1,rparams=rparams,rlines=False,boundary=False)
	## region 2
	func2 = xh.mf.hayward(R=R2,l=l)
	reg2 = xh.reg.EFreg(func2,rparams=rparams,rlines=False,boundary=False)
	## horizon radii
	Rh1, Rh2 = reg1.blocks[-1].rj[0], reg2.blocks[-1].rj[0]
	## passive slice of reg1
	r1 = reg1.blocks[-1].r_of_uv(np.array([[u1],[v1]]))[0]
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=r1, u0=u1)
	## plot passive slice
	if False:
		rmin, rmax, npoints = 0., 2., 5001
		rr = np.linspace(rmin,rmax,npoints)
		plt.figure()
		plt.title('REG1 PSLICE U(r,v0)=red, V(r,u0)=blue')
		plt.xlabel('r')
		plt.ylabel('U(r),V(r) on slice')
		plt.grid()
		plt.plot(rr, pslice.U_of_r_at_v0(rr), 'r-', marker='x')
		plt.plot(rr, pslice.V_of_r_at_u0(rr), 'b-', marker='x')
	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=pslice.r0, u0=u2, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0)
	## print slice params
	print "\n"
	print "                          %22r, %22r, %22r, %22r, %22r"%('Rh', 'r', 'r/Rh', 'u', 'v')
	print "Pslice Inputs:"
	print "Rh1, r1, r1/Rh1, u1, v1 = %22r, %22r, %22r, %22r, %22r"%(Rh1, r1, r1/Rh1, u1, v1)
	print "From PSlice:"
	print "Rh1, r1, r1/Rh1, u1, v1 = %22r, %22r, %22r, %22r, %22r"%(Rh1, pslice.r0, pslice.r0/Rh1, pslice.u0, pslice.v0)
	print "Aslice Inputs:"
	print "Rh2, r2, r2/Rh2, u2, v2 = %22r, %22r, %22r, %22r, %22r"%(Rh2, pslice.r0, pslice.r0/Rh2, u2, '?')
	print "From ASlice:"
	print "Rh2, r2, r2/Rh2, u2, v2 = %22r, %22r, %22r, %22r, %22r"%(Rh2, aslice.r0, aslice.r0/Rh2, aslice.u0, aslice.v0)
	print "\n"
	## reglist
	reglist = [reg1,reg2]
	## draw diagram
	if False:
		## add lines
		xh.evap.colorlines(reglist)
		xh.evap.boundarylines(reglist)
		xh.evap.s0_lines(reglist, sty=dict(lw=.5))
		## draw
		plt.figure()
		plt.xlim(-3,3)
		plt.gca().set_aspect('equal')
		for reg in reglist:
			reg.rplot()
		## fill
		xh.evap.fill_by_R(reglist, cm=plt.cm.prism)


	## show plots
	plt.show()
	##

go()
















