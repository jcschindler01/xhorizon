
import numpy as np
import matplotlib.pyplot as plt

import xhorizon as xh



def go(st='s'):
	## params
	R1 = 1.
	R2 = .8
	l = .1
	rparams = dict(s0=10.)
	u1 = 22.
	v1 = 0.
	u2 = 1.*u1
	## funcs
	if False:
		func1 = xh.mf.hayward(R=R1,l=l)
		func2 = xh.mf.hayward(R=R2,l=l)
	else:
		func1 = xh.mf.schwarzschild(R=R1)
		func2 = xh.mf.schwarzschild(R=R2)
	## regs
	reg1 = xh.reg.EFreg(func1,rparams=rparams,rlines=False,boundary=False)
	reg2 = xh.reg.EFreg(func2,rparams=rparams,rlines=False,boundary=False)
	## horizon radii
	Rh1, Rh2 = reg1.blocks[-1].rj[0], reg2.blocks[-1].rj[0]
	## passive slice of reg1
	r1 = reg1.blocks[-1].r_of_uv(np.array([[u1],[v1]]))[0]
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=r1, u0=u1)

	## active slice of reg2
	aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=pslice.r0, u0=u2, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0, r_refs=[pslice.reg.metfunc.r_ref])
	
	## plot slice
	if True:
		print 'REG1 PSLICE U(r,v0)=red, V(r,u0)=blue'
		print 'REG2 ASLICE U(r,v0)=red, V(r,u0)=blue'
		rmin, rmax, npoints = 0., 2., 10001
		rr = np.linspace(rmin,rmax,npoints)
		plt.figure()
		plt.title('REG1 PSLICE U(r,v0)=red, V(r,u0)=blue\nREG2 ASLICE U(r,v0)=mag, V(r,u0)=cyan')
		plt.xlabel('r')
		plt.ylabel('U(r),V(r) on slice')
		plt.grid()
		plt.xlim(.5,1.5)
		## passive slice
		style1 = dict(marker='x', markersize=8, lw=1, ls='none')
		plt.plot(pslice.r, pslice.UV_v0[0], c='r', **style1)
		plt.plot(pslice.r, pslice.UV_u0[1], c='b', **style1)
		## active slice
		style2 = dict(marker='o', markersize=6, lw=1, ls='none')
		plt.plot(aslice.r, aslice.UV_v0[0], c='m', **style2)
		plt.plot(aslice.r, aslice.UV_u0[1], c='c', **style2)
		## junction radius
		plt.plot([pslice.r0, pslice.r0], [-1.,1.], 'k-' )
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
	reglist = []

	## split regions
	reglist1 = xh.evap.split_reg_abcd(reg1,  abcd='bcd', u0=pslice.u0, v0=pslice.v0)
	reglist2 = xh.evap.split_reg_abcd(reg2,  abcd='a',   u0=aslice.u0, v0=aslice.v0)

	## update transformations for reg2
	for i in range(len(reglist2)):
		reglist2[i].U_of_udl = aslice.U_of_udl_at_v0
		reglist2[i].V_of_vdl = aslice.V_of_vdl_at_u0

	## add to reglist
	reglist += reglist1
	reglist += reglist2

	## print block uvbounds
	print "\n"
	print '%22r, %22r, %22r, %22r, %22r, %22r'%('R', 'b.rj', 'vmin', 'vmax', 'umin', 'umax')
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in reg.blocks:
				print '%22s, %22s, %22r, %22r, %22r, %22r'%(reg.metfunc.fparams['R'], b.rj, b.uvbounds['vmin'], b.uvbounds['vmax'], b.uvbounds['umin'], b.uvbounds['umax'])
	print "\n"

	## draw diagram
	if True:
		## add lines
		xh.evap.colorlines(reglist, sty=dict(ls='-'))
		xh.evap.boundarylines(reglist)
		xh.evap.s0_lines(reglist, sty=dict(lw=3))
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
















