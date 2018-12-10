"""
This module contains helpers to draw with evap code.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

import xhorizon as xh
from xhorizon.shell_junction import interpolators as interp


############## plot reglist ####################

def drawreg(reglist, chainparams):
	"""
	Plot all regions in reglist.
	"""
	## styles
	rline_zero_sty      = dict(lw=0.9, c='0.5', zorder=11000)
	singularity_sty     = dict(ls='dashed', dashes=(2,2))
	rline_inf_sty       = dict(c='0.5', zorder=10000, lw=1.3)
	rline_hor_sty       = dict(lw=0.2, c='k', zorder=9999)
	acc_shells_sty      = dict(lw=0.2, ls='dashed', dashes=(4,4), c='0.75', zorder=2000)
	evap_shells_out_sty = acc_shells_sty
	evap_shells_in_sty  = dict(lw=0.4, ls='dashed', dashes=(.5,1), c='0.75', zorder=2000)
	fill_horizons_sty   = dict(fc='none', ec='k', lw=0, hatch='.', zorder=9990)
	fill_density_sty    = dict(zorder=100)

	## lines
	rline_zero(reglist, sty=rline_zero_sty, sty2=singularity_sty, npoints=5001)
	rline_inf_col(reglist, sty=rline_inf_sty, npoints=5001)
	rline_hor(reglist, sty=rline_hor_sty)
	acc_shells(reglist, chainparams, sty=acc_shells_sty, inf=100., npoints=5001)
	evap_shells_out(reglist, chainparams, sty=evap_shells_out_sty, inf=100., npoints=5001)
	evap_shells_in(reglist, chainparams, sty=evap_shells_in_sty, inf=100., npoints=5001)
	make_rlines(reglist, chainparams, l=.05)
	vticks(reglist)
	uticks(reglist)
	## plot
	for reg in reglist:
		reg.rplot()
	## fill
	fill_horizons(reglist, sty=fill_horizons_sty)
	fill_density(reglist, sty=fill_density_sty)




#############################################################


## being actually used ############

def rline_zero(reglist, npoints=5001, sty={}, sty2={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			if b.j==0:
				## style
				style = dict(lw=0.9, c='0.5', zorder=11000)
				style.update(sty)
				## singularity style
				if b.sgnf<0.:
					style.update(sty2)
				## curve
				cv = xh.cm.rstarlines_special_2([0.], b.uvbounds, c=0., sty=style, inf=100., npoints=1001, eps=1e-12)
				b.add_curves_uv(cv)



def rline_inf_col(reglist, npoints=5001, inf=1000., sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			## outermost blocks only
			if b.j==len(b.master.metfunc.rj)-2:
				## style
				style = dict(c='0.5', zorder=10000)
				style.update(sty)
				## rstar value
				rstar = b.master.metfunc.rstar_ref[-1]
				## mass
				m = 0.
				if 'R' in b.master.metfunc.fparams.keys():
					m = 0.5 * b.master.metfunc.fparams['R']
				## linewidth
				lw0 = 1.*style['lw']
				lw = lw0 * (1.+2.*m)
				style.update(dict(lw=lw))
				## curve
				b.add_curves_uv(xh.cm.block_boundary_2(b, sty=style, inf=100., npoints=npoints)[-1:])
				#b.add_curves_uv(xh.cm.uvlines([-inf,inf], uv='uv', uvbounds=b.uvbounds, sty=style, c=0., inf=inf, npoints=npoints))


def rline_hor(reglist, npoints=5001, sty={}):
	"""
	Add horizon lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			## trapped blocks only
			if b.sgnf<0.:
				## style
				style = dict(lw=0.2, c='k', zorder=9999)
				style.update(sty)
				## rstar value
				rstar = b.master.metfunc.rstar_ref[-1]
				## curve
				cvs = xh.cm.block_boundary_2(b, sty=style, inf=100., npoints=npoints)
				if b.j==0:
					cvs = cvs[-1:]
				## add
				b.add_curves_uv(cvs)


def acc_shells(reglist, chainparams, sty={}, inf=50., npoints=5001):
	"""
	Add boundary lines to regions in reglist.
	"""
	Rh = chainparams['Rh']
	v0 = chainparams['fs_v0']-1e-15
	for i in range(len(reglist[:-1])):
		reg = reglist[i]
		## accretion only
		if Rh[i+1]>Rh[i]:
			##
			for b in reg.blocks:
				## style
				style = dict(lw=0.2, ls='dashed', dashes=(4,4), c='0.65', zorder=2000)
				style.update(sty)
				## v value
				vv = v0[i:i+1]
				## curve
				b.add_curves_uv(xh.cm.uvlines(vv, uv='v', uvbounds=b.uvbounds, sty=style, c=0., inf=inf, npoints=npoints))

def evap_shells_out(reglist, chainparams, sty={}, inf=100., npoints=5001):
	"""
	Add boundary lines to regions in reglist.
	"""
	Rh = chainparams['Rh']
	u0 = chainparams['fs_u0']-1e-15
	for i in range(len(reglist[:-1])):
		reg = reglist[i]
		## evap only
		if Rh[i+1]<Rh[i]:
			## outgoing
			for b in reg.blocks:
				if not np.isfinite(b.uvbounds['vmax']):
					## style
					style = dict(lw=0.2, ls='dashed', dashes=(4,4), c='0.65', zorder=2000)
					style.update(sty)
					## u value
					uu = u0[i:i+1]
					## curve
					b.add_curves_uv(xh.cm.uvlines(uu, uv='u', uvbounds=b.uvbounds, sty=style, c=0., inf=inf, npoints=npoints))


def evap_shells_in(reglist, chainparams, sty={}, inf=5., npoints=5001):
	"""
	Add boundary lines to regions in reglist.
	"""
	Rh = chainparams['Rh']
	u0 = chainparams['fs_u0']
	v0 = chainparams['fs_v0']-1e-15
	r0 = chainparams['fs_r0']
	for i in range(len(reglist[:-1])):
		reg = reglist[i]
		## evap only
		if Rh[i+1]<Rh[i]:
			## uv value
			uu = u0[i:i+1]
			vv = v0[i:i+1]
			rr = r0[i:i+1]
			## ingoing
			for b in reg.blocks:
				## style
				style = dict(lw=0.9, ls=':', c='0.65', zorder=2000)
				style.update(sty)
				## inner blocks
				if b.rj[1]<rr:
					## values
					u = np.linspace(-inf,inf,npoints)
					v = vv+0.*u
					## make curve
					cv = xh.curve()
					cv.uv = np.array([u,v])
					cv.sty.update(style)
					## add
					b.add_curves_uv([cv])
				## outer block
				elif b.rj[0]<rr and not np.isfinite(b.uvbounds['umax']):
					## values
					u = np.linspace(uu,inf,npoints)
					v = vv+0.*u
					## make curve
					cv = xh.curve()
					cv.uv = np.array([u,v])
					cv.sty.update(style)
					## add
					b.add_curves_uv([cv])




def fill_horizons(reglist, sty={}):
	##
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in reg.blocks:
			if b.sgnf < 0.:
				style = dict(fc='none', ec='k', lw=0, hatch='.', zorder=9990)
				style.update(sty)
				b.fill(sty=style)


def fill_density(reglist, sty={}):
	## regions
	for reg in reglist:
		## linespacing
		x = np.arange(0.,10.1,.5)
		## default to zero
		density = 0.*x
		## vac spacetimes
		vac = ['Schwarzschild', 'Minkowski']
		if reg.metfunc.info['Type'] in vac:
			for b in reg.blocks:
				b.fill(dstyle(0.))
		## hayward
		if reg.metfunc.info['Type'] == 'Hayward':
			## get params
			l = reg.metfunc.fparams['l']
			R = reg.metfunc.fparams['R']
			## rcore and r
			rcore = (R*l**2.)**(1./3.)
			r = rcore * x
			rmax = 100.
			print r
			r = np.concatenate([r, np.array([rmax])])
			print r
			## density 
			density = 1./(1.+x**3.)**2.   # rho/rho0 = 1/(1+x^3)^2
			## blocks
			for b in reg.blocks:
				## loop over radii
				for n in range(len(x)):
					## rr and rho
					rr = r[n:n+2]
					rho = density[n]
					## style
					style = dstyle(rho)
					## fill
					rvals=b.fill_between_r(rr, sty=style, npoints=1001, inf=50.)



def dstyle(rho, cm=plt.cm.Oranges, sty={}, ovr={}):
	"""
	Density rho in (0,1), output fill style.
	Since all Hayward have same l, all regions same scale.
	"""
	## style
	style = dict(fc='none', ec='none', lw=0, zorder=900)
	style.update(sty)
	## cnorm
	cmin, cmax = .1, .5
	cnorm = cmin + (cmax-cmin)*rho
	## anorm
	amin, amax = .3, .7
	anorm = amin + (amax-amin)*rho
	## color and alpha
	style.update(fc=cm(cnorm), alpha=.47)
	## override
	style.update(ovr)
	## return
	return style.copy()
		
			

def make_rlines(reglist, chainparams, l=.05):
	"""
	"""
	## params
	R = np.max(chainparams['Rh'])
	l = 1.*l
	rcore = (R*l**2)**(1./3.)
	print "%.3f, %.3f, %.3f"%(l, rcore, R)
	## general
	style = dict(ls='-', zorder=10, lw=.4, alpha=.14)
	## l scale
	x = np.arange(0.,10.01,.5)
	scale = l
	style.update(dict(c='c'))
	rlines(reglist, scale*x, sty=style, inf=10., npoints=5001.)
	## R scale
	x = np.arange(0.,25.01,.5)
	scale = R
	style.update(dict(c='m' ))
	rlines(reglist, scale*x, sty=style, inf=25., npoints=5001.)


def vticks(reglist, dv=1., inf1=100., inf2=50.5):
	"""
	Remainder lets dv carry across regions.
	"""
	## init
	remainder = 0.
	## regions
	for reg in reglist:
		## blocks
		for b in reg.blocks:
			## last blocks only
			if b.j == reg.blocks[-1].j:
				## outer last blocks only
				if not np.isfinite(b.uvbounds['umin']):
					## get min and max
					vmin = np.max([b.uvbounds['vmin'], -inf2])
					vmax = np.min([b.uvbounds['vmax'],  inf2])
					## make array
					vv = remainder + np.arange(vmin, vmax, dv)
					## get new remainder
					remainder = dv - (vmax - vv[-1])
					## make and add curve
					crv = xh.curve()
					crv.uv = np.array([-inf1+0.*vv, 1.*vv])
					style = dict(c='.5',alpha=.3, ls='none', marker=(2,0,45), markersize=7, zorder=100)
					crv.sty.update(style)
					b.add_curves_uv([crv])
					## print
					print "dv = %s"%(dv)



def uticks(reglist, du=1., inf1=100., inf2=50.5):
	"""
	Remainder lets dv carry across regions.
	"""
	## init
	remainder = 0.
	## regions
	for reg in reglist:
		## blocks
		for b in reg.blocks:
			## last blocks only
			if b.j == reg.blocks[-1].j:
				## outer last blocks only
				if not np.isfinite(b.uvbounds['umax']):
					## get min and max
					umin = np.max([b.uvbounds['umin'], -inf2])
					umax = np.min([b.uvbounds['umax'],  inf2])
					## make array
					uu = remainder + np.arange(umin, umax, du)
					## get new remainder
					remainder = du - (umax - uu[-1])
					## make and add curve
					crv = xh.curve()
					crv.uv = np.array([1.*uu, inf1+0.*uu])
					style = dict(c='.5',alpha=.3, marker=(2,0,-45), markersize=7, zorder=100)
					crv.sty.update(style)
					b.add_curves_uv([crv])
					## print
					print "dv = %s"%(du)



##############################

################# tools for drawing ########################################



def rlines(reglist, rvals, sty={}, npoints=2001, inf=25.):
	"""
	Add colorscaled lines of constant radius to region.
	Useful to check for matching.
	"""
	rvals = 1.*rvals
	for reg in reglist:
		for b in reg.blocks:
			for i in range(len(rvals)):
				if (b.rj[0]<rvals[i]) and (rvals[i]<b.rj[1]):
					style = dict(c='c', ls='-', lw=.5, zorder=1500)
					style.update(sty)
					rstarvals = 1.*b.master.metfunc.F(rvals[i:i+1])
					b.add_curves_uv(xh.cm.rstarlines_special_2(rstarvals, b.uvbounds, c=b.master.rparams['c'], sty=style, inf=1.*inf, npoints=1.*npoints))
	return reglist


def colorlines(reglist, rmin=1.1, rmax=25., dr=.5, sty={}, cm=plt.cm.hsv_r, npoints=2001, inf=25.):
	"""
	Add colorscaled lines of constant radius to region.
	Useful to check for matching.
	"""
	rvals = np.arange(rmin,rmax,dr)
	colvals = (rvals - np.min(rvals)) / (np.max(rvals) - np.min(rvals))
	cols = cm(colvals)
	for reg in reglist:
		for b in reg.blocks:
			for i in range(len(rvals)):
				if (b.rj[0]<rvals[i]) and (rvals[i]<b.rj[1]):
					style = dict(c=cols[i], ls='-', markersize=1, lw=.5, zorder=1500)
					style.update(sty)
					rstarvals = 1.*b.master.metfunc.F(rvals[i:i+1])
					b.add_curves_uv(xh.cm.rstarlines_special_2(rstarvals, b.uvbounds, c=b.master.rparams['c'], sty=style, inf=2.*inf, npoints=npoints))
	return reglist

def s0_lines(reglist, sty={}, npoints=1001, inf=50., eps=1e-24):
	"""
	"""
	for reg in reglist:
		for b in reg.blocks:
			x = np.array([(2.*reg.rparams['s0'] - eps)])
			style1 = dict(c='m', alpha=.5, lw=3, ls=':', zorder=5000)
			style2 = dict(c='c', alpha=.5, lw=3, ls=':', zorder=5000)
			style1.update(sty)
			style2.update(sty)
			b.add_curves_uv(xh.cm.uvlines( x, uv='uv', uvbounds=b.uvbounds, sty=style1, c=0., inf=inf, npoints=npoints))
			b.add_curves_uv(xh.cm.uvlines(-x, uv='uv', uvbounds=b.uvbounds, sty=style2, c=0., inf=inf, npoints=npoints))
	return reglist

##########################################################################################





if __name__=='__main__':
	xh.evap.demo()

