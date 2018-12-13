
"""
This module provides method for making forming and evaporation BH diagrams.
This module imports the entire xhorizon package. 
It is meant for a higher level usage than the other subpackages, none of the
guts of xhorizon rely on this.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy, pprint
import scipy.optimize as opt

import xhorizon as xh
from helpers import *



###############################################################################################################3



def funclist_chain(funclist, seed=0, du=None, dv=None, r0p=None, r0f=None, u0=None, v0=None, ps_matchmode=None, fs_matchmode=None):
	"""
	Create a chain of matched regions, starting at seed region which is unmodified.
	Each region except ends has two slices through it, a future slice fslice and past slice pslice.
	Each fslice and pslice can be either active or passive, but there can only be one active slice per region.
	The index i refers to each region in the sequence for all variables.

	Inputs:
		funclist = list of func objects, in order, to chain together
		seed = index value for seed region (seed region has trivial transforms to target coords)
		du = list of du values so that du[i] will always be size of region[i]
		dv = list of du values so that du[i] will always be size of region[i]
		r0p = list of r0 values for past slice so that r0p will always be ps_r0 when pslice is active
		r0f = list of r0 values for future slice so that r0f will always be fs_r0 when fslice is active
		u0 = list of offset values for range of u values in slice, defaults to zero
		v0 = list of offset values for range of v values in slice, defaults to zero
		ps_matchmode = list of strings, each either 'ru' or 'rv', to determine how past slice is sliced when pslice is active
		ps_matchmode = list of strings, each either 'ru' or 'rv', to determine how future slice is sliced when fslice is active
	"""
	## init default values
	if u0==None:
		u0 = np.zeros(len(funclist))
	if v0==None:
		v0 = np.zeros(len(funclist))
	if ps_matchmode==None:
		ps_matchmode = ['rv' for func in funclist]
	if fs_matchmode==None:
		fs_matchmode = ['rv' for func in funclist]
	## set irrelevant first and last du and dv values to zero
	du[0], du[-1] = 0., 0.
	dv[0], dv[-1] = 0., 0.	
	## init internal variables
	reglist = [xh.reg.EFreg(funcx, boundary=False, rlines=False) for funcx in funclist]
	pslice  = [None for funcx in funclist]
	fslice  = [None for funcx in funclist]
	Rh      = [funcx.rj[-2] for funcx in funclist]
	ps_r0   = [np.nan for funcx in funclist]
	ps_u0   = [np.nan for funcx in funclist]
	ps_v0   = [np.nan for funcx in funclist]
	fs_r0   = [np.nan for funcx in funclist]
	fs_u0   = [np.nan for funcx in funclist]
	fs_v0   = [np.nan for funcx in funclist]
	i0 = range(len(funclist))[1*seed]
	ps_matchpop = [mp(mmm) for mmm in ps_matchmode]
	fs_matchpop = [mp(mmm) for mmm in fs_matchmode]
	## seed region
	i = 1*i0
	for i in [1*i0]:
		###### past passive slice
		## past passive slice input params (not mutually consistent)
		ps_u0[i]  = u0[i] - 0.5*du[i]
		ps_v0[i]  = v0[i] - 0.5*dv[i]
		ps_r0[i]  = 1.*r0p[i]
		## get past passive slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(ps_matchpop[i])
		print "i=%s pslice loc: %s"%(i,sliceloc)
		## execute past passive slice at sliceloc
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update past passive slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		#### future passive slice
		## future passive slice input params (not mutually consistent)
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*r0f[i]
		## get future passive slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(fs_matchpop[i])
		print "i=%s fslice loc: %s"%(i,sliceloc)
		## execute future passive slice at sliceloc
		fslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		fs_u0[i], fs_v0[i], fs_r0[i]  = 1.*fslice[i].u0, 1.*fslice[i].v0, 1.*fslice[i].r0
	## forward regions
	i = 1*i0 + 1
	while i < len(reglist):
		###### past active slice
		## past active slice input params (not mutually consistent)
		ps_u0[i]  = u0[i] - 0.5*du[i]
		ps_v0[i]  = v0[i] - 0.5*dv[i]
		ps_r0[i]  = 1.*fs_r0[i-1]
		## get past active slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(ps_matchpop[i])
		print "i=%s pslice loc: %s"%(i,sliceloc)
		## execute past active slice at sliceloc
		pslice[i] = xh.junc.aslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), U0=fslice[i-1].U_of_r_at_v0, V0=fslice[i-1].V_of_r_at_u0, r_refs=[fslice[i-1].reg.metfunc.r_ref], **sliceloc)
		## update past active slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		#### modify transformations
		## adjust transformations
		reglist[i].U_of_udl = pslice[i].U_of_udl_at_v0
		reglist[i].V_of_vdl = pslice[i].V_of_vdl_at_u0
		#### future passive slice
		## future passive slice input params (not mutually consistent)
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*r0f[i]
		## get past active slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(fs_matchpop[i])
		print "i=%s fslice loc: %s"%(i,sliceloc)
		## execute future passive slice at sliceloc
		fslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		fs_u0[i], fs_v0[i], fs_r0[i]  = 1.*fslice[i].u0, 1.*fslice[i].v0, 1.*fslice[i].r0
		##### iterate
		## iterate
		i += 1
	## backward regions
	i = 1*i0 - 1
	while i>=0:
		###### future active slice
		## past active slice input params (not mutually consistent)
		fs_u0[i]  = u0[i] - 0.5*du[i]
		fs_v0[i]  = v0[i] - 0.5*dv[i]
		fs_r0[i]  = 1.*ps_r0[i+1]
		## get future active slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(fs_matchpop[i])
		print "i=%s fslice loc: %s"%(i,sliceloc)
		## execute future active slice at sliceloc
		fslice[i] = xh.junc.aslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), U0=pslice[i+1].U_of_r_at_v0, V0=pslice[i+1].V_of_r_at_u0, r_refs=[pslice[i+1].reg.metfunc.r_ref], **sliceloc)
		## update future active slice location to true values
		fs_u0[i], fs_v0[i], fs_r0[i]  = 1.*fslice[i].u0, 1.*fslice[i].v0, 1.*fslice[i].r0
		#### modify transformations
		## adjust transformations
		reglist[i].U_of_udl = fslice[i].U_of_udl_at_v0
		reglist[i].V_of_vdl = fslice[i].V_of_vdl_at_u0
		#### past passive slice
		## past passive slice input params (not mutually consistent)
		ps_u0[i]  = 1.*fs_u0[i] - 1.*du[i]
		ps_v0[i]  = 1.*fs_v0[i] - 1.*dv[i]
		ps_r0[i]  = 1.*r0p[i]
		## get past passive slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(ps_matchpop[i])
		print "i=%s pslice loc: %s"%(i,sliceloc)
		## execute past passive slice at sliceloc
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		##### iterate
		## iterate
		i -= 1
	## make sliceparams dict
	chainparams = dict(Rh=1.*np.array(Rh), ps_u0=1.*np.array(ps_u0), ps_v0=1.*np.array(ps_v0), ps_r0=1.*np.array(ps_r0), fs_u0=1.*np.array(fs_u0), fs_v0=1.*np.array(fs_v0), fs_r0=1.*np.array(fs_r0), i0=1*i0, ps_matchmode=ps_matchmode, fs_matchmode=fs_matchmode, funclist=funclist)
	##
	print "\n"
	pprint.pprint(chainparams)
	print "\n"
	## return
	return reglist, chainparams






def chain_masker(reglist, chainparams):
	"""
	"""
	##
	for i in range(len(reglist)):
		## mask interior blocks
		for b in reglist[i].blocks[:-1]:
			## past
			if i>0:
				b.uvbounds.update(dict(vmin=chainparams['ps_v0'][i]))
			## future
			if i<len(reglist)-1:
				b.uvbounds.update(dict(vmax=chainparams['fs_v0'][i]))
		## mask final blocks for part that is always there
		for b in reglist[i].blocks[-1:]:
			## past
			if i>0:
				b.uvbounds.update(dict(vmin=chainparams['ps_v0'][i], umin=chainparams['ps_u0'][i]))
			## future
			if i<len(reglist)-1:
				b.uvbounds.update(dict(vmax=chainparams['fs_v0'][i]))
		## copy final block for parts which depend on radius change values
		for b in reglist[i].blocks[-1:]:
			## copies
			ba = xh.block(b.master, b.j, b.bparams)
			bb = xh.block(b.master, b.j, b.bparams)
			bc = xh.block(b.master, b.j, b.bparams)
			## mask a=top b=bottom c=right
			ba.uvbounds.update(dict(vmin=chainparams['fs_v0'][i], vmax= np.inf, umin=chainparams['ps_u0'][i], umax=chainparams['fs_u0'][i]))
			bb.uvbounds.update(dict(vmin=chainparams['ps_v0'][i], vmax=chainparams['fs_v0'][i], umin=-np.inf, umax=chainparams['ps_u0'][i]))
			bc.uvbounds.update(dict(vmin=chainparams['fs_v0'][i], vmax= np.inf, umin=-np.inf, umax=chainparams['ps_u0'][i]))
		## add bottom if increasing from past
		if i>0 and chainparams['Rh'][i-1]<chainparams['Rh'][i]:
			reglist[i].blocks += [bb]
		## add top if decreasing to future
		if i<len(reglist)-1 and chainparams['Rh'][i+1]<chainparams['Rh'][i]:
			reglist[i].blocks += [ba]
		## add right if both
		if i>0 and i<len(reglist)-1 and chainparams['Rh'][i-1]<chainparams['Rh'][i] and chainparams['Rh'][i+1]<chainparams['Rh'][i]:
			reglist[i].blocks += [bc]
	## add masses to chainparams
	chainparams.update(dict(m=getmm(reglist)))
	## return
	return reglist, chainparams





def dR_function(dR, func, dv=1., l=.01, A=10., Rmax=1.):
	"""
	This function is zero for valid values of dR, given the parameters.
	A is delta_u lifetime of BH.
	"""
	## params from func
	R0, Rh0, F = 1.*func.fparams['R'], 1.*func.rj[-2], func.F
	## du1 from du ~ R^2 dR
	du1 = 3.*A*(R0**2)*dR / (Rmax**3)
	## du2 from dv-du=2*drstar
	du2 = dv - 2. * ( F(Rh0+dR+l) - F(Rh0+l) )
	## du1=du2=du so du1-du2 = zero
	z = np.abs(du1 - du2)**2
	## return
	return 1.*z


def shellparams_from_func(func, dv=1., l=.01, A=10., Rmax=1.):
	"""
	Suppose you are given a region with metfunc func.
	This region has a known mass R0 and horizon radius Rh0, both given by func.
	Assume that this region has a future slice at r0 = Rh+l.
	Assume that this region will be attached to another region with mass R0+dR.
	Assume that the time difference between the two slices has a fixed value dv.
	Assume that the past slice will be at a radius Rh0+dR+l.
	Assume that this region has a time difference du = 3*A*(R0**2)*dR.
	Assume that Rh ~= R.
	This routine determines what value of dR and du makes the assumptions possible.
	"""
	## define dR_function
	dR_f = lambda dR: dR_function(dR, func, dv=1.*dv, l=1.*l, A=1.*A, Rmax=1.*Rmax)
	## plot dR_f
	if False:
		ds = np.linspace(-2,2,8001)
		plt.plot(ds, dR_f(ds), 'k-')
		plt.plot([0,2],[0,0],'k--')
		plt.grid()
		plt.ylim(-10,10)
		plt.show()
	## find root of dR_f
	dR = opt.minimize_scalar(dR_f, bounds=(0., 3.), method='bounded').x
	## find du in terms of R and dR
	R0 = 1.*func.fparams['R']
	du = 3.*A*(R0**2)*dR / (Rmax**3)
	## define shellparams
	shellparams = dict(func=copy.deepcopy(func), Rself=1.*func.fparams['R'], dR=1.*dR, du=1.*du, dv=1.*dv, l=1.*l, A=1.*A)
	shellparams.update(dict(Rnext=1.*shellparams['Rself']+1.*shellparams['dR']))
	## print
	pprint.pprint(shellparams)
	## return
	return shellparams.copy()



def shellparams_list(Rmin=.1, Rmax=1., dv=1., l=.1, A=10., functype=xh.mf.schwarzschild, fparams=dict()):
	"""
	"""
	## init
	R = 1.*Rmin
	## get shellparams
	sp = []
	while R<=Rmax:
		print R
		func = functype(R=1.*R, **fparams)
		sp += [shellparams_from_func(func, dv=1.*dv, l=1.*l, A=1.*A, Rmax=1.*Rmax)]
		R = 1.*sp[-1]['Rnext']
	## return
	return sp


def sp_transpose(sp_list):
	"""
	"""
	##
	sp = sp_list[::-1]
	## outs
	funclist = [sp[i]['func'] for i in range(len(sp))]
	RR = np.array([sp[i]['Rself'] for i in range(len(sp))])
	du = np.array([sp[i]['du'] for i in range(len(sp))])
	dv = np.array([sp[i]['dv'] for i in range(len(sp))])
	dR = np.array([sp[i]['dR'] for i in range(len(sp))])
	A = np.array([sp[i]['A'] for i in range(len(sp))])
	l = np.array([sp[i]['l'] for i in range(len(sp))])
	## make
	uu = np.cumsum(du)
	vv = np.cumsum(dv)
	## dict
	spt = dict(funclist=copy.deepcopy(funclist), R=1.*RR, dR=1.*dR, du=1.*du, dv=1.*dv, A=1.*A, l=1.*l, uu=1.*uu, vv=1.*vv)
	## return
	return spt.copy()


def cp_from_fdudv(funclist, du=None, dv=None, l=None, uoff=0., voff=0., ueta=1., veta=1.):
	"""
	"""
	## init
	funclist = funclist
	reglist = [xh.reg.EFreg(funcx, boundary=None, rlines=None) for funcx in funclist]
	Rh = np.array([funclist[i].rj[-2] for i in range(len(funclist))])
	du  = 1.*du
	dv  = 1.*dv
	r0f = 1.*Rh + 1.*l
	r0p = 1.*np.roll(r0f,1)
	u0  = 1.*ueta*np.cumsum(du-du[0]) + 1.*uoff
	v0  = 1.*veta*np.cumsum(dv-dv[0]) + 1.*voff
	ps_matchmode = ['rv' for i in range(len(funclist))]
	fs_matchmode = ['rv' for i in range(len(funclist))]
	## iterator
	ii = range(len(funclist))
	## get rinf
	rinf = np.nan * Rh
	for i in ii:
		ia, ib = max(0, i-1), min(i+2, len(ii))
		rinf[i] = get_rinf_uv0(reglist[ia:ib], v0=1.*v0)
	print rinf	
	## correct first and last r0 values
	r0p[0]  = 1.*rinf[0]
	r0f[-1] = 1.*rinf[-1]
	## correct r0 values for formation and evaporation
	for i in ii:
		## past
		if i>0:
			## accretion
			if Rh[i]>=Rh[i-1]:
				r0p[i] = 1.*rinf[i]
			## evaporation
			if Rh[i]< Rh[i-1]:
				r0p[i] = 1.*Rh[i-1] + 1.*l
		## future
		if i<len(ii)-1:
			## accretion
			if Rh[i]<=Rh[i+1]:
				r0f[i] = 1.*rinf[i]
			## evaporation
			if Rh[i]> Rh[i+1]:
				r0f[i] = 1.*Rh[i] + 1.*l
	## make cp
	cp = dict(du=1.*du, dv=1.*dv, r0p=1.*r0p, r0f=1.*r0f, u0=1.*u0, v0=1.*v0, ps_matchmode=ps_matchmode, fs_matchmode=fs_matchmode)
	# ## return
	return cp.copy()




def formevap_input(Rmin=.1, Rmax=1., dv_evap=1., l=.01, A=1., B=1., Naccrete=5, uoff=0., voff=0., ueta=1., veta=1., functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.schwarzschild, fparams1=dict()):
	"""
	Build inputs in reverse order starting from far future.

	funclist, seed=0, du=None, dv=None, r0p=None, r0f=None, u0=None, v0=None, ps_matchmode=None, fs_matchmode=None
	"""
	## init
	funclist = []
	du  = []
	dv  = []
	## final region
	funclist += [functype0(**fparams0)]
	du  += [0.]
	dv  += [0.]
	## evap
	sp = shellparams_list(Rmin=1.*Rmin, Rmax=1.*Rmax, dv=1.*dv_evap, l=1.*l, A=1.*A, functype=functype1, fparams=fparams1)
	for i in range(len(sp)):
		funclist += [sp[i]['func']]
		du  += [sp[i]['du']]
		dv  += [sp[i]['dv']]
	## max radius
	Rmax = sp[-1]['Rself']
	## accrete params
	RR = np.linspace(Rmax,0.5*Rmax, Naccrete)[1:]
	for R in RR:
		funclist += [functype1(R=1.*R, **fparams1)]
		du += [0.]
		dv += [B/float(Naccrete-1)]
	## first region
	funclist += [functype0(**fparams0)]
	du  += [0.]
	dv  += [0.]
	## prep for output
	funclist = funclist[::-1]
	du = np.array(du[::-1])
	dv = np.array(dv[::-1])
	l = 1.*l
	## get chain params
	cp = cp_from_fdudv(funclist, du=1.*du, dv=1.*dv, l=1.*l, uoff=1.*uoff, voff=1.*voff, ueta=1.*ueta, veta=1.*veta)
	##
	pprint.pprint(cp)
	## return
	return funclist, cp



def create_evap(params, seed=0):
	"""
	Takes input parameters of the form:

	"""
	## 
	import pprint
	## print
	pprint.pprint("params = %s"%(params))
	pprint.pprint("seed = %s"%(seed))
	## formevap_input
	print "inputs"
	funclist, cp = xh.evap.formevap_input(**params)
	## funclist_chain
	print "chain"
	reglist, chainparams = xh.evap.funclist_chain(funclist, seed=seed, **cp)
	## chain_masker
	print "mask"
	reglist, chainparams = xh.evap.chain_masker(reglist, chainparams)
	## print
	pprint.pprint(chainparams)
	## return
	return reglist, chainparams


def evapsave(path="temp/temp", params=None, chainparams=None, seed=None, sfp=dict(), temp_only=False, massplot=False):
	"""
	Save figure with timestamp and txt notes.
	"""
	##
	import shutil
	import time
	import pprint
	import matplotlib.pyplot as plt
	## get path with timestamp
	ts = str(time.time()).replace(".","")
	## save figure
	print( "save...")
	plt.figure(1)
	sfpp = dict(dpi=400)
	sfpp.update(sfp)
	plt.savefig("%s_%s.png"%(path,ts), **sfpp)
	print( "save done")
	##save text
	print( "save txt")
	ff = open("%s_%s.txt"%(path,ts), 'w')
	ff.write("%s_%s\n"%(path,ts))
	ff.write('\n')
	ff.write('Input:\nparams=\n%s\nseed=\n%s\n'%(pprint.pformat(params),seed))
	ff.write('\n')
	ff.write('Output:\nchainparams=\n%s\n'%(pprint.pformat(chainparams)))
	ff.close()
	##save massplot
	if massplot==True:
		print( "save massplot...")
		xh.evap.massplot.massplotrc()
		plt.figure(99)
		plt.savefig("%s_%s_mass.png"%(path,ts), **sfpp)
		print( "save done")
	## copy to temp
	print( "copy...")
	## copy normally
	if temp_only==False:
		shutil.copy("%s_%s.png"%(path,ts), path+"_temp.png")
		shutil.copy("%s_%s.txt"%(path,ts), path+"_temp.txt")
		shutil.copy("%s_%s_mass.png"%(path,ts), path+"_temp_mass.png")
	if temp_only==True:
		shutil.move("%s_%s.png"%(path,ts), path+"_temp.png")
		shutil.move("%s_%s_mass.txt"%(path,ts), path+"_temp.txt")
		shutil.move("%s_%s_mass.png"%(path,ts), path+"_temp_mass.png")
	print( "copy done")




if __name__=='__main__':
	pass








##################################################################################################################



