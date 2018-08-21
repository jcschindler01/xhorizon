
"""
This module provides method for making forming and evaporation BH diagrams.
This module imports the entire xhorizon package. 
It is meant for a higher level usage than the other subpackages, none of the
guts of xhorizon rely on this.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy, pprint

import xhorizon as xh
from helpers import *




###############################################################################################################3



def funclist_chain(funclist, seed=0, u0=None, v0=None, du=None, dv=None, eta=0., matchmode='ru'):
	"""
	Create a chain of matched regions, starting at seed region which is unmodified.
	Each region except ends has two slices through it, a future slice fslice and past slice pslice.
	Each fslice and pslice can be either active or passive, but there can only be one active slice per region.
	The index i refers to each region in the sequence for all variables.
	
	Slive u0 and v0 values determined as follows:
	If eta=0 then every region goes from u0-du/2 to u0+du/2, equally spaced about u0 with total size du.
	If eta=1  then the seed region is exactly as it would otherwise be, but the other regions are matched to it
	sequentially so that each begins where the last left off.
	This either applies to u or v, depending on matchmode setting.

	"""
	## init
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
	matchpop = mp(matchmode)
	## default u0 and v0 values
	if u0==None:
		[0. for funcx in funclist]
	if v0==None:
		[0. for funcx in funclist]
	## seed region
	i = 1*i0
	for i in [1*i0]:
		###### past passive slice
		## past passive slice input params (mutually consistent)
		ps_u0[i]  = u0[i] - 0.5*du[i]
		ps_v0[i]  = v0[i] - 0.5*dv[i]
		ps_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[ps_u0[i]],[ps_v0[i]]]))[0]
		## get past passive slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(matchpop)
		print "i=%s pslice loc: %s"%(i,sliceloc)
		## execute past passive slice at sliceloc
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update past passive slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		#### future passive slice
		## future passive slice input params (mutually consistent)
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[fs_u0[i]],[fs_v0[i]]]))[0]
		## get future passive slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(matchpop)
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
		ps_u0[i]  = (eta) * ( fs_u0[i-1] ) + (1.-eta) * ( u0[i] - 0.5*du[i] )
		ps_v0[i]  = (eta) * ( fs_v0[i-1] ) + (1.-eta) * ( v0[i] - 0.5*dv[i] )
		ps_r0[i]  = 1.*fs_r0[i-1]
		## get past active slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(matchpop)
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
		## future passive slice input params (mutually consistent)
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[fs_u0[i]],[fs_v0[i]]]))[0]
		## get past active slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(matchpop)
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
		fs_u0[i]  = (eta) * ( ps_u0[i+1] ) + (1.-eta) * ( u0[i] - 0.5*du[i] )
		fs_v0[i]  = (eta) * ( ps_v0[i+1] ) + (1.-eta) * ( v0[i] - 0.5*dv[i] )
		fs_r0[i]  = 1.*ps_r0[i+1]
		## get future active slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(matchpop)
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
		## past passive slice input params (mutually consistent)
		ps_u0[i]  = 1.*fs_u0[i] - 1.*du[i]
		ps_v0[i]  = 1.*fs_v0[i] - 1.*dv[i]
		ps_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[fs_u0[i]],[fs_v0[i]]]))[0]
		## get past passive slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(matchpop)
		print "i=%s pslice loc: %s"%(i,sliceloc)
		## execute past passive slice at sliceloc
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		##### iterate
		## iterate
		i -= 1
	## make sliceparams dict
	chainparams = dict(Rh=1.*np.array(Rh), ps_u0=1.*np.array(ps_u0), ps_v0=1.*np.array(ps_v0), ps_r0=1.*np.array(ps_r0), fs_u0=1.*np.array(fs_u0), fs_v0=1.*np.array(fs_v0), fs_r0=1.*np.array(fs_r0), i0=1*i0, matchmode=matchmode)
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
			ba = copy.deepcopy(b)
			bb = copy.deepcopy(b)
			bc = copy.deepcopy(b)
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
	## return
	return reglist, chainparams













##################################################################################################################



