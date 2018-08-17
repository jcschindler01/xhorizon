
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



def funclist_chain(funclist, seed=0, u0=0., v0=0., du=None, dv=None, mu=0.):
	"""
	Create a chain of matched regions, starting at seed region which is unmodified.
	Each region except ends has two slices through it, a future slice fslice and past slice pslice.
	Each fslice and pslice can be either active or passive, but there can only be one active slice per region.
	Return reglist.
	The index i refers to each region in the sequence for all variables.
	If mu=0 then all blocks are centered about u=0 symettrically. If mu=1 then u,v are connected from block to block.
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
	## seed region
	i = 1*i0
	ps_u0[i]  = 1.*mu*u0 - 0.5*(1.-mu)*du[i]
	ps_v0[i]  = 1.*mu*v0 - 0.5*(1.-mu)*dv[i]
	ps_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[ps_u0[i]],[ps_v0[i]]]))[0]
	fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
	fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
	fs_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[fs_u0[i]],[fs_v0[i]]]))[0]
	pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), r0=1.*ps_r0[i], u0=1.*ps_u0[i])
	fslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), r0=1.*fs_r0[i], u0=1.*fs_u0[i])
	## forward slices
	i = 1*i0 + 1
	while i < len(reglist):
		## past active slice
		ps_u0[i]  = 1.*mu*fs_u0[i-1] - 0.5*(1.-mu)*du[i]
		ps_v0[i]  = 1.*mu*fs_v0[i-1] - 0.5*(1.-mu)*dv[i]
		ps_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[ps_u0[i]],[ps_v0[i]]]))[0]
		pslice[i] = xh.junc.aslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), r0=1.*ps_r0[i], u0=1.*ps_u0[i], U0=fslice[i-1].U_of_r_at_v0, V0=fslice[i-1].V_of_r_at_u0, r_refs=[fslice[i-1].reg.metfunc.r_ref])
		## adjust transformations
		reglist[i].U_of_udl = pslice[i].U_of_udl_at_v0
		reglist[i].V_of_vdl = pslice[i].V_of_vdl_at_u0
		## future passive slice
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[fs_u0[i]],[fs_v0[i]]]))[0]
		fslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), r0=1.*fs_r0[i], u0=1.*fs_u0[i])
		## iterate
		i += 1
	## backward slices
	i = 1*i0 - 1
	while i >= 0:
		## future active slice
		fs_u0[i]  = 1.*mu*ps_u0[i+1] + 0.5*(1.-mu)*du[i]
		fs_v0[i]  = 1.*mu*ps_v0[i+1] + 0.5*(1.-mu)*dv[i]
		fs_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[fs_u0[i]],[fs_v0[i]]]))[0]
		fslice[i] = xh.junc.aslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), r0=1.*fs_r0[i], u0=1.*fs_u0[i], U0=pslice[i+1].U_of_r_at_v0, V0=pslice[i+1].V_of_r_at_u0, r_refs=[pslice[i+1].reg.metfunc.r_ref])
		## adjust transformations
		reglist[i].U_of_udl = fslice[i].U_of_udl_at_v0
		reglist[i].V_of_vdl = fslice[i].V_of_vdl_at_u0
		## past passive slice
		ps_u0[i]  = 1.*fs_u0[i] - 1.*du[i]
		ps_v0[i]  = 1.*fs_v0[i] - 1.*dv[i]
		ps_r0[i]  = 1.*reglist[i].blocks[-1].r_of_uv(np.array([[ps_u0[i]],[ps_v0[i]]]))[0]
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), r0=1.*ps_r0[i], u0=1.*ps_u0[i])
		## iterate
		i -= 1
	## print
	print ps_u0
	print fs_u0		
	## return
	return reglist






##################################################################################################################



