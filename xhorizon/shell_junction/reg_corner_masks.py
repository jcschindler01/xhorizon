
"""
This module provides functions for applying uvbounds masks to correct blocks in a given region.
These functions should directly edit and return the input region object.
"""

import numpy as np



def EFreg2a(reg, abcd=None, u0=None, v0=None):
	"""
	For Hayward-like regions (aka N=2 and f(0)>0).
	Assumes slice in outermost block.
		reg = The region to mask.
		abcd = Which corner junction label is this piece, with standard label scheme. 
			Should be a single letter string 'a' 'b' 'c' or 'd'.
		u0, v0 = Location to slice.
	Returns the same region that was input, after editing.

	Region is masked by removing unwanted blocks and setting uvbounds on remaining blocks.

	Block order is [inner, middle, outer] (this is arbitrary, just happens to be how EFreg2a is generated).
	Setting any of the uvbounds to np.nan causes all inequalities to fail, which causes an error unless block is removed.
	"""
	## initialize
	uvb = [dict(), dict(), dict()]
	## set proper uvbounds for each block based on abcd type
	## case a
	if abcd=='a':
		uvb = [dict(vmin=v0), dict(vmin=v0), dict(vmin=v0, umin=u0)]
	## case b
	if abcd=='b':
		uvb = [dict(vmin=np.nan), dict(vmin=np.nan), dict(vmax=v0,umax=u0)]
	## case c
	if abcd=='c':
		uvb = [dict(vmin=np.nan), dict(vmin=np.nan), dict(vmin=v0,umax=u0)]
	## case d
	if abcd=='d':
		uvb = [dict(vmax=v0), dict(vmax=v0), dict(vmax=v0, umin=u0)]
	## update blocks uvbounds
	for i in range(len(reg.blocks)):
		reg.blocks[i].uvbounds.update(uvb[i])
	## keep blocks only if no nan in uvbounds values
	keep = []
	for b in reg.blocks:
		if not np.nan in b.uvbounds.values():
			keep += [b]
	reg.blocks = keep
	## return
	return reg


def MAXreg2a(reg, abcd=None, u0=None, v0=None):
	"""
	For Schwarzschild-like regions (aka N=2 and f(0)>0).
		reg = The region to mask.
		abcd = Which corner junction label is this piece, with standard label scheme. 
			Should be a single letter string 'a' 'b' 'c' or 'd'.
		u0, v0 = Location to slice.
	Returns the same region that was input, after editing.

	Region is masked by removing unwanted blocks and setting uvbounds on remaining blocks.

	Block order is [top, right, bottom, left] (this is arbitrary, just happens to be how MAXreg2a is generated).
	Setting any of the uvbounds to np.nan causes all inequalities to fail, which causes an error unless block is removed.
	"""
	## initialize
	uvb = [dict(), dict(), dict(), dict()]
	## set proper uvbounds for each block based on abcd type
	## case a
	if abcd=='a':
		uvb = [dict(vmin=v0), dict(vmin=v0,umin=u0), dict(vmin=np.nan), dict(vmin=np.nan)]
	## case b
	if abcd=='b':
		uvb = [dict(vmin=np.nan), dict(vmax=v0,umax=u0), dict(umax=u0), dict(vmin=np.nan)]
	## case c
	if abcd=='c':
		uvb = [dict(vmin=np.nan), dict(vmin=v0,umax=u0), dict(vmin=np.nan), dict(vmin=np.nan)]
	## case d
	if abcd=='d':
		uvb = [dict(vmax=v0), dict(vmax=v0,umin=u0), dict(umin=u0), dict()]
	## update blocks uvbounds
	for i in range(len(reg.blocks)):
		reg.blocks[i].uvbounds.update(uvb[i])
	## keep blocks only if no nan in uvbounds values
	keep = []
	for b in reg.blocks:
		if not np.nan in b.uvbounds.values():
			keep += [b]
	reg.blocks = keep
	## return
	return reg








