

import numpy as np
import matplotlib.pyplot as plt


"""
Consider an evaporating black hole with maximum mass M.
Suppose it evaporates by emitting Nevap shells of Hawking radiation.
This is described by a sequence of 
    Nevap+1 regions (first with m=M, last with m=0) with parameters m_j, du_j, dv_j, rp_j, rf_j
    Nevap shells (future boundary of each region) with parameters dmi, r_i, u_i, v_i
Suppose each is Schwarz so that
 	F(r) = 2m*( (r/2m-1) + log(r/2m-1) ).
Evap length is le.

j
M  ...      ... 0
0 | 1 | 2 | ...
  0   1   2  ...
  i

Constraints:
	R_j  = 2m_j
	dm_i = m_i+1 - m_i
	du_j = u_j - u_j-1
	dv_j = v_j - v_j-1

	dm_i < 0
	du_j = -3T/M m_j^2 dm_j-1 > 0
	dv_j - du_j = 2( F(rf_j) - F(rp_j) )

	rf_j = Rh(m_j)+le
	rp_j = Rh(m_j-1)+le

"""
Nevap = 10
T = 10.
M = 1.
le = .1

def F(r,m):
	r[r<=2.*m] = np.nan
	return r + 2.*m * np.log(np.abs(r/(2.*m) - 1.))

def Rh(m):
	return 2.*m


m = np.nan * np.ones(Nevap+1)
du = np.nan*m
dv = np.nan*m

dv[1:-1] = .1*T
du[1.]





print "m  = %s"%(m)
print "du = %s"%(du)
print "dv = %s"%(dv)

