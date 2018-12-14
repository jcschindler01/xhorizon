

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


def F(r,m):
	return r + 2.*m * np.log(np.abs(r/(2.*m) - 1))

def R(m):
	return 2.*m

def zero(m, mprev=1., dv=1., T=1., M=1., l=1.):
	dv = 1.*dv
	du = -(3.*T/M**3) * m**2 * (m-mprev)
	drstar = F(R(m)+l,m) - F(R(mprev)+l,m)
	#return dv-du-2.*drstar
	return dv-2.*drstar, du


Nevap = 10
T = 100.
M = .1
l = .1

## u
u = np.linspace(1e-9,T-1e-9,Nevap)

## m(u)
m = M * (1. - u/T)**(1./3.)

## r(u)
r = 2.*m + l

## v(u)
v = u + 2.*F(r,m)

## du dv dm
du = u[1:]-u[:-1]
dv = v[1:]-v[:-1]
dm = m[1:]-m[:-1]

## mself mprev
mself = m[1:]
mprev = m[:-1]

## drstar
drstar = F(R(mself)+l,mself) - F(R(mprev)+l,mself)


print np.mean(drstar)


