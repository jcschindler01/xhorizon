

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
	r = r*np.array([1.])
	if np.any(r<=2.*m):
		out = None
	else:
		out = r + 2.*m * np.log(np.abs(r - 2.*m))
	return 1.*out

def R(m):
	return 2.*m

def zero(m, mprev=1., dv=1., T=1., M=1., l=1.):
	dv = 1.*dv
	du = -(3.*T/M**3) * m**2 * (m-mprev)
	drstar = F(R(m)+l,m) - F(R(mprev)+l,m)
	#return dv-du-2.*drstar
	return dv-2.*drstar, du



m = np.nan * np.ones(Nevap+1)
du = np.nan*m
dv = np.nan*m


m[0] = 1.*M
dv[1:-1] = .1*T



for mprev in np.linspace(1.,0.,5)[:-1]:

	M = .01

	T  = 50.   *M
	l  = .01   *M
	dv = .5   *M
	mp = mprev*M

	params = dict(mprev=mp, dv=dv, T=T, M=M, l=l)
	mm = np.linspace(0,mp,1001)
	zz = zero(mm, **params)

	plt.plot(mm, zz[0], 'r-')
	plt.plot(mm, zz[1], 'b-')

	plt.plot(mm, 0.*mm, 'k-')
	plt.xlim(0,M)
	plt.ylim(np.array([-1.,1.])*np.max(np.abs(plt.ylim())))
	plt.show()









print "m  = %s"%(m)
print "du = %s"%(du)
print "dv = %s"%(dv)

