
"""
Provides methods for drawing forming and evaporating black hole diagrams.
"""

import numpy as np
import matplotlib.pyplot as plt

import xhorizon as xh
from evap import *
from accrete import *
from helpers import *




def formevap(accrete_R=[1.], accrete_v=[-1.], evap_R=[0.], evap_u=[1.], l=0.1, rparams=dict(c=0., s0=10.)):
	"""
	Form and evaporate a black hole by specifying accretion radius as a function of time and evaporation
	radius as a function of time, where time is measured by a faraway observer.
	"""
	## reglist
	reglist = []
	## accretion
	reglist += accretion(accrete_R, accrete_v, l=l, rparams=rparams)
	## evaporation
	reglist += evaporation(evap_R, evap_u, reglist.pop(), l=l, rparams=rparams)
	## return
	return reglist




