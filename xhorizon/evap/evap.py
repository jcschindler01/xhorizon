
"""
This module provides method for making forming and evaporation BH diagrams.
"""

import numpy as np
import matplotlib.pyplot as plt


def accrete(reg1, v1=0., v2=0., R2=1.):
	"""
	Accrete a Hayward of outer radius R2 onto the region reg1 on a line of constant v by 
	slicing reg1 at v1 and the new region reg2 at v2.

	Inputs:
		reg1 = existing region to accrete onto
		v1   = value of v to slice reg1
		v2   = value of v to slice new region
		R2   = Hayward outer radius of new region

	Returns:
		reg2 = New Hayward region

	This method actively edits the object reg1, so reg1 should be different after running the 
	method even though it is not returned.
	"""
	return 0





