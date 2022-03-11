

import numpy as np
import matplotlib.pyplot as plt
import xhorizon as xh


func = xh.mf.schwarzschild()
reg = xh.reg.MAXreg(func)

xh.newfig()
reg.rplot()
plt.show()

