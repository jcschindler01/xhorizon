
"""
Example of generating a Penrose diagram with xhorizon.
"""

# Imports.
import numpy as np
import matplotlib.pyplot as plt
import xhorizon as xh


# Set up a nicely styled new figure.
xh.newfig()

# Create a "metfunc" object. 
# Try adjusting the charge parameter Q to a lower value or increasing the mass M.
f = xh.mf.reissner_nordstrom(M=0.5, Q=0.3)

# Create a "region" object.
# Try changing the boolean keyword inputs to see effect. Or, try changing MAXreg (maximally extended) to EFreg (Eddington-Finklestein).
reg = xh.reg.MAXreg(f, rlines=True, boundary=True, uvgrid=False)

# Manually add additional curves to the plot. Easy to turn off by changing "if True" to "if False".
if False:
  # Radius and time values in Schwarzschild (t,r) coordinates, parameterized by s. 
  s = np.linspace(-1,3,5001)
  tt = 2.*np.sinh(s)
  rr = 1.1*np.cosh(s)
  # Create new "curve" object, load (t,r) values, set style.
  c = xh.curve()
  # Load (t,r) Schwarzschild coordinate values.
  c.tr = np.array([tt,rr])
  # Set style.
  c.sty = dict(c='r', zorder=700, lw=1.5)
  # Add curve to all blocks in region. Automatically not added if curve not in region, so the above curve will only appear in exterior regions for R=1 Schwarzschild. Input to add_curves_tr must be a list of curve objects.
  for b in reg.blocks:
    b.add_curves_tr([c])

# Manually add additional curves to the plot.
## Initialize list of curves to add.
crvlist = []
## Let's add some lines of t=constant in the exterior region.
T = np.linspace(-10,10,15)
for t0 in T:
  ## Create new "curve" object.
  crv = xh.curve()
  ## The Schwarzchild coordinates (t,r) of each curve will be parameterized by a parameter s.
  s = np.linspace(1e-9,10,1001)
  ## Define t and r coordinate arrays.
  t = t0 + 0.*s            ## t is a constant array of correct length
  r = f.rj[-2] + 1.*s      ## f.rj[-2] is the outer horizon radius
  ## Add coordinates to the curve object.
  crv.tr = np.array([t,r])
  ## Set curve display style. In this case red, .5 units wide, and behind rlines.
  crv.sty.update(dict(c='r', lw=.5, zorder=300))
  ## Add curve to curve list.
  crvlist += [crv]
## Add curves to all blocks in region. Will be rejected from incompatible blocks.
for b in reg.blocks:
  b.add_curves_tr(crvlist)

# Fill trapped regions (that is, blocks where sign of f is -1) with gray background.
for b in reg.blocks:
  if b.sgnf==-1.:
    b.fill(sty=dict(fc='k', alpha=0.2))

# Tell the region object to plot itself. Always do this *after* adding any additional curves to the region.
reg.rplot()

# Add an annotation with the metric function info.
# Note that the label is made by extracting strings from the dictionary f.info and separating with linebreaks.
nice_label = '\n'.join([f.info['Type'], f.info['Metric Function'], f.info['Parameters']])
plt.annotate(s=nice_label, xy=(.035,.03), xycoords='axes fraction', ha='left', va='bottom', size=8)

# Save the plot.
plt.savefig("out.png", dpi=300)




