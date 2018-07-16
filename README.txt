
-------------------------------------------------
XHORIZON
Jan 2018
Joseph C Schindler
jcschind@ucsc.edu
https://github.com/jcschindler01/xhorizon
-------------------------------------------------

-------------------------------------------------
XHORIZON is a python package for the explicit computation and drawing of Penrose diagrams.

Diagrams can be generated numerically and plotted for any spacetime in the following classes:
 (A) "Strongly spherically symmetric" spacetimes, which have a metric of the form  
     ds^2 = - f(r) dt^2 + f(r)^{-1} dr^2 + r^2 d\Omega^2.
 (B) Spacetimes which are composed piecewise from an arbitrary finite number of SSS regions
     joined across null shell junctions.
This includes the Minkowski, Schwarzschild, Reissner-Nordstrom, 2D Kerr, de Sitter, Anti de 
Sitter, and Schwarzschild-dS metrics, among many others.

The method is based on the article:

"Algorithms for the explicit computation of Penrose diagrams"
JC Schindler, A Aguirre
Jan 2018
arXiv:1802.02263

https://arxiv.org/abs/1802.02263
https://doi.org/10.1088/1361-6382/aabce2
-------------------------------------------------

-------------------------------------------------
INSTALLATION AND USAGE

To install:
 - Download a .zip copy of (the desired version of) the directory containing this file from 
         https://github.com/jcschindler01/xhorizon
 - Unpack the directory to any location.
 - Make the directory available to your python path 
         (e.g. by editing your PYTHONPATH environment variable).
 - Done! The xhorizon tools are now available for import like any other package.

To use:
 - Import the package in python using standard import commands:
         $ python
         >>> import xhorizon
 - All set! Call xhorizon's tools from the package as necessary.
-------------------------------------------------

-------------------------------------------------
DOCUMENTATION AND TUTORIAL

Docs:
There is currently no documentation external to the source code.
However, all modules contain detailed commentary and documentation in the source.
We hope to provide external documentation soon!

Tutorial:
A tutorial and set of examples should be provided soon. Stay tuned!
-------------------------------------------------

-------------------------------------------------
DEVELOPMENT INFO

Developed under:

PYTHON 2.7.9
NUMPY 1.9.1
MATPLOTLIB 1.4.2
SCIPY 0.15.0

OS: WINDOWS 7 64 BIT
PYTHON DISTRO: WINPYTHON 64BIT 2.7.9.2
-------------------------------------------------

-------------------------------------------------
COMPATIBILITY NOTES

Plot commands are styled in MatPlotLib 1 style.
For use with MatPlotLib 2, please use
	plt.style.use('classic')
and be aware of possible color compatibility issues.
-------------------------------------------------

-------------------------------------------------
MODULE DEPENDENCY TREE

See xhorizon.__init__.py.
-------------------------------------------------








