

"""
This file imports subroutines and modules into the xhorizon namespace.
Only functions and modules needed by the user should be imported here.
"""

"""
MODULE DEPENDENCY TREE

xhorizon2
	diagram_tools
		region_factory < curve_class
		region_factory < curvemakers
		region_factory < diagram_classes
		diagram_classes < coord_transf
		diagram_classes < block_masks
		diagram_classes < block_fill
		block_fill < curvemakers
		block_fill < curve_class
		block_fill < block_masks
		curvemakers < block_masks
		curvemakers < curve_class
		block_masks < curve_class
		coord_transf
		curve_class
	tortoise
		metfunc_factory < metfunc_class
		metfunc_factory < tortoise
		metfunc_factory < math_util
		metfunc_class
		tortoise
		math_util
		metfunc_tests
	shell_junction
		slicers < interpolators
		reg_corner_masks
		interpolators
"""

from diagram_tools.plotstyle import *

from diagram_tools.curve_class import curve
from diagram_tools.diagram_classes import diagram, region, block

from tortoise.metfunc_tests import func_test

import diagram_tools.region_factory as reg
import tortoise.metfunc_factory as mf

import diagram_tools.block_masks as mask
import diagram_tools.curvemakers as cm
import diagram_tools.block_fill as fill

import shell_junction as junc
import shell_junction.reg_corner_masks as cornermask

import evap


