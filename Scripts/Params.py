import string


# CELL PARAMETERS

DURATION = {
    'orbit:bx': 3564,
    'orbit': 3564*25,
    'bx': 25.,
    'tdc': 25./30
}

XCELL = 42. # cell width in mm
ZCELL = 13. # cell height in mm

WIRE_DIAM     = 0.050 # in mm
PLANE_WIDTH   = 1.5 # in mm
IBEAM_WIDTH   = 1.3  # in mm FROM https://github.com/cms-sw/cmssw/blob/master/Geometry/DTGeometry/src/DTTopology.cc
IBEAM_WING    = 6.35  # in mm FROM https://github.com/cms-sw/cmssw/blob/master/Geometry/DTGeometry/src/DTTopology.cc

TM         = 15.5
TDRIFT     = TM*DURATION['bx']    # drift time in ns
VDRIFT     = XCELL*0.5 / TDRIFT   # drift velocity in mm/ns
VDRIFTMMBX = XCELL*0.5 / TM       # drift velocity in mm/BX
VHRATIO    = XCELL*0.5/TM/ZCELL

## number of cells
NLAYERS = 4 # number of layers
NWIRES  = 4 # numbers of cells per layer
# Starting from a single cell, new columns are added alternating on the right and on the left,
# while new rows are added alternating on the bottom and on the top

## shifts
pos_shift_z  = [ZCELL*(-(NLAYERS-1-(NLAYERS%2))/2 + i) for i in range(NLAYERS)] # [..., -1.5*ZCELL, -0.5*ZCELL, 0.5*ZCELL, 1.5*ZCELL, ...]
                                                                                # If NLAYERS is odd, the additional layer in on top
is_shifted_right = [(i+1+(NLAYERS//2)%2)%2 for i in range(NLAYERS)]             # 1 if layer is shifted to the right, 0 otherwise
                                                                                # For NLAYERS=4 it is [1, 0, 1, 0]

# mappings

WIRE_MAP = {i:string.ascii_uppercase[i-1] for i in range(1, NWIRES+1)}
SIDE_MAP = {-1:'L',1:'R'}

LAYER_MAP = [-99] + [i-1 for i in range(NLAYERS, 0, -1)]


cell_ineff = 0.06