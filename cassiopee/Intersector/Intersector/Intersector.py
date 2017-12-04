"""Grid generation module.
"""
__version__ = '2.5'
__author__ = "Sam Landier, Christophe Benoit, Stephanie Peron, Luis Bernardos"
# 
# Python Interface to create arrays defining meshes
#
import intersector
import Generator as G
import Converter as C

def conformUnstr(a1, a2=None, tol=0., left_or_right=0, itermax=10):
    """Conformize a1 (optionally with a2).
    If a2 is specified, the third argument "left_or_right" tells wheter the ouput contains only a1 modified (0), a2 modified (1) or both (2).
    Usage: conformUnstr(a1, a2, tol, left_or_right)"""
    c = intersector.conformUnstr(a1, a2, tol, left_or_right, itermax)
    return G.close(c)

def intersection(a1, a2, tol=0.):
    """Compute the intersection between two input closed surfaces.
    Usage: intersection(a1, a2, tol)"""
    try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
    except: pass
    c = intersector.booleanIntersectionBorder(a1, a2, tol, 1, 1, 0) #last 3 args are dummy for now
    return G.close(c)

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the intersection between two input closed entities.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right, agg_mode)"""
    if a1[3] != 'NGON' and a2[3] != 'NGON':
      try:
          import Converter
          a1 = Converter.convertArray2Tetra(a1)
          a2 = Converter.convertArray2Tetra(a2)
          a1 = G.close(a1); a2 = G.close(a2)
      except: pass
    c = intersector.booleanIntersection(a1, a2, tol, preserve_right, solid_right, agg_mode)
    return G.close(c)

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, extrude_pgs=[]): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the union of two input closed surfaces.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""
    if a1[3] != 'NGON' and a2[3] != 'NGON':
      try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
      except: pass
      c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, extrude_pgs)
      return G.close(c)
    else: 
      c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, extrude_pgs)
      return c #close is done inside

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the difference between the two input closed surfaces.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    if a1[3] != 'NGON' and a2[3] != 'NGON':
      try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
      except: pass
    c = intersector.booleanMinus(a1, a2, tol, preserve_right, solid_right, agg_mode)
    return G.close(c)
    
def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    c = intersector.booleanModifiedSolid(a2, solid, tol, 1, preserve_solid, agg_mode)
    return G.close(c)
    
#==============================================================================
# XcellN
# IN: coords: 3D structured or unstructured mesh
# OUT: returns the cellnfields, 0 for fully inside, 1 for fully outside, in between when intersecting
#==============================================================================
def XcellN(coords, cellnfields, maskingMesh, wall_pgl=[], ghost_pgl=[]):
    """Compute the acurate cellN based on overlapped volumes.
    Usage: XcellN(coords, cellnfields, maskingMesh)"""
    cellnt = []
    #C.convertArrays2File([maskingMesh], "mask.plt")
    #print pgl
    for i in xrange(len(coords)):
      #print 'coords : %d / %d' %(i+1, len(coords))
      #C.convertArrays2File([coords[i]], "bloc%d.plt"%(i))
      #print pgl
      #print coords[i]
      cn = intersector.XcellN(coords[i], cellnfields[i], maskingMesh, wall_pgl, ghost_pgl)
      #C.convertArrays2File([cn], "walls%d.plt"%(i))
      cellnt.append(cn)
    return cellnt
    
#==============================================================================
# P1ConservativeChimeraCoeffs
# IN: aR : receiver mesh
# IN: cellnR : receiver cellN (only cells with value equal to 2 will be considered)
# IN: aD : donor (source) mesh
# OUT: (indices, coeffs, delimiter, receiver original ids)
#==============================================================================
def P1ConservativeChimeraCoeffs(aR, cellnR, aD):
    return intersector.P1ConservativeChimeraCoeffs(aR, cellnR, aD)

#==============================================================================
# triangulateExteriorFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(mesh):
    return intersector.triangulateExteriorFaces(mesh)
    
#==============================================================================
# convexifyFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the faces made convex
#==============================================================================
def convexifyFaces(mesh, convexity_TOL = 1.e-8):
    return intersector.convexifyFaces(mesh, convexity_TOL)

#==============================================================================
# prepareCellsSplit
# IN : mesh         : 3D NGON mesh
# IN : PH_set       : PH to process. 0 for concave cells or 1 for non-centroid-star_shaped cells
# IN : split_policy : 0 : convexify concave pgs. 1 : starify concave pgs. 2 : starify any pgs at concave-chains ends.
# OUT: returns a 3D NGON Mesh with some face polygon splits : 
#      split (convexify, starify) some targeted polygons on targeted cells 
#      (typically bad polyhedra -concaves, non-centroid-star-shaped-)
#      to prepare the split of those bad cells.
#==============================================================================
def prepareCellsSplit(mesh, PH_set = 1, split_policy = 0, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8):
    return intersector.prepareCellsSplit(mesh, PH_set, split_policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)
   
 #==============================================================================
# XXX
#==============================================================================
def splitNonStarCells(mesh, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8):
    return intersector.splitNonStarCells(mesh, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)

#==============================================================================
# simplifyCells : agglomerate superfluous polygons that overdefine cells
# IN: mesh: 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifyCells(mesh, treat_externals, angular_threshold = 1.e-12):
    return intersector.simplifyCells(mesh, treat_externals, angular_threshold)

#==============================================================================
# agglomerateSmallCells : XXX
#==============================================================================
def agglomerateSmallCells(mesh, vmin=0., vratio=1000.):
    return intersector.agglomerateSmallCells(mesh, vmin, vratio)

#==============================================================================
# agglomerateSmallCells : XXX
#==============================================================================
#def agglomerateUncomputableCells(mesh):
#    return intersector.agglomerateUncomputableCells(mesh)

#==============================================================================
# agglomerateSmallCells : XXX
#==============================================================================
def collapseUncomputableFaces(mesh):
    return intersector.collapseUncomputableFaces(mesh)

    
#==============================================================================
# extractUncomputables : XXX
#==============================================================================
def extractUncomputables(mesh, neigh_level=2):
    return intersector.extractUncomputables(mesh, neigh_level)

#==============================================================================
# extractPathologicalCells : XXX
#==============================================================================
def extractPathologicalCells(mesh, concave_threshold=1.e-13, neigh_level=0):
    return intersector.extractPathologicalCells(mesh, concave_threshold, neigh_level)

#==============================================================================
# extractOuters : XXX
#==============================================================================
def extractOuterLayers(mesh):
    return intersector.extractOuterLayers(mesh)

#==============================================================================
# XXX : XXX
#==============================================================================
def extractNthCell(mesh, nth):
    return intersector.extractNthCell(mesh, nth)

#==============================================================================
# XXX : XXX
#==============================================================================
def extractNthFace(mesh, nth):
    return intersector.extractNthFace(mesh, nth)

#==============================================================================
# XXX : XXX
#==============================================================================
def removeNthCell(mesh, nth):
    return intersector.removeNthCell(mesh, nth)

#==============================================================================
# diffMesh : returns the diff between 2 meshes as 2 zones
#==============================================================================
def diffMesh(mesh1, mesh2):
    return intersector.diffMesh(mesh1, mesh2)

#==============================================================================
# checkCellsClosure : returns any open cells (those with free edges) in the input mesh
#==============================================================================
def checkCellsClosure(mesh):
    return intersector.checkCellsClosure(mesh)

#==============================================================================
# selfX : returns any intersecting cells in the input mesh
#==============================================================================
def selfX(mesh):
    return intersector.selfX(mesh)

def extrudeUserDefinedBC(mesh, extrude_pgs=[]):
    return intersector.extrudeUserDefinedBC(mesh, extrude_pgs)

def reorientExternalFaces(mesh):
  return intersector.reorientExternalFaces(mesh)


#~ def conservativeTransfer(a1, flowsol, a2, tol=0., reconstruction_type=0):
    #~ c = intersector.conservative_transfer(a1, flowsol, a2, tol, reconstruction_type)
    #~ return c
    #~ 
#~ def totalMass(a1, flowsol):
    #~ intersector.total_mass(a1, flowsol)
    #~ return a1
    
