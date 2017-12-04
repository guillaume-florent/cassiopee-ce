"""Grid generation module.
"""
# 
# Python Interface to create PyTrees defining meshes
#
import Intersector as XOR
import intersector

__version__ = XOR.__version__

import numpy

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
except:
    raise ImportError("Intersector.PyTree: requires Converter.PyTree module.")
        
#=============================================================================
# Ajout du cellN pour une zone qui ne le contient pas. celln=1
#=============================================================================
def addCellN__(t, loc='centers'):
    tp = Internal.copyRef(t)
    _addCellN__(tp, loc)
    return tp

def _addCellN__(t, loc='centers'):
    if loc == 'centers': var = loc+':cellN'
    else: var = 'cellN'
    C._addVars(t, var)
    return None

def nb_cells(a):
  import Converter.Internal as I
  ncellsTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
      dim = I.getZoneDim(z)
      np = dim[1]
      ncells = dim[2]
      ncellsTot += ncells
  return ncellsTot

#=============================================================================
# Concatenation des PointList d un type de BC donne dans une liste de zones
#=============================================================================
def concatenateBC(bctype, zones, wallpgs, cur_shift):
    i=0
    for z in zones:
      c = C.getFields(Internal.__GridCoordinates__, z)

      if (c == []): continue

      #print ' -- zone : %d / %d' %(i+1, len(zones))
      i=i+1
      bnds = Internal.getNodesFromType(z, 'BC_t')
      #print " -- this zone has %d boundaries"%(len(bnds))
      #print ' -- cur shift %d' %(cur_shift)

      # GET THE WALL PGS FROM THE POINTLISTS
      for bb in bnds :
        #print bb
        if (Internal.isValue(bb, bctype) == False) : continue
          
        #print bb[1]#[2][1]
        wpgs = bb[2][1][1][0] # POINTLIST NUMPY
        #print wpgs
        # SYNC THE POINTLIST BEFORE APPENDING  : SHIFT WITH THE CURRENT NB OF STORED POLYGONS
        id2 = numpy.empty(len(wpgs), numpy.int32)
        id2[:] = wpgs[:] + cur_shift
        #print id2
        wallpgs.append(id2)

      c = c[0]
      #print c
      #z_nb_pts= len(c[1][0])
      z_nb_pgs= c[2][0][0]
      #print z_nb_pts
      #print z_nb_pgs
      cur_shift += z_nb_pgs
    return (wallpgs, cur_shift)
  
#------------------------------------------------------------------------------
# Conformisation d'une soupe de TRI ou de BAR
#------------------------------------------------------------------------------
def conformUnstr(surface1, surface2=None, tol=0., left_or_right=0, itermax=10):
    """Conformize a TRI or BAR soup (surface1) with taking into account surface2 if it's provided.
    Usage: conformUnstr(s1, s2, tol, left_or_right, itermax)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    if surface2 is not None:
        s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    else: s2 = None
    s = XOR.conformUnstr(s1, s2, tol, left_or_right, itermax)
    return C.convertArrays2ZoneNode('conformized', [s])

def intersection(surface1, surface2, tol=0.):
    """Compute the intersection between two input closed surfaces.
    Usage: intersection(s1, s2, tol)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    s = XOR.intersection(s1, s2, tol)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the intersection between two input closed entities.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanIntersection(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the union between two input closed entities.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]

    cur_shift=0
    extrudepgs=[]
    if (solid_right == 1) :
        zones = Internal.getZones(a2)
        (extrudepgs, cur_shift) = concatenateBC('BCUserDefined', zones, extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    #print "nb of pgs to pass : %s" %(len(extrudepgs))

    s = XOR.booleanUnion(s1, s2, tol, preserve_right, solid_right, agg_mode, extrudepgs)
    return C.convertArrays2ZoneNode('union', [s])

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the difference between the two input closed surfaces.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanMinus(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('minus', [s])
    
def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanModifiedSolid(a1, a2, tol, preserve_right, solid_right)"""
    sld = C.getFields(Internal.__GridCoordinates__, solid)[0]
    operand = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanModifiedSolid(operand, sld, tol, preserve_solid, agg_mode)
    return C.convertArrays2ZoneNode('modified_solid', [s])
    
#==============================================================================
# XcellN
# IN: t: background Mesh (NGON 3D)
# IN: prioritaryMesh: hiding Mesh (NGON 3D)
# IN: blankingMatrix
# OUT: returns the cellnfields, between 0 (fully hidden) and 1 (fully visible)
#==============================================================================
def XcellN(t, prioritaryMesh, blankingMatrix=[]):
    try: import Transform as T
    except: raise ImportError("XcellN: requires Transform module.")

    nb = -1
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    a = addCellN__(a, loc=loc)
    bases = Internal.getBases(a)
    if blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(prioritaryMesh)), numpy.int32)
    for b in bases:
        nb += 1
        #print 'bgm base : %d / %d' %(nb+1, len(bases))
        coords = C.getFields(Internal.__GridCoordinates__, b)
        if coords == []: continue

        coords = Converter.convertArray2NGon(coords)

        if loc == 'centers': cellN = C.getField('centers:cellN', b)
        else: cellN = C.getField('cellN', b)

        bc = []
        wallpgs = [] # LIST OF BODY WALLS IDS USED TO IGNORE BGM CELLS INSIDE THEM 
        ghostpgs = [] # LIST OF BOUNDARIES TO EXTRUDE TO PREVENT UNECESSARY X COMPUTATIONS
        cur_shift=0
        for nb2 in xrange(len(prioritaryMesh)):
            blanking = blankingMatrix[nb, nb2]
            #if (prioritaryMesh[nb2] == []): print 'empty'
            if (prioritaryMesh[nb2] == []): continue
            
            #print 'hiding base : %d / %d' %(nb2+1, len(prioritaryMesh))
            zones = Internal.getZones(prioritaryMesh[nb2])
            i=0
            for z in zones:
                c = C.getFields(Internal.__GridCoordinates__, z)

                if c == []: continue

                #print ' -- hiding base %d zone : %d / %d' %(nb2+1, i+1, len(zones))

                c = c[0]
                bc.append(c)

            (wallpgs, cur_shift_new) = concatenateBC('BCWall', zones, wallpgs, cur_shift)

            (ghostpgs, cur_shift_new) = concatenateBC('BCUserDefined', zones, ghostpgs, cur_shift)
            cur_shift=cur_shift_new

        if (wallpgs != []) : wallpgs = numpy.concatenate(wallpgs) # create a single list
        #print "nb of wall pgs %s"%(len(wallpgs))
        if (ghostpgs != []) : ghostpgs = numpy.concatenate(ghostpgs) # create a single list
        #print "nb of ghost pgs %s"%(len(ghostpgs))
                
        if bc == []:
            #print 'Warning : no xcelln to compute for base %d'%(nb)
            continue

        bc = Converter.convertArray2NGon(bc); bc = T.join(bc);

        cellN = XOR.XcellN(coords, cellN, bc, wallpgs, ghostpgs)
        bc = None
        coords = None
        C.setFields(cellN, b, loc, False)
    return a

#==============================================================================
# triangulateExteriorFaces
# IN: mesh: 3D NGON mesh
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(t):
    """Triangulate the exterior faces of a mesh.
    Usage: triangulateExteriorFaces(t)"""
    return C.TZA(t, 'nodes', 'nodes', XOR.triangulateExteriorFaces, None)

def _triangulateExteriorFaces(t):
    return C._TZA(t, 'nodes', 'nodes', XOR.triangulateExteriorFaces, None)

#==============================================================================
# convexifyFaces
# IN: mesh: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the faces made convex
#==============================================================================
def convexifyFaces(mesh, convexity_TOL=1.e-8):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.convexifyFaces(m)
    return C.convertArrays2ZoneNode('allPGconvex', [m])
    
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
def prepareCellsSplit(mesh, PH_set=1, split_policy=0, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold = 1.e-8):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.prepareCellsSplit(m, PH_set, split_policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)
    return C.convertArrays2ZoneNode('cvxSplitReady', [m])
    
    
#==============================================================================
# simplifyCells : agglomerate superfluous polygons that overdefine cells
# IN: mesh: 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifyCells(mesh, treat_externals, angular_threshold = 1.e-12):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.simplifyCells(m, treat_externals, angular_threshold)
    return C.convertArrays2ZoneNode('simplifiedCells', [m])


#==============================================================================
# XXX
#==============================================================================
def splitNonStarCells(mesh, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.splitNonStarCells(m, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)
    return C.convertArrays2ZoneNode('simplifiedCells', [m])

#==============================================================================
# agglomerateSmallCells : XXX
#==============================================================================
def agglomerateSmallCells(mesh, vmin=0., vratio=1000.):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]

    res = XOR.agglomerateSmallCells(m, vmin, vratio)
    #print "NB ZONES %d"%(len(res))

    z = C.convertArrays2ZoneNode('agglomeratedCells', [res[0]])

    debug = False

    if (debug == False) : return z

    zones = []
        
    nb_zones = len(res)-1
    if (nb_zones == 0) : return z
    #print nb_zones

    for i in xrange(nb_zones):
        zones.append(C.convertArrays2ZoneNode('agg', [res[i+1]]))

    C.convertPyTree2File(zones, 'agglo.cgns')

    return z

# def agglomerateSmallCells(mesh, vmin=0., vratio=1000.):
#     m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
#     print "one"
#     res = XOR.agglomerateSmallCells(m, vmin, vratio)
#     print "NB ZONES %d"%(len(res))

#     z = C.convertArrays2ZoneNode('agglomeratedCells', [res[0]])

#     debug = True

#     if (debug == False) : return z

#     nb_aggs = res[1]
#     nb_cels = nb_cells(z);
#     nb_points = len(res[0][1][0])
    
#     #print "NB AGG OVER NB CELLS : %d / %d "%(nb_aggs, nb_cels)

#     z= C.initVars(z, 'centers:cellN', 0)

#     cellN = Internal.getNodesFromName(z, 'cellN')[0][1]
#     cellN[0:nb_aggs] = numpy.arange(1,nb_aggs+1)

#     C.convertPyTree2File(z, 'agglo.cgns')

#     return z

#==============================================================================
# agglomerateSmallCells : XXX
#==============================================================================
# def agglomerateUncomputableCells(mesh):
#     m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
#     m = XOR.agglomerateUncomputableCells(m)
#     return C.convertArrays2ZoneNode('agglomeratedCells', [m])

#==============================================================================
# agglomerateSmallCells : XXX
#==============================================================================
def collapseUncomputableFaces(mesh):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.collapseUncomputableFaces(m)
    return C.convertArrays2ZoneNode('agglomeratedCells', [m])
    
    
#==============================================================================
# extractUncomputables : XXX
#==============================================================================
def extractUncomputables(mesh,neigh_level=2):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    res = XOR.extractUncomputables(m, neigh_level)
    
    zones = []
    nb_zones = len(res)
    
    if (nb_zones == 1) :
      zones.append(C.convertArrays2ZoneNode('remaining', [res[0]]))
    else:
      zones.append(C.convertArrays2ZoneNode('upgs', [res[0]]))
      zones.append(C.convertArrays2ZoneNode('uphs', [res[1]]))
      zones.append(C.convertArrays2ZoneNode('uphs_wv1', [res[2]]))
      zones.append(C.convertArrays2ZoneNode('remaining', [res[3]]))

    return zones
    
#==============================================================================
# extractPathologicalCells : XXX
#==============================================================================
def extractPathologicalCells(mesh, concave_threshold=1.e-13, neigh_level=0):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    res = XOR.extractPathologicalCells(m, concave_threshold, neigh_level)
    
    zones = []
        
    nb_zones = len(res)
    #print nb_zones
    if (nb_zones == 1) :
      zones.append(C.convertArrays2ZoneNode('okstar', [res[0]]))
    else:

      if (len(res[0][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('open_cells', [res[0]]))
      if (len(res[1][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('show_stopper', [res[1]]))
      if (len(res[2][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('w_uncomputable_pgs', [res[2]]))
      if (len(res[3][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('non_star', [res[3]]))
      if (len(res[4][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('neighbors', [res[4]]))
      if (len(res[5][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('good', [res[5]]))

    return zones

#==============================================================================
# extractOuters : XXX
#==============================================================================
def extractOuterLayers(mesh):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    res = XOR.extractOuterLayers(m)
    
    zones = []
        
    nb_zones = len(res)
    #print nb_zones
    if (nb_zones == 1) :
      zones.append(C.convertArrays2ZoneNode('remaining', [res[0]]))
    else:
      zones.append(C.convertArrays2ZoneNode('outers', [res[0]]))
      zones.append(C.convertArrays2ZoneNode('remaining', [res[1]]))

    return zones

#==============================================================================
# extractNthCell : XXX
#==============================================================================
def extractNthCell(mesh, nth):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.extractNthCell(m, nth)
    return C.convertArrays2ZoneNode('cell_%d'%(nth), [m])

#==============================================================================
# extractNthFace : XXX
#==============================================================================
def extractNthFace(mesh, nth):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.extractNthFace(m, nth)
    return C.convertArrays2ZoneNode('face_%d'%(nth), [m])

#==============================================================================
# removeNthCell : XXX
#==============================================================================
def removeNthCell(mesh, nth):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.removeNthCell(m, nth)
    return C.convertArrays2ZoneNode('mes_wo_%d'%(nth), [m])

#==============================================================================
# diffMesh : returns the diff between 2 meshes as 2 zones
#==============================================================================
def diffMesh(mesh1, mesh2):
    m1 = C.getFields(Internal.__GridCoordinates__, mesh1)[0]
    m2 = C.getFields(Internal.__GridCoordinates__, mesh2)[0]
    
    res = XOR.diffMesh(m1, m2)
    
    zones = []
    nb_zones = len(res)

    if (nb_zones == 0) : 
        print "No difference." ; return zones
    
    zones.append(C.convertArrays2ZoneNode('z1', [res[0]]))
    zones.append(C.convertArrays2ZoneNode('z2', [res[1]]))
    
    return zones

#==============================================================================
# selfX : returns any intersecting cells in the input mesh
#==============================================================================
def selfX(mesh):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.selfX(m)
    return C.convertArrays2ZoneNode('selfX', [m])

#==============================================================================
# checkCellsClosure : XXX
#==============================================================================
def checkCellsClosure(mesh):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    return XOR.checkCellsClosure(m)

def extrudeUserDefinedBC(mesh):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    cur_shift=0
    extrudepgs=[]
    zones = Internal.getZones(mesh)
    print "nb of zones %d"%(len(zones))
    (extrudepgs, cur_shift) = concatenateBC('BCUserDefined', [zones], extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    print "nb of pgs to pass : %s" %(len(extrudepgs))

    mo = XOR.extrudeUserDefinedBC(m, extrudepgs)

    return C.convertArrays2ZoneNode('union', [mo])

def reorientExternalFaces(mesh):
    m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
    m = XOR.reorientExternalFaces(m)
    return C.convertArrays2ZoneNode('oriented', [m])

    
#~ def conservativeTransfer(a1, a2, tol=0., reconstruction_type=0):
    #~ 
    #~ s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    #~ flowsol1 = C.getFields(Internal.__FlowSolutionCenters__, a1)[0]
    #~ s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
#~ 
    #~ flowsol2 = XOR.conservativeTransfer(s1, flowsol1, s2, tol, reconstruction_type)
    #~ 
    #~ C._deleteFlowSolutions__(a2)
    #~ return C.setFields([flowsol2], a2, 'centers')
    #~ 
#~ def totalMass(a1):
#~ 
    #~ s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    #~ flowsol1 = C.getFields(Internal.__FlowSolutionCenters__, a1)[0]
    #~ 
    #~ return XOR.totalMass(s1, flowsol1)
    #~ 
#~ 
#~ def normL1(a1,a2):
    #~ 
    #~ s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    #~ flowsol1 = C.getFields(Internal.__FlowSolutionCenters__, a1)[0]
    #~ 
    #~ s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    #~ flowsol2= C.getFields(Internal.__FlowSolutionCenters__, a2)[0]
#~ 
    #~ return generator.normL1(s1, flowsol1, s2, flowsol2)
#~ 
