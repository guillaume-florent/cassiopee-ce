"""Grid generation module.
"""
# 
# Python Interface to create PyTrees defining meshes
#
import Generator
import generator
__version__ = Generator.__version__

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
except:
    raise ImportError("Generator.PyTree: requires Converter.PyTree module.")

def cart(Xo, H, N):
    """Create a structured cartesian mesh.
    Usage: cart((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    a = generator.cart(Xo, H, N, 2)
    return C.convertArrays2ZoneNode('cart', [a])

def cartHexa(Xo, H, N):
    """Create a hexahedral cartesian mesh.
    Usage: cartHexa((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    a = generator.cartHexa(Xo, H, N, 2)
    return C.convertArrays2ZoneNode('cartHexa', [a])

def cartTetra(Xo, H, N):
    """Create a tetrahedrical cartesian mesh.
    Usage: cartTetra((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    a = generator.cartTetra(Xo, H, N, 2)
    return C.convertArrays2ZoneNode('cartTetra', [a])

def cartPenta(Xo, H, N):
    """Create a prismatic cartesian mesh.
    Usage: cartPenta((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    a = generator.cartPenta(Xo, H, N, 2)
    return C.convertArrays2ZoneNode('cartPenta', [a])

def cartPyra(Xo, H, N):
    """Create a pyramidal cartesian mesh.
    Usage: cartPyra((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    a = generator.cartPyra(Xo, H, N, 2)
    return C.convertArrays2ZoneNode('cartPyra', [a])

def cartNGon(Xo, H, N):
    """Create a NGON cartesian mesh.
    Usage: cartNGon((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    a = generator.cartNGon(Xo, H, N, 2)
    return C.convertArrays2ZoneNode('cartNGon', [a])

#------------------------------------------------------------------------------
# Generation d'un quadtree en 2D ou octree en 3D a partir d'une liste
# de contours ou surfaces
#------------------------------------------------------------------------------
def octree(surfaces, snearList, dfar=5., balancing=0, levelMax=1000, ratio=2):
    """Generate an octree (or a quadtree) mesh starting from a list of TRI
    (or BAR) arrays defining bodies, a list of corresponding snears,
    and the extension dfar of the mesh.
    Usage: octree(surfaces, snearList, dfar, balancing, levelMax, ratio)"""
    stlArrays = C.getFields(Internal.__GridCoordinates__, surfaces)
    stlArrays = Converter.convertArray2Tetra(stlArrays)
    a = Generator.octree(stlArrays, snearList, dfar, balancing, levelMax,ratio)
    return C.convertArrays2ZoneNode('octree', [a])

#==============================================================================
def conformOctree3(o):
    """Conformize an octree3.
    Usage: conformOctree3(octree)"""
    a = C.getAllFields(o, 'nodes')[0]
    a2 = Generator.conformOctree3(a)
    return C.convertArrays2ZoneNode(o[0], [a2])
    
#------------------------------------------------------------------------------
# Conformisation d'une soupe de TRI ou de BAR
#------------------------------------------------------------------------------
def conformUnstr(surface1, surface2=None, tol=0., left_or_right=0):
    """Conformizes a TRI or BAR soup (surface1) with taking into account surface2 if it's provided.
    Usage: conformUnstr(s1, s2, tol)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    if surface2 is not None:
        s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    else: s2 = None
    s = Generator.conformUnstr(s1, s2, tol, left_or_right)
    return C.convertArrays2ZoneNode('conformized', [s])

#------------------------------------------------------------------------------
# Conversion d un maillage octree en ensemble de grilles cartesiennes
#------------------------------------------------------------------------------
def octree2Struct(o, vmin=15, ext=0, optimized=1, merged=1, AMR=0,
                  sizeMax=1000000000):
    """Generates a structured set of regular Cartesian grid starting from an
    octree HEXA or QUAD mesh. vmin is the number of minimum points per grid,
    and can be specified for each level;
    ext is the extension of grids in all the directions. If optimized=1,
    then the extension can be reduced for minimum overlapping.
    merged=1 means that Cartesian grids are merged when possible.
    If AMR = 1, a list of AMR grids is generated.
    Usage: octree2Struct(a, vmin, ext, optimized, merged, AMR)"""
    if ext == 1: ext = 2
    dim = Internal.getZoneDim(o)
    if (dim[0] != 'Unstructured'): raise ValueError("octree2Struct: zone must be unstructured.")
    if (dim[3] == 'QUAD'): dimPb = 2
    elif (dim[3] == 'HEXA'): dimPb = 3
    else: raise ValueError("octree2Struct: zone must be QUAD or HEXA.")

    eps = 1.e-6
    a = C.getFields(Internal.__GridCoordinates__, o)[0]

    # Conversion en structure
    cartzones = Generator.octree2Struct(a, vmin, ext, optimized, merged,
                                        AMR, sizeMax)
    
    # Creation des zones du pyTree
    c = 1; zones = []
    for mc in cartzones:
        zone = C.convertArrays2ZoneNode('cart'+str(c), [mc])        
        zones.append(zone); c += 1
    if (AMR == 1): return zones
    if (ext == 0):
        try: import Connector.PyTree as X
        except: 
            print 'Warning: octree2Struct requires Connector module. No grid connectivity built.'
            return zones
        
        if dimPb == 3: ratios = [[2,2,2],[4,4,4],[8,8,8],[16,16,16]]
        else: ratios = [[2,2,1],[4,4,1],[8,8,1],[16,16,1]]
        zones = X.connectMatch(zones, dim=dimPb)
        for ratio0 in ratios:
            zones = X.connectNearMatch(zones,ratio=ratio0,dim=dimPb)
        return zones 
    else:
        #-----------------------
        # Creation des BCOverlap
        #-----------------------
        # determination de la bounding box de la grille
        bbox0 = bbox(zones)
        xmin = bbox0[0]; ymin = bbox0[1]; zmin = bbox0[2]
        xmax = bbox0[3]; ymax = bbox0[4]; zmax = bbox0[5]
        noz = 0
        for z in zones:
            [x1,y1,z1,x2,y2,z2] = bbox(z)
            if (x1 > xmin+eps): z=C.addBC2Zone(z,'overlap1','BCOverlap','imin')
            if (x2 < xmax-eps): z=C.addBC2Zone(z,'overlap2','BCOverlap','imax')
            if (y1 > ymin+eps): z=C.addBC2Zone(z,'overlap3','BCOverlap','jmin')
            if (y2 < ymax-eps): z=C.addBC2Zone(z,'overlap4','BCOverlap','jmax')
            if (z1 > zmin+eps): z=C.addBC2Zone(z,'overlap5','BCOverlap','kmin')
            if (z2 < zmax-eps): z=C.addBC2Zone(z,'overlap6','BCOverlap','kmax')
            zones[noz] = z; noz += 1
    return zones

def _adaptOctree(a,indicator="indicator",balancing=1,ratio=2):
    indicator = indicator.split(':')
    if len(indicator) == 2: indicator = indicator[1]
    else: indicator = indicator[0]
    zones = Internal.getZones(a)
    hexa = C.getFields(Internal.__GridCoordinates__,zones)
    indic = C.getFields(Internal.__FlowSolutionCenters__,zones)
    indic = Converter.extractVars(indic, [indicator])
    C._deleteFlowSolutions__(a)
    for noz in xrange(len(zones)):        
        res = Generator.adaptOctree(hexa[noz], indic[noz], balancing, ratio)
        C.setFields([res],zones[noz],'nodes',writeDim=True)
    return None

def adaptOctree(a, indicator='indicator', balancing=1, ratio=2):
    """Adapt the octree with respect to the field 'indicator' located at centers.
    Usage: adaptOctree(a, indicator, balancing, ratio)"""
    tp = Internal.copyRef(a)
    _adaptOctree(tp,indicator=indicator,balancing=balancing,ratio=ratio)
    return tp

def expandLayer(o, level=0, corners=0, balancing=0):
    """Expand layer of level l of an unstructured quadtree/octree.
    Usage: expandLayer(o,level,corners,balancing)"""
    tp = Internal.copyRef(o)
    _expandLayer(tp,level=level,corners=corners,balancing=balancing)
    return tp

def _expandLayer(o, level=0, corners=0, balancing=0):
    zones = Internal.getZones(o)
    hexa = C.getFields(Internal.__GridCoordinates__,zones)
    for noz in xrange(len(zones)):
        res = Generator.expandLayer(hexa[noz], level, corners, balancing)
        C.setFields([res],zones[noz],'nodes',writeDim=True)
    return None

def cylinder(Xo, R1, R2, tetas, tetae, H, N):
    """Create a portion of regular cylindrical grid.
    Usage: cylinder((xo,yo,zo), R1, R2, tetas, tetae, H, (ni,nj,nk))"""
    a = Generator.cylinder(Xo, R1, R2, tetas, tetae, H, N)
    return C.convertArrays2ZoneNode('cylinder', [a])

def cylinder2(Xo, R1, R2, tetas, tetae, H, dR, dteta, dZ):
    """Create a portion of regular cylindrical grid.
    Usage: cylinder2((xo,yo,zo), R1, R2, tetas, tetae, H, dR, dteta, dZ)"""
    ar = C.getFields(Internal.__GridCoordinates__, dR)[0]
    at = C.getFields(Internal.__GridCoordinates__, dteta)[0]
    az = C.getFields(Internal.__GridCoordinates__, dZ)[0]
    a = Generator.cylinder2(Xo, R1, R2, tetas, tetae,
                            H, ar, at, az)
    return C.convertArrays2ZoneNode('cylinder', [a])

def cylinder3(distxz, tetas, tetae, dteta):
    """Create a portion of cylindrical grid.
    Usage: cylinder3(arrayxz, tetas, tetae, dteta)"""
    axz = C.getFields(Internal.__GridCoordinates__, distxz)[0]
    at = C.getFields(Internal.__GridCoordinates__, dteta)[0]
    a = Generator.cylinder3(axz, tetas, tetae, at)
    return C.convertArrays2ZoneNode('cylinder', [a])

def delaunay(a, tol=1.e-10, keepBB=0):
    """Create a delaunay mesh given a set of points defined by array.
    Usage: delaunay(mesh, tol, keepBB)"""
    contour = C.deleteFlowSolutions__(a, 'centers')
    return C.TZGC(contour, 'nodes', Generator.delaunay, tol, keepBB)

def checkDelaunay(contour, tri):
    """Check if the Delaunay triangulation defined by tri is inside the contour.
    Usage: checkDelaunay(contour, tri)"""
    contour = C.deleteFlowSolutions__(contour, 'centers')
    tria = C.getFields(Internal.__GridCoordinates__, tri)[0]
    return C.TZA(contour, 'nodes', 'nodes',
                 Generator.checkDelaunay, None, tria)

def constrainedDelaunay(contour, tol=1.e-10, keepBB=0):
    """Create a constrained-Delaunay mesh starting from a BAR-array defining
    the contour.
    Usage: constrainedDelaunay( contour, tol, keepBB )"""
    contour = C.deleteFlowSolutions__(contour, 'centers')
    return C.TZGC(contour, 'nodes',
                  Generator.constrainedDelaunay, tol, keepBB)

def plaster(contours, surfaces, side=0):
    """Create a sticky plaster around contours surrounded by surfaces.
    Usage: plaster(contours, surfaces)"""
    surfn = C.getAllFields(surfaces, 'nodes')
    cs = C.getAllFields(contours, 'nodes')
    pa = Generator.plaster(cs, surfn, side)
    return C.convertArrays2ZoneNode('plaster', [pa])

def fittingPlaster(contour, bumpFactor=0.):
    """Create a sticky plaster around a contour and pump it.
    Usage: plaster(contour, bumpFactor)"""
    c = C.getAllFields(contour, 'nodes')[0]
    pa = Generator.fittingPlaster(c, bumpFactor)
    return C.convertArrays2ZoneNode('plaster', [pa])

def T3mesher2D(a, triangulateOnly=0):
    """Create a delaunay mesh given a set of points defined by a.
    Usage: T3mesher2D(a, triangulateOnly) """
    c = C.getAllFields(a, 'nodes')[0]
    c = Generator.T3mesher2D(c, triangulateOnly)
    return C.convertArrays2ZoneNode('tri', [c]) 

def tetraMesher(a, maxh=-1., grading=0.4, triangulateOnly=0, 
                remeshBoundaries=0, algo=1):
    """Create a TRI/TETRA mesh given a set of BAR or surfaces in a.
    Usage: tetraMesher(a, fineness, grading)"""
    c = C.getAllFields(a, 'nodes')
    c = Generator.tetraMesher(c, maxh, grading, triangulateOnly,
                              remeshBoundaries, algo)
    return C.convertArrays2ZoneNode('mesh', [c])

def gapfixer(contour, cloud, hardPoints=None, refine=1):
    """Fix a gap defined by a contour bar and a point cloud representing the gap surface.
    Some hard points can be sepcified to force the constructed surface to pass by.
    If the optional refine argument is set to 0, the resulting surface will be a contrained triangulation of the contour [and the additional hard nodes].
    Usage: gapFixer(contour, cloud, hardPoints, refine)"""
    clouda = C.getFields(Internal.__GridCoordinates__, cloud)[0]
    c = C.getAllFields(contour, 'nodes')[0]
    hp = None
    if hardPoints is not None:
        hp = C.getFields(Internal.__GridCoordinates__, hardPoints)[0]
    c = Generator.gapfixer(c, clouda, hp, refine)
    return C.convertArrays2ZoneNode('tri', [c])

def gapsmanager(components, mode=0, refine=0, coplanar=0):
    """Fix a gap between several component surfaces (list of arrays).
    The mode sets the case (0 for POST with a center mesh, 1 for POST with a nodal mesh, 2 otherwise).
    The planar argument tells whether the components are coplanar or not (1 for coplanar case).
    Usage: gapsmanager(components, mode, refine, coplanar)"""
    stlArrays = C.getAllFields(components, loc='nodes')
    stlArrays = Converter.convertArray2Tetra(stlArrays)
    a = Generator.gapsmanager(stlArrays, mode, refine, coplanar)
    zones = []
    for z in a:
        zones.append(C.convertArrays2ZoneNode('gap', [z]))
    return zones

def front2Hexa(a, surf, h, hf, hext, density=50):
    """Generates an hexa grid starting from a front a, a surface surf,
    and h, hf, hext the height of the mesh, of the first and last cells."""
    arrayFront = C.getFields(Internal.__GridCoordinates__, a)[0]
    arraySurf = C.getFields(Internal.__GridCoordinates__, surf)[0]
    zone = Generator.front2Hexa(arrayFront, arraySurf, h, hf, hext, density)
    zone = C.convertArrays2ZoneNode('fill', [zone])
    return zone

def front2Struct(front, surf, distrib, Vmin):
    """Generates struct grids starting from a front, a surface surf,
    and a point distribution."""
    arrayFront = C.getFields(Internal.__GridCoordinates__, front)[0]
    arraySurf = C.getFields(Internal.__GridCoordinates__, surf)[0]
    arrayDistrib = C.getFields(Internal.__GridCoordinates__, distrib)[0]
    zones = Generator.front2Struct(arrayFront, arraySurf, arrayDistrib, Vmin)
    for z in zones:
        z = C.addBC2Zone(z, 'wall', 'BCWall', 'kmin')
        z = C.addBC2Zone(z, 'overlap', 'BCOverlap', 'kmax')
    return zones

#------------------------------------------------------------------------------
# Projection de t sur des surfaces surfs
# optimized=0,1,2 (front exterieur, front optimise, front avec projection sur les contours)
# Les arguments suivants sont requis si optimized=2 (ignores sinon):
# step: pas de discretisation pour un remaillage regulier des surfaces
# (ce pas doit correspondre a peu pres a la taille des mailles de t projetees)
# angle: angle suivant lequel les surfaces sont decoupees
#------------------------------------------------------------------------------
def snapFront(t, surfs, optimized=1):
    """Adapt t to a given surface (cellN defined in t). 
    Usage: snapFront(t, surfs, step, angle, optimized)"""
    arrays = C.getFields(Internal.__GridCoordinates__, surfs)
    return C.TZA(t, 'nodes', 'nodes', Generator.snapFront, None,
                 arrays, optimized)

def _snapFront(t, surfs, optimized=1):
    arrays = C.getFields(Internal.__GridCoordinates__, surfs)
    return C._TZA(t, 'nodes', 'nodes', Generator.snapFront, None,
                  arrays, optimized)

#------------------------------------------------------------------------------
# Deplacement de points de t sur ceux des surfaces surfs discretisees
#------------------------------------------------------------------------------
def snapSharpEdges(t, surfs, step=None, angle=30.):
    """Adapt t to a given surface. 
    Usage: snapSharpEdges(t, surfs, step)"""
    arrays = C.getFields(Internal.__GridCoordinates__, surfs)
    return C.TZA(t, 'nodes', 'nodes', Generator.snapSharpEdges, None, 
                 arrays, step, angle)
    
def _snapSharpEdges(t, surfs, step=None, angle=30.):
    arrays = C.getFields(Internal.__GridCoordinates__, surfs)
    return C._TZA(t, 'nodes', 'nodes', Generator.snapSharpEdges, None, 
                  arrays, step, angle)

def check(t):
    """Check a mesh for regularity, orthogonality...
    Usage: check( t )"""
    a = C.getFields(Internal.__GridCoordinates__, t)
    for i in a: Generator.check(i)

def bbox(t):
    """Returns the bounding box of a pytree.
    Usage: bbox(t)"""
    A = C.getFields(Internal.__GridCoordinates__, t)
    return Generator.bbox(A)

def BB(t, method='AABB', weighting=0):
    """Return the bounding box of a pyTree as a pyTree.
    Usage: b = BB(a, method, weighting)"""
    tp = Internal.copyRef(t)
    C._deleteFlowSolutions__(tp)
    C._deleteZoneBC__(tp)
    C._deleteGridConnectivity__(tp)
    C._TZGC(tp, 'nodes', Generator.BB, method, weighting)
    return tp

def _BB(t, method='AABB', weighting=0):
    """Return the bounding box of a pyTree as a pyTree.
    Usage: b = BB(a, method, weighting)"""
    C._deleteFlowSolutions__(t)
    C._deleteZoneBC__(t)
    C._deleteGridConnectivity__(t)
    C._TZGC(t, 'nodes', Generator.BB, method, weighting)
    return None

def barycenter(t, weight='None'):
    """Get the barycenter of a pyTree.
    Usage: barycenter(t)"""
    A = C.getFields(Internal.__GridCoordinates__, t)
    if (weight != 'None'):
        W = C.getField(weight, t)
        return Generator.barycenter(A, W)
    else: return Generator.barycenter(A)

def CEBBIntersection(a1, a2, tol=1.e-10):
    """Get the Cartesian Elements bounding box intersection of 2 meshes."""
    a1 = Internal.getZones(a1)
    a2 = Internal.getZones(a2)
    if (len(a1) != 1 or len(a2) != 1):
        print 'Warning: CEBBIntersection applied on one zone.'
        return 0
    
    m1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    m2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    return Generator.CEBBIntersection(m1, m2, tol)

        
def bboxIntersection(z1, z2, tol=1.e-6, isBB=False, method='AABB'):
    """Return 1 if bounding boxes of z1 and z2 intersect."""       
    if not isBB:
        z1 = BB(z1,method)
        z2 = BB(z2,method)
    if method == 'AABB':  # Computes the intersection between 2 AABB
        return generator._bboxIntersectionZ(z1, z2, tol, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)
    elif method == 'OBB':  # Computes the intersection between 2 OBB
        return generator._obboxIntersectionZ(z1, z2, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)
    elif method == 'AABBOBB': # Intersection between AABB and OBB
        return generator._crossIntersectionZ(z1, z2, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)
    else:
        print 'Warning: bboxIntersection: method',method,'not implemented, switching to AABB.'
        return generator._bboxIntersectionZ(z1, z2, tol,
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)          

def _bboxIntersection(z1, z2, tol=1.e-6, isBB=False, method='AABB'):
    """Return 1 if bounding boxes of z1 and z2 intersect."""       
    if not isBB:
        _BB(z1,method)
        _BB(z2,method)
    if method == 'AABB':  # Computes the intersection between 2 AABB
        return generator._bboxIntersectionZ(z1, z2, tol, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)
    elif method == 'OBB':  # Computes the intersection between 2 OBB
        return generator._obboxIntersectionZ(z1, z2, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)
    elif method == 'AABBOBB': # Intersection between AABB and OBB
        return generator._crossIntersectionZ(z1, z2, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)
    else:
        print 'Warning: bboxIntersection: method',method,'not implemented, switching to AABB.'
        return generator._bboxIntersectionZ(z1, z2, tol, 
                                            Internal.__GridCoordinates__, 
                                            Internal.__FlowSolutionNodes__, 
                                            Internal.__FlowSolutionCenters__)          

def checkPointInCEBB(t, P):
    """Check if point P is in the Cartesian Elements Bounding Box of a mesh."""
    m = C.getFields(Internal.__GridCoordinates__,t)[0]
    return Generator.checkPointInCEBB(m, P)

def bboxOfCells(t):
    """Compute the bounding box of all cells of a mesh.
    Usage: getBBoxOfCells(t)"""
    return C.TZGC(t, 'centers', Generator.bboxOfCells)

def _bboxOfCells(t):
    return C._TZGC(t, 'centers', Generator.bboxOfCells)

def getVolumeMap(t):
    """Return the volume map in an array.
    Usage: getVolumeMap(t)"""
    return C.TZGC(t, 'centers', Generator.getVolumeMap)

def _getVolumeMap(t):
    """Return the volume map in an array.
    Usage: _getVolumeMap(t)"""
    return C._TZGC(t, 'centers', Generator.getVolumeMap)

def getNormalMap(t):
    """Return the map of surface normals in an array.
    Usage: getNormalMap(t)""" 
    return C.TZGC(t, 'centers', Generator.getNormalMap)

def _getNormalMap(t):
    return C._TZGC(t, 'centers', Generator.getNormalMap)

def getSmoothNormalMap(t, niter=2, eps = 0.4):
    """Return the map of smoothed and non-normalized surface normals in an array.
    eps is the smoothing factor.
    Usage: getSmoothNormalMap(t, niter, eps)"""
    return C.TZGC(t, 'nodes', Generator.getSmoothNormalMap, niter, eps)

def _getSmoothNormalMap(t, niter=2, eps = 0.4):
    return C._TZGC(t, 'nodes', Generator.getSmoothNormalMap, niter, eps)

def getCellPlanarity(t):
    """Return the cell planarity of a surface mesh in an array.
    Usage: getCellPlanarity(t)"""
    return C.TZGC(t, 'centers', Generator.getCellPlanarity)

def _getCellPlanarity(t):
    return C._TZGC(t, 'centers', Generator.getCellPlanarity)

def getCircumCircleMap(t):
    """Return the map of circum circle radius of a 'TRI' array.
    Usage: getCircumCircleMap(t)""" 
    return C.TZGC(t, 'centers', Generator.getCircumCircleMap)

def _getCircumCircleMap(t):
    return C._TZGC(t, 'centers', Generator.getCircumCircleMap)

def getInCircleMap(t):
    """Return the map of inscribed circle radius of a 'TRI' array.
    Usage: getInCircleMap(t)""" 
    return C.TZGC(t, 'centers', Generator.getInCircleMap)

def _getInCircleMap(t):
    return C._TZGC(t, 'centers', Generator.getInCircleMap)

def getEdgeRatio(t):
    """Computes the ratio between the max and min lengths of all the edges of
    cells in an array.
    Usage: getEdgeRatio(t)"""
    return C.TZGC(t,'centers', Generator.getEdgeRatio)

def _getEdgeRatio(t):
    return C._TZGC(t,'centers', Generator.getEdgeRatio)


def getMaxLength(t):
    """Computes the max length of all the edges of cells in a zone.
    Usage: getMaxLength(t)"""
    return C.TZGC(t,'centers', Generator.getMaxLength)

def _getMaxLength(t):
    return C._TZGC(t,'centers', Generator.getMaxLength)

def enforceX(a, x0, enforcedh, N, add=0):
    """Enforce a x0-centered line in a distribution defined by an array.
    Usage: enforceX(a, x0, enforcedh, supp, add) -or-
    Usage: enforceX(a, x0, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforceX, x0, enforcedh, N, add )
    dir = 1
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforceY(a, y0, enforcedh, N, add=0):
    """Enforce a y0-centered line in a distribution defined by an array.
    Usage: enforceY(a, y0, enforcedh, supp, add) -or-
    Usage: enforceY(a, y0, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforceY, y0, enforcedh, N, add )
    dir = 2
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforceZ(a, z0, enforcedh, N, add=0):
    """Enforce a y0-centered line in a distribution defined by an array.
    Usage: enforceZ(a, z0, enforcedh, supp, add) -or-
    Usage: enforceZ(a, z0, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforceZ, z0, enforcedh, N, add )
    dir = 3
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforcePlusX(a, enforcedh, N, add=0):
    """Enforce the first X-line in a distribution defined by an array.
    (one sided distribution, right).
    Usage: enforcePlusX(array, enforcedh, supp, add) -or-
    Usage: enforcePlusX(array, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforcePlusX, enforcedh, N, add )
    dir = 1
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforcePlusY(a, enforcedh, N, add=0):
    """Enforce the first Y-line in a distribution defined by an array.
    (one sided distribution, right).
    Usage: enforcePlusY(array, enforcedh, supp, add) -or-
    Usage: enforcePlusY(array, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforcePlusY, enforcedh, N, add )
    dir = 2
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforcePlusZ(a, enforcedh, N, add=0):
    """Enforce the first Z-line in a distribution defined by an array.
    (one sided distribution, right).
    Usage: enforcePlusZ(array, enforcedh, supp, add) -or-
    Usage: enforcePlusZ(array, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforcePlusZ, enforcedh, N, add )
    dir = 3
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforceMoinsX(a, enforcedh, N, add=0):
    """Enforce the last X-line in a distribution (one sided, left).
    Usage: enforceMoinsX(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsX(array, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforceMoinsX, enforcedh, N, add )
    dir = 1
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforceMoinsY(a, enforcedh, N, add=0):
    """Enforce the last Y-line in a distribution (one sided, left).
    Usage: enforceMoinsY(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsY(array, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforceMoinsY, enforcedh, N, add )
    dir = 2
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforceMoinsZ(a, enforcedh, N, add=0):
    """Enforce the last Z-line in a distribution (one sided, left).
    Usage: enforceMoinsZ(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsZ(array, enforcedh, (supp,add))"""
    dims = Internal.getZoneDim(a)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    a = C.deleteFlowSolutions__(a)
    a = C.TZGC(a, 'nodes', Generator.enforceMoinsZ, enforcedh, N, add)
    dir = 3
    a = modifyBC__(dir, ni0, nj0, nk0, a)
    return a

def enforceLine(a, line, enforcedh, N):
    """Enforce a line in a distribution.
    Usage: enforceLine(array, line, enforcedh, (supp,add))"""
    a = C.deleteAllBCAndSolutions__(a)
    l = C.getFields(Internal.__GridCoordinates__,line)[0]
    return C.TZGC(a,  'nodes', Generator.enforceLine, l, enforcedh, N)

def enforcePoint(a, x0):
    """Enforce a point in a distribution.
    Usage: enforcePoint(a, x0)"""
    a = C.deleteAllBCAndSolutions__(a)
    return C.TZGC(a, 'nodes', Generator.enforcePoint, x0)

def enforceCurvature(distrib, curvature, power=0.5):
    """Enforce curvature of a curve in a distribution.
    The distribution is defined by distrib, and the curvature by curvature
    Usage: enforceCurvature( distrib, curvature, power )"""
    distrib = C.deleteAllBCAndSolutions__(distrib)
    ac = C.getFields(Internal.__GridCoordinates__, curvature)[0]
    return C.TZGC(distrib, 'nodes', Generator.enforceCurvature, ac, power)

def enforceCurvature2(distrib, curve, alpha=1.e-2):
    """Enforce a 1D distribution wrt the curvature radius
    alpha is the factor of stretch for point of maximum curvature.
    Usage: enforceCurvature2(distrib, curve, alpha)"""
    distrib = C.deleteAllBCAndSolutions__(distrib)
    ac = C.getFields(Internal.__GridCoordinates__, curve)[0]
    return C.TZGC(distrib, 'nodes', Generator.enforceCurvature2, ac, alpha)

def addPointInDistribution(a, ind):
    """Add a point in a distribution defined by a.
    Usage: addPointInDistribution(a, ind)"""
    a = C.deleteAllBCAndSolutions__(a)
    return C.TZGC(a, 'nodes', Generator.addPointInDistribution, ind)

def close(a, tol=1.e-12):
    """Close a mesh defined by an array gathering points nearer than tol.
    Usage: close(array, tol)"""
    t = Internal.copyRef(a)
    fields = C.getAllFields(t, 'nodes')
    fields = Generator.close(fields, tol)
    C.setFields(fields, t, 'nodes')
    return t

def _close(t, tol=1.e-12):
    fields = C.getAllFields(t, 'nodes')
    fields = Generator.close(fields, tol)
    C.setFields(fields, t, 'nodes')
    return None

def pointedHat(a, (x,y,z)):
    """Create a structured surface defined by a contour and a point (x,y,z).
    Usage: pointedHat(array, (x,y,z))"""
    a = C.deleteFlowSolutions__(a, 'both')
    return C.TZA(a, 'nodes', 'nodes', Generator.pointedHat, None, (x,y,z))

def stitchedHat(a, (offx,offy,offz), tol=1.e-6, tol2=1.e-5):
    """Create a structured surface defined by a contour and an offset (dx,dy,dz).
    Usage: stitchedHat(a, (dx,dy,dz))"""
    a = C.deleteFlowSolutions__(a, 'both')
    return C.TZGC(a, 'nodes', Generator.stitchedHat, \
                  (offx,offy,offz), tol, tol2)

def selectInsideElts(a, curvesList):
    """Select elements whose center is in the surface delimited by curves.
    Usage: selectInsideElts(array, curvesList)"""
    a = C.deleteFlowSolutions__(a, 'centers')
    curves = C.getFields(Internal.__GridCoordinates__, curvesList)
    return C.TZA(a, 'nodes', 'nodes', Generator.selectInsideElts, None,
                 curves)

def grow(t, vector):
    """Grow a surface array of one layer by moving points of vector.
    Usage: grow( t, vector )"""
    tp = Internal.copyRef(t)
    if (len(vector) != 3): raise ValueError("grow: 3 variables are required.")
    nodes = Internal.getZones(tp)
    for z in nodes:
        fa = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if (fa != []):
            a = Converter.extractVars(fa, vector)
        else:
            print "Warning: grow: variables not found in zone."
            a = []
        if (a != []):
            nodes = C.getAllFields(z, 'nodes')[0]
            nodes = Generator.grow(nodes, a)
            C.setFields([nodes], z, 'nodes')
    tp = Internal.addOneLayer2BC(tp, 3)
    return tp

def stack(t1, t2):
    """Stack two meshes (with same nixnj) into a single mesh.
    Usage: stack(a1, a2)"""
    a2 = C.getAllFields(t2, 'nodes')[0]
    return C.TZA(t1, 'nodes', 'nodes', Generator.stack, None, a2)

#================================================================
# traitement des BC seulement definies comme fenetres completes
# dir: direction dans laquelle le maillage a ete modifie
# ni0,nj0,nk0: dimension du maillage avant modification
# z: zone apres modification.
#================================================================
def modifyBC__(dir, ni0, nj0, nk0, z):
    z = C.rmBCOfType(z, 'BCMatch'); z = C.rmBCOfType(z, 'BCNearMatch')
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured': return z
    ni = dims[1]; nj = dims[2]; nk = dims[3]
    wins = Internal.getNodesFromType(z, 'BC_t')
    # BC
    for w in wins:
        (parent, d) = Internal.getParentOfNode(z, w)
        w0 = w[2][0][1]
        i1 = w0[0,0]; j1 = w0[1,0]; k1 = w0[2,0]
        i2 = w0[0,1]; j2 = w0[1,1]; k2 = w0[2,1] 
        if dir == 1:
            if i1 == 1 and i2 == ni0:
                range0 = [1,ni,j1,j2,k1,k2]
                r2 = Internal.window2Range(range0)
                parent[2][d][2][0][1] = r2
            elif i1 == i2 and i1 == ni0:
                range0 = [ni,ni,j1,j2,k1,k2]
                r2 = Internal.window2Range(range0)
                parent[2][d][2][0][1] = r2
            elif i1 != i2: del parent[2][d]
        elif dir == 2:
            if j1 == 1 and j2 == nj0:
                range0 = [i1,i2,1,nj,k1,k2]
                r2 = Internal.window2Range(range0)
                parent[2][d][2][0][1] = r2
            elif j1 == j2 and j1 == nj0:
                range0 = [i1,i2,nj,nj,k1,k2]
                r2 = Internal.window2Range(range0)
                parent[2][d][2][0][1] = r2
            elif j1 != j2: del parent[2][d]
        else:
            if k1 == 1 and k2 == nk0:
                range0 = [i1,i2,j1,j2,1,nk]
                r2 = Internal.window2Range(range0)
                parent[2][d][2][0][1] = r2
            elif k1 == k2 and k1 == nk0:
                range0 = [i1,i2,j1,j2,nk,nk]
                r2 = Internal.window2Range(range0)
                parent[2][d][2][0][1] = r2
            elif k1 != k2: del parent[2][d]
    # Grid Connectivity
    connect = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    for cn in connect:
        wins = Internal.getNodesFromName2(cn, 'PointRange')
        for w in wins:
            (parent, d) = Internal.getParentOfNode(cn, w)
            w0 = w[1]
            i1 = w0[0,0]; j1 = w0[1,0]; k1 = w0[2,0]
            i2 = w0[0,1]; j2 = w0[1,1]; k2 = w0[2,1] 
            if dir == 1:
                if i1 == 1 and i2 == ni0:
                    range0 = [1,ni,j1,j2,k1,k2]
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
                elif i1 == i2 and i1 == ni0:
                    range0 = [ni,ni,j1,j2,k1,k2]
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
                elif i1 != i2: z=C.rmNodes(z, parent[0])
            elif dir == 2:
                if j1 == 1 and j2 == nj0:
                    range0 = [i1,i2,1,nj,k1,k2]
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
                elif j1 == j2 and j1 == nj0:
                    range0 = [i1,i2,nj,nj,k1,k2]
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
                elif j1 != j2: z=C.rmNodes(z, parent[0])
            else:
                if k1 == 1 and k2 == nk0:
                    range0 = [i1,i2,j1,j2,1,nk]
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
                elif k1 == k2 and k1 == nk0:
                    range0 = [i1,i2,j1,j2,nk,nk]
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
                elif k1 != k2: z=C.rmNodes(z, parent[0])
    return z

def map(z, d, dir=0):
    """Map a distribution d on a curve defined by zone z.
    Usage: map( z, d )"""
    dims = Internal.getZoneDim(z)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    z = C.deleteFlowSolutions__(z)
    dist = C.getFields(Internal.__GridCoordinates__, d)[0]
    z = C.TZGC(z, 'nodes', Generator.map, dist, dir)
    z = modifyBC__(dir, ni0, nj0, nk0, z)
    return z

def mapCurvature(z, N, power, dir):
    """Remesh a mesh following curvature.
    Usage: mapCurvature(z, N, power, dir)"""
    dims = Internal.getZoneDim(z)
    ni0 = 1; nj0 = 1; nk0 = 1
    if dims[0] == 'Structured': ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
    z = C.deleteFlowSolutions__(z)
    z = C.TZGC(z, 'nodes', Generator.mapCurvature, N, power, dir)
    z = modifyBC__(dir, ni0, nj0, nk0, z)
    return z

def refineBCRanges__(r0, ni0, nj0, nk0, ni, nj, nk, dir, factor):
    if not isinstance(factor,int):
        raise ValueError("refineBCRanges__: factor must be an integer.")
    if factor == 1: return r0

    alp1 = 1; alp2 = 1; alp3 = 1
    if dir == 1 or dir == 0: alp1 = factor
    if dir == 2 or dir == 0: alp2 = factor
    if dir == 3 or dir == 0: alp3 = factor

    w0 = Internal.range2Window(r0)
    i1 = w0[0]; j1 = w0[2]; k1 = w0[4]
    i2 = w0[1]; j2 = w0[3]; k2 = w0[5]
    i1N = i1; i2N = i2; j1N = j1; j2N = j2; k1N = k1; k2N = k2
    shift = factor-1

    # nouveaux indices
    if dir == 1:            
        if i1 == 1: i1N = 1
        elif i1 == ni0: i1N = ni
        else: i1N = alp1*i1-shift
        if i2 == 1: i2N = 1
        elif i2 == ni0: i2N = ni
        else: i2N = alp1*i2-shift

    elif dir == 2:
        if j1 == 1: j1N = 1
        elif j1 == nj0: j1N = nj
        else: j1N = alp2*j1-shift
        if j2 == 1: j2N = 1
        elif j2 == nj0: j2N = nj
        else: j2N = alp2*j2-shift
            
    elif dir == 3:
        if k1 == 1: k1N = 1
        elif k1 == nk0: k1N = nk
        else: k1N = alp3*k1-shift
        if k2 == 1: k2N = 1
        elif k2 == nk0: k2N = nk
        else: k2N = alp3*k2-shift

    else: # dir = 0
        if i1 == 1: i1N = 1
        elif i1 == ni0: i1N = ni
        else: i1N = alp1*i1-shift
        if i2 == 1: i2N = 1
        elif i2 == ni0: i2N = ni
        else: i2N = alp1*i2-shift
        if j1 == 1: j1N = 1
        elif j1 == nj0: j1N = nj
        else: j1N = alp2*j1-shift
        if j2 == 1: j2N = 1
        elif j2 == nj0: j2N = nj
        else: j2N = alp2*j2-shift
        if k1 == 1: k1N = 1
        elif k1 == nk0: k1N = nk
        else: k1N = alp3*k1-shift
        if k2 == 1: k2N = 1
        elif k2 == nk0: k2N = nk
        else: k2N = alp3*k2-shift

    w0 = [i1N,i2N, j1N,j2N,k1N,k2N]
    return Internal.window2Range(w0)

def _refineBC__(z, ni0, nj0, nk0, factor, dir):
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured': return None
    ni = dims[1]; nj = dims[2]; nk = dims[3]
    wins = Internal.getNodesFromType2(z, 'BC_t')
    # BC
    for w in wins:
        r0 = Internal.getNodeFromName1(w, 'PointRange')
        r0[1] = refineBCRanges__(r0[1], ni0, nj0, nk0, ni, nj, nk, dir, factor)
    # Connectivite
    connect = Internal.getNodesFromType2(z, 'ZoneGridConnectivity_t')
    for cn in connect:
        prs = Internal.getNodesFromName2(cn, 'PointRange')
        for pr in prs:
            pr[1] = refineBCRanges__(pr[1], ni0, nj0, nk0, ni, nj, nk, dir, factor)
        prs = Internal.getNodesFromName3(cn, 'PointRangeDonor')
        for pr in prs:
            pr[1] = refineBCRanges__(pr[1], ni0, nj0, nk0, ni, nj, nk, dir, factor)
    return None

def refine(t, power, dir):
    tp = Internal.copyRef(t)
    _refine(tp, power, dir)
    return tp

def _refine(t, power, dir):
    """Refine a mesh keeping original point distribution.
    Usage: refine(z, power, dir)"""
    zones = Internal.getZones(t)
    for z in zones:
        dims = Internal.getZoneDim(z)
        if dims[0] != 'Structured':
            print "Warning: refine: zone must be structured."
        else:
            ni0 = dims[1]; nj0 = dims[2]; nk0 = dims[3]
            C._deleteFlowSolutions__(z)
            C._TZGC(z, 'nodes', Generator.refine, power, dir)
            # delete Chimera data : OversetHoles and InterpolationData
            C._deleteChimeraInfo__(z)
            # modify BCs
            factor = int(power)
            if float(factor)-power == 0.:
                if dir != 0: C._rmBCOfType(z,'BCMatch'); C._rmBCOfType(z,'BCNearMatch')
                else: _refineBC__(z, ni0, nj0, nk0, factor, dir)
            else: C._deleteZoneBC__(z); C._deleteGridConnectivity__(z)
    return None

def densify(z, h):
    """Return zone with densified mesh.
    Usage: densify(z, h)"""
    z = C.deleteFlowSolutions__(z)
    return C.TZA(z, 'nodes', 'nodes', Generator.densify, None, h)

def _densify(z, h):
    C._deleteFlowSolutions__(z)
    return C._TZA(z, 'nodes', 'nodes', Generator.densify, None, h)

def TFI(a):
    """Generate a transfinite interpolation mesh from boundaries.
    Usage: TFI(a)"""
    m = []
    for ai in a:
        m.append(C.getFields(Internal.__GridCoordinates__, ai)[0])
    r = Generator.TFI(m)
    return C.convertArrays2ZoneNode('tfi', [r])

def TFITri(a1, a2, a3):
    """Generate a transfinite interpolation mesh from 3 input curves.
    Usage: TFITri(a)"""
    a1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    a2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    a3 = C.getFields(Internal.__GridCoordinates__, a3)[0] 
    r1,r2,r3 = Generator.TFITri(a1, a2, a3)
    return [C.convertArrays2ZoneNode('tfi1', [r1]),
            C.convertArrays2ZoneNode('tfi2', [r2]),
            C.convertArrays2ZoneNode('tfi3', [r3])]

def TFIO(a):
    """Generate a transfinite interpolation mesh from 1 input curve.
    Usage: TFIO(a)"""
    a = C.getFields(Internal.__GridCoordinates__, a)[0]
    r1,r2,r3,r4,r5 = Generator.TFIO(a)
    return [C.convertArrays2ZoneNode('tfi1', [r1]),
            C.convertArrays2ZoneNode('tfi2', [r2]),
            C.convertArrays2ZoneNode('tfi3', [r3]),
            C.convertArrays2ZoneNode('tfi4', [r4]),
            C.convertArrays2ZoneNode('tfi5', [r5])]

def TFIHalfO(a1, a2):
    """Generate a transfinite interpolation mesh from 2 input curves.
    Usage: TFIHalfO(a)"""
    a1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    a2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    r1,r2,r3,r4 = Generator.TFIHalfO(a1, a2)
    return [C.convertArrays2ZoneNode('tfi1', [r1]),
            C.convertArrays2ZoneNode('tfi2', [r2]),
            C.convertArrays2ZoneNode('tfi3', [r3]),
            C.convertArrays2ZoneNode('tfi4', [r4])]

def TFIMono(a1, a2):
    """Generate a transfinite interpolation mesh from 2 input curves.
    Usage: TFIMono(a)"""
    a1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    a2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    r = Generator.TFIMono(a1, a2)
    return [C.convertArrays2ZoneNode('tfi', [r[0]])]

def TTM(a, niter=100):
    """Smooth a mesh using Thomson-Mastin elliptic generator.
    Usage: TTM(a, niter)"""
    m = C.getFields(Internal.__GridCoordinates__,a)[0]
    m = Generator.TTM(m, niter)
    return C.convertArrays2ZoneNode('ttm', [m])

def hyper2D(t, distrib, type):
    """Generate an hyperbolic mesh. 
    Usage: hyper2D(t, distrib, type)"""
    d = C.getFields(Internal.__GridCoordinates__, distrib)[0]
    return C.TZGC(t, 'nodes', Generator.hyper2D, d, type)

def hyper2D2(t, distrib, type, alpha):
    """Generate an hyperbolic mesh with a constant alpha angle.
    Usage: hyper2D2(array, arrayd, type, alpha)"""
    d = C.getFields(Internal.__GridCoordinates__,distrib)[0]
    return C.TZGC(t, 'nodes', Generator.hyper2D2, d, type, alpha)

def hyper2D3(t, distrib, type, alpha1, alpha2):
    """Generate an hyperbolic mesh with boundary alpha angles.
    Usage: hyper2D3(array, arrayd, type, alpha1, alpha2)"""
    d = C.getFields(Internal.__GridCoordinates__,distrib)[0]
    return C.TZGC(t, 'nodes', Generator.hyper2D3, d, type, alpha1, alpha2)

def hyper2D4(t, distrib, type):
    """Generate an hyperbolic mesh.
    Usage: hyper2D4(array, arrayd, type)"""
    d = C.getFields(Internal.__GridCoordinates__, distrib)[0]
    return C.TZGC(t, 'nodes', Generator.hyper2D4, d, type)

def addNormalLayers(t, distrib, check=0, niter=0):
    """Generate N layers to a surface following normals. Distrib is the 
    height of each layer.
    If niter = 0, the normal are not smoothed; else niter is the number of
    smoothing iterations applied to normals.
    Usage: addNormalLayers(surface, distrib, check, niter)"""
    d = C.getFields(Internal.__GridCoordinates__, distrib)[0]
    tp = Internal.copyRef(t)
    C._deleteZoneBC__(tp)
    C._deleteFlowSolutions__(tp, 'centers')
    C._deleteGridConnectivity__(tp)
    coords = C.getAllFields(tp, 'nodes')
    coords = Generator.addNormalLayers(coords, d, check, niter)
    return C.setFields(coords, tp, 'nodes')

#===============================================================================
# builds an extension starting from a contour c using the normals to surfaces
# and a point distribution dh
#===============================================================================
def buildExtension(c, surfaces, dh, niter=0):
    """Build an extension zone starting from contour c with respect to 
    normals (smoothed niter times) to surfaces.
    Usage: buildExtension(c,surfaces,dh,niter)"""
    cp = Internal.copyRef(c)
    contours = C.getFields(Internal.__GridCoordinates__,cp)[0]
    surfacesA = C.getFields(Internal.__GridCoordinates__,surfaces)
    dhj = C.getFields(Internal.__GridCoordinates__,dh)[0]
    coords = G.buildExtension(contours,surfacesA, dhj, niter)
    return C.setFields([coords], cp, 'nodes')

#==============================================================================
# Walk on surface defined by t starting from a contour 
#==============================================================================
def surfaceWalk(t, contour, distrib, constraints=[], niter=0,
                alphaRef=180., check=0, toldist=1.e-6):
    """Generate a surface mesh by a walk on a list of surfaces, starting from 
    a contour c and following constraints. niter is the number of smoothings.
    if check=1, walk stops before negative cells appear.
    Usage: surfaceWalk(t, contour, distrib, constraints, niter, alphaRef, check)"""
    tp = Internal.copyRef(t)
    surfaces = C.getFields(Internal.__GridCoordinates__, tp)
    c = C.getFields(Internal.__GridCoordinates__, contour)[0]
    d = C.getFields(Internal.__GridCoordinates__, distrib)[0]
    if constraints !=[]: cons = C.getFields(Internal.__GridCoordinates__,constraints)
    else: cons = []
    m = Generator.surfaceWalk(surfaces, c, d, constraints=cons, niter=niter,
                              alphaRef=alphaRef, check=check, toldist=toldist)
    return C.convertArrays2ZoneNode('surface', [m])

#=============================================================================
# Create volume collar grid(s) starting from s1 and s2,
# type: 'union' or 'difference': assembly type between s1 and s2
# return the collar grids(s) and s0 the surface resulting from the
# boolean operation between s1 and s2 as: collar1,collar2,s0
#=============================================================================
def collarMesh(s1, s2, distribj, distribk,
               niterj=100, niterk=100, ext=10, alphaRef=180., type='union',
               contour=[], constraints1=[], constraints2=[], toldist = 1.e-10):
    """Generates a collar mesh starting from s1 and s2 surfaces, distributions
    along the surfaces and along the normal direction, with respect to the assembly type between grids.
    Usage: collarMesh(s1, s2, distribj, distribk, niterj, niterk, ext, alphaRef, type, contour,constraints1,constraints2,toldist)"""
    import Collar
    surf1 = C.getFields(Internal.__GridCoordinates__,s1)
    surf2 = C.getFields(Internal.__GridCoordinates__,s2)
    dj = C.getFields(Internal.__GridCoordinates__,distribj)[0]
    dk = C.getFields(Internal.__GridCoordinates__,distribk)[0]
    if constraints1 != []: constraints11 = C.getFields(Internal.__GridCoordinates__,constraints1)
    else: constraints11 = []
    if constraints2 != []: constraints22 = C.getFields(Internal.__GridCoordinates__,constraints2)
    else: constraints22 = []
    if contour != []: contoura = C.getFields(Internal.__GridCoordinates__,contour)[0]
    else: contoura = contour
    infos = Collar.createCollarMesh__(surf1, surf2, dj, dk, niterj, niterk, 
                                      ext, alphaRef, type,
                                      contoura, constraints11, constraints22, toldist)
    zones = []
    for info in infos:
        ranges = info[1:]
        z = C.convertArrays2ZoneNode('Collar', [info[0]])
        for r in ranges: z = C.addBC2Zone(z,'wall','BCWall',r)
        # ajout des BCOverlaps 
        if type == 'union':
            C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z, 'jmax', trirac=[1,2,3]) 
            C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', trirac=[1,2,3]) 
        elif type =='difference':
            C._addBC2Zone(z, 'match', 'BCMatch', 'imin', z, 'imax', trirac=[1,2,3]) 
            C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', trirac=[1,2,3]) 
        z = C.fillEmptyBCWith(z, 'overlap', 'BCOverlap')               
        zones += [z]
        
    return zones

#=============================================================================
# Generation de grilles cartesiennes coincidentes
#=============================================================================
def gencartmb(t, h, Dfar, nlvl):
    """Generate a Cartesian multiblock set of grids, whose sizes are based
    on the proximity to the body.
    Usage: gencartmb(A, h, Dfar, nlvl)"""
    bodies = C.getFields(Internal.__GridCoordinates__, t)
    cartzones = Generator.gencartmb(bodies, h, Dfar, nlvl)
    c = 1
    zones = []
    for mc in cartzones:
        zone = C.convertArrays2ZoneNode('cart'+str(c), [mc])        
        zones.append(zone)
        c += 1
    return zones

#=============================================================================
# Generation d'un maillage a partir d'une polyline definie par une zone
# Retourne une liste de zones
#=============================================================================
def polyLineMesher(z, h, yplus, density):
    """Generate a multiple mesh for a polyline defined by z.
    Usage: polyLineMesher(z, h, yplus, density)"""
    try:
        import Connector.PyTree as X
        import PolyLine as GP
    except:
        raise ImportError("polyLineMesher: requires Polyline and Connector modules.")
    name = z[0]
    coord = C.getFields(Internal.__GridCoordinates__,z)[0]
    res = GP.polyLineMesher(coord, h, yplus, density)
    allmeshes = res[0]; allwalls = res[1]; hout = res[2]; dout = res[3]
    # creation des maillages
    noz = 1
    zones = []
    for m in allmeshes:
        suffz = '_'+str(noz)
        zone = C.convertArrays2ZoneNode(name+suffz,[m])
        walls = allwalls[noz-1]
        now = 1
        for range in walls :
            bndName = 'wall'+str(noz)+'_'+str(now)
            zone = C.addBC2Zone(zone, bndName, 'BCWall', range)
            now = now + 1
        noz += 1
        zones.append(zone)
    tb = C.newPyTree(['Base']); tb[2][1][2] += zones
    tb = X.connectMatch(tb,dim=2)
    tb = C.fillEmptyBCWith(tb, 'overlap', 'BCOverlap', dim = 2) 
    zones = tb[2][1][2]
    return [zones, hout, dout]

#=============================================================================
# Generation d'un maillage a partir d'une courbe polyC1 definie par une zone
# Retourne une liste de zones
#=============================================================================
def polyC1Mesher(z, h, yplus, density, splitCrit=10., dalpha=5., depth=1):
    """Generate a multiple mesh for a polyC1-curve defined by z.
    Usage: polyC1Mesher(z, h, yplus, density, splitCrit, dalpha, depth)"""
    try: import Connector.PyTree as X; import PolyC1 as GP
    except: raise ImportError("polyC1Mesher: requires PolyC1 and Connector.PyTree modules.")
    name = z[0]
    coord = C.getFields(Internal.__GridCoordinates__,z)[0]
    res = GP.polyC1Mesher(coord, h, yplus, density, splitCrit, dalpha, depth)
    allmeshes = res[0]; allwalls = res[1]; hout = res[2]; dout = res[3]

    # creation des maillages
    noz = 1
    zones = []
    for m in allmeshes:
        suffz = '_'+str(noz)
        zone = C.convertArrays2ZoneNode(name+suffz,[m])
        walls = allwalls[noz-1]
        now = 1
        for range in walls :
            bndName = 'wall'+str(noz)+'_'+str(now)
            zone = C.addBC2Zone(zone, bndName, 'BCWall', range)
            now += 1
        noz += 1
        zones.append(zone)
    tb = C.newPyTree(['Base']); tb[2][1][2] += zones
    tb = X.connectMatch(tb,dim=2)
    tb = C.fillEmptyBCWith(tb, 'overlap', 'BCOverlap', dim=2)
    zones = Internal.getZones(tb)
    return [zones, hout, dout]

#=============================================================================
# Generation d'un maillage a partir d'un polyQuad
# Retourne une liste de zones
#=============================================================================
def polyQuadMesher(z, h, hf, density, next):
    """Generate a multiple mesh for a polyquad defined by z.
     Usage: polyQuadMesher(z, h, hf, density, next)"""
    try:
        import Connector.PyTree as X
        import PolyQuad as GP
    except:
        raise ImportError("polyQuadMesher: requires PolyQuad and Connector modules.")
    name = z[0]
    coord = C.getFields(Internal.__GridCoordinates__,z)[0]
    res = GP.polyQuadMesher(coord, h, hf, density, next)
    allmeshes = res[0]; allwalls = res[1]; hout = res[2]; dout = res[3]
    # creation des maillages
    noz = 1
    zones = []
    for m in allmeshes:
        suffz = '_'+str(noz)
        zone = C.convertArrays2ZoneNode(name+suffz,[m])
        walls = allwalls[noz-1]
        now = 1
        for range in walls :
            bndName = 'wall'+str(noz)+'_'+str(now)
            zone = C.addBC2Zone(zone, bndName, 'BCWall', range)
            now = now + 1
        noz = noz + 1
        zones.append(zone)
    tb = C.newPyTree(['Base']); tb[2][1][2] += zones
    tb = X.connectMatch(tb)
    tb = C.fillEmptyBCWith(tb, 'overlap', 'BCOverlap') 
    zones = tb[2][1][2]
    return [zones, hout, dout]

#=============================================================================
# Generation d'un maillage a partir d'un polyTri
# Retourne une liste de zones
#=============================================================================
def polyTriMesher(z, h, hf, density, next):
    """Generate a multiple mesh for a polytri defined by z.
    Usage: polyTriMesher(z, h, hf, density, next)"""
    try:
        import Connector.PyTree as X
        import PolyTri as GP
    except:
        raise ImportError("polyTriMesher: requires PolyTri and Connector modules.")
    name = z[0]
    coord = C.getFields(Internal.__GridCoordinates__,z)[0]
    res = GP.polyTriMesher(coord, h, hf, density, next)
    allmeshes = res[0]; allwalls = res[1]; hout = res[2]; dout = res[3]
    # creation des maillages
    noz = 1
    zones = []
    for m in allmeshes:
        suffz = '_'+str(noz)
        zone = C.convertArrays2ZoneNode(name+suffz,[m])
        walls = allwalls[noz-1]
        now = 1
        for range in walls :
            bndName = 'wall'+str(noz)+'_'+str(now)
            zone = C.addBC2Zone(zone, bndName, 'BCWall', range)
            now += 1
        noz += 1
        zones.append(zone)
    tb = C.newPyTree(['Base']); tb[2][1][2] += zones
    tb = X.connectMatch(tb)
    tb = C.fillEmptyBCWith(tb, 'overlap', 'BCOverlap', dim=3) 
    zones = tb[2][1][2]
    return [zones, hout, dout]

#------------------------------------------------------------------------------
# Split a i-array and map a distribution on the splitted i-array
#------------------------------------------------------------------------------
def mapSplit(z, d, split_crit, dens_max=1000):
    """Split a i-array defined by a zone and map a distribution d on a curve defined by splitted i-arrays.
    Usage: mapSplit(t, d, split_crit, dens_max)"""
    z = C.deleteFlowSolutions__(z)
    a = C.getFields(Internal.__GridCoordinates__, z)[0]
    name = z[0]
    dist = C.getFields(Internal.__GridCoordinates__, d)[0]
    A = Generator.mapSplit(a, dist, split_crit, dens_max)
    c = 1; zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode(name+str(c),[i])        
        zones.append(zone); c += 1
    return zones

def intersection(surface1, surface2, tol=0.):
    """Compute the intersection between two input closed surfaces.
    Usage: intersection(s1, s2, tol)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    s = Generator.intersection(s1, s2, tol)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the intersection between two input closed entities.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = Generator.booleanIntersection(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the union between two input closed entities.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = Generator.booleanUnion(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('union', [s])

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the difference between the two input closed surfaces.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = Generator.booleanMinus(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('minus', [s])
    
def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Compute the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    sld = C.getFields(Internal.__GridCoordinates__, solid)[0]
    operand = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = Generator.booleanModifiedSolid(operand, sld, tol, preserve_solid, agg_mode)
    return C.convertArrays2ZoneNode('modified_solid', [s])

#------------------------------------------------------------------------------
# Calcul la carte d'orhogonalite d'une grille
# 1D: retourne un tableau d'angles alpha constants egaux a 90 degres
# 2D: retourne un tableau d'angles alpha
# 3D: retourne 3 tableaux d'angles alpha (dans l'ordre alpha_IJ, alpha_IK and alpha_JK pour les grilles structures)
#------------------------------------------------------------------------------
def getOrthogonalityMap(t):
    """Return the orthogonality map in an array.
    Usage: getOrthogonalityMap(t)"""
    return C.TZGC(t, 'centers', Generator.getOrthogonalityMap)

def _getOrthogonalityMap(t):
    return C._TZGC(t, 'centers', Generator.getOrthogonalityMap)

#------------------------------------------------------------------------------
# Calcul de la regularite (ratio entre des mailles adjacentes) d'une grille
# 1D: retourne un champ "reg"
# 2D: retourne deux champs "reg_i", "reg_j"
# 3D: retourne 3 champs (dans l'ordre "reg_i", "reg_j", "reg_k" pour les grilles structures)
#------------------------------------------------------------------------------
def getRegularityMap(t):
    """Return the regularity map in an array.
    Usage: getRegularityMap(t)"""
    return C.TZGC(t, 'centers', Generator.getRegularityMap)

def _getRegularityMap(t):
    return C._TZGC(t, 'centers', Generator.getRegularityMap)
    
#------------------------------------------------------------------------------
# Calcul la qualite pour un maillage TRI (0. triangle degenere, 1. equilateral)
#------------------------------------------------------------------------------
def getTriQualityMap(t):
    """Return the quality map of a TRI array (0. for a degenerated triangle, 1. for an equilateral one).
    Usage: getTriQualityMap(t)"""
    return C.TZGC(t, 'centers', Generator.getTriQualityMap)

def _getTriQualityMap(t):
    return C._TZGC(t, 'centers', Generator.getTriQualityMap)
