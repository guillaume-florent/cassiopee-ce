# Ce module permet de manipuler les SkeletonTree, les PartialTree et les
# BBoxTree
# - Skeleton tree: arbre complet avec les numpy float de taille > 100 non
# charge et valant None
# - Partial tree: arbre partiel ne contenant que les zones chargees
# - BBox tree: arbre identique a un pyTree mais ou les zones sont les
# BBox des zones du pyTree.

import Converter
import Internal
import numpy
import PyTree

#==============================================================================
# Lit un arbre squelette
# Warning: pour l'instant limite a hdf et adf
#==============================================================================
def convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5, 
                             maxDepth=-1):
    """Read a file and return a skeleton tree."""
    return PyTree.convertFile2PyTree(
        fileName, format, skeletonData=[maxFloatSize,maxDepth])

#==============================================================================
# Determine si une zone est une zone squelette
#==============================================================================
def isZoneSkeleton__(z):
    cx = Internal.getNodeFromType1(z, 'GridCoordinates_t')
    if cx is not None:
        t1 = Internal.getNodesFromType1(cx, 'DataArray_t')
        for d in t1:
            if d[1] is None: return True
        return False
    cx = Internal.getNodesFromType1(z, 'FlowSolution_t')
    for x in cx:
        t1 = Internal.getNodesFromType1(x, 'DataArray_t')
        for d in t1:
            if d[1] is None: return True
    return False

#==============================================================================
# Converti un arbre en arbre squelette
# i.e supprime les noeuds DataArray_t et les remplace par None
#==============================================================================
def convert2SkeletonTree(t):
    """Convert a tree to a skeleton tree."""
    tp = Internal.copyRef(t)
    _convert2SkeletonTree(tp)
    return tp

def _convert2SkeletonTree(t):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType(z, 'DataArray_t')
        for n in nodes: n[1] = None
    return None

#==============================================================================
# Converti un arbre squelette charge en arbre partiel
# rank=-1: enleve les zones squelettes
# rank>=0: enleve les zones de proc != rank
#==============================================================================
def convert2PartialTree(t, rank=-1):
    """Convert a tree to a partial tree."""
    tp = Internal.copyRef(t)
    _convert2PartialTree(tp, rank)
    return tp

def _convert2PartialTree(t, rank=-1):
    zones = Internal.getZones(t)
    for z in zones:
        crit = False
        if rank == -1: crit = isZoneSkeleton__(z)
        else:
            proc = Internal.getNodeFromName2(z, 'proc')
            if proc is not None: 
                proc = Internal.getValue(proc) 
                crit = (rank != proc)
        if crit:
            (p, c) = Internal.getParentOfNode(t, z)
            if (Internal.isStdNode(t) == 0 and id(p) == id(t)): del p[c]
            else: del p[2][c]
    return None

#==============================================================================
# Lit des zones et les remplace dans l'arbre
# IN: tree: arbre
# IN: rank: zones a charger
# IN: ou zoneNames: une liste des noms des zones a charger (['Base/Zone'])
# Warning: limite a adf et hdf
# Note: si un noeud Proc existe dans la zone remplacee, il est recopie
# dans la nouvelle zone
#==============================================================================
def readZones(t, fileName, format=None, rank=None, zoneNames=None):
    """Read some zones in a skeleton tree (by rank or name)."""
    tp = Internal.copyRef(t)
    _readZones(tp, fileName, format, rank, zoneNames)
    return tp

def _readZones(t, fileName, format=None, rank=None, zoneNames=None):
  if (zoneNames is None and rank is None): return None
  bases = Internal.getBases(t)
  if rank is not None: # load zones by rank
      # Chemins des zones a remplacer
      paths = []
      for b in bases:
          zones = Internal.getNodesFromType1(b, 'Zone_t')
          for z in zones:
              nproc = Internal.getNodeFromName2(z, 'proc')
              if nproc is not None:
                  nproc = Internal.getValue(nproc)
                  if nproc == rank: paths.append('/'+b[0]+'/'+z[0])
      
  else: # by zone names
      paths = zoneNames[:]
      for c in xrange(len(paths)):
          if paths[c][0] != '/': paths[c] = '/'+paths[c]

  #print 'Reading '+fileName+' '+str(paths)+'...',
  print 'Reading %s [%d zones]...'%(fileName,len(paths)),
  if format is None:
      format = Converter.convertExt2Format__(fileName)
      
  loadedZones = Converter.converter.readPyTreeFromPaths(fileName, paths,
                                                        format)
  
  # Replace/add now loaded zones
  m = 0
  for p in paths:
      z = Internal.getNodeFromPath(t, p)
      if z is not None:
          if rank is not None: # recopie les solveurs data
              nproc = Internal.getNodeFromName2(z, 'proc')
              if nproc is not None:
                  nproc = Internal.getValue(nproc)
              else: nproc = 0
          (p, c) = Internal.getParentOfNode(t, z)
          if (Internal.isStdNode(t) == 0 and id(t) == id(p)):
              p[c] = loadedZones[m]; zr = p[c]
          else: p[2][c] = loadedZones[m]; zr = p[2][c]
          if rank is not None:
              node = Internal.getNodeFromName(zr, '.Solver#Param')
              if node is not None: param = node
              else:
                  param = ['.Solver#Param', None, [], 'UserDefinedData_t']
                  zr[2].append(param)
              v = numpy.zeros((1,1), numpy.int32); v[0,0] = nproc
              node = Internal.getNodeFromName(param, 'proc')
              if node is not None:
                  node[1] = v
              else:
                  a = ['proc', v, [], 'DataArray_t']
                  param[2].append(a)
          m += 1
  print 'done.'
  return None

#==============================================================================
# Ecrit des zones dans un fichier deja cree
# Warning: limite a adf et hdf
#==============================================================================
def writeZones(t, fileName, format=None, proc=None, zoneNames=None):
    """Write some zones in an existing file (adf or hdf)."""
    if (zoneNames is None and proc is None): return None
    tp, ntype = Internal.node2PyTree(t)
    bases = Internal.getBases(tp)
    if proc is not None: # write zones by proc
        # Chemins des noeuds a remplacer (Zone et IntegralData)
        paths = []; nodes = []
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t') + Internal.getNodesFromType1(b, 'IntegralData_t')
            for z in zones:
                nproc = Internal.getNodeFromName2(z, 'proc')
                if nproc is not None:
                    nproc = Internal.getValue(nproc)
                    if nproc == proc:
                        paths.append('/%s/%s'%(b[0],z[0]))
                        nodes.append(z)
                else: # write nevertheless
                    paths.append('/%s/%s'%(b[0],z[0]))
                    nodes.append(z)
    else: # by zone names
        paths = zoneNames[:]
        nodes = []
        for p in paths:
            n = Internal.getNodeFromPath(tp, p)
            nodes.append(n)
        for c in xrange(len(paths)):
            if paths[c][0] != '/': paths[c] = '/'+paths[c]
            
    #print 'Writing '+fileName+' '+str(paths)+'...',
    print 'Writing %s [%d zones]...'%(fileName,len(paths)),
    if format is None:
        format = Converter.convertExt2Format__(fileName)
    Converter.converter.writePyTreePaths(fileName, nodes, paths, format)
    print 'done.'
    return None

#==============================================================================
# zones est une liste de zones, mais triees par procs
# Ex: [ [zonesDeProc0], [zonesDeProc1], ...]
# IN: liste des zones a remplacer dans t
# IN: size: nbre de processeurs
# Met les zones dans l'arbre par identification des noms
#==============================================================================
def setZonesInTree(t, zones):
    tp = Internal.copyRef(t)
    _setZonesInTree(tp, zones)
    return tp

def _setZonesInTree(t, zones):
    size = len(zones)
    for i in xrange(size):
        for j in xrange(len(zones[i])):
            zone = zones[i][j]
            zoneName = zone[0]
            z = Internal.getNodeFromName2(t, zoneName)
            if z is not None:
                (p, nb) = Internal.getParentOfNode(t, z)
                p[2][nb] = zone
            else:
                # append to given base
                bases = Internal.getBases(t)
                bases[0][2] += [zone]
    return None

#==============================================================================
# update le graph si zoneName sur proc est a envoyer a popp
#==============================================================================
def updateGraph__(graph, proc, popp, zoneName):
    if popp != proc:
        if proc not in graph: graph[proc] = {popp:[zoneName]}
        else:
            if popp not in graph[proc]: graph[proc][popp] = [zoneName]
            else:
                if zoneName not in graph[proc][popp]:
                    graph[proc][popp].append(zoneName)
    return graph

#==============================================================================
# Retourne le proc de z si z est local
#==============================================================================
def getProcLocal__(z, procDict):
    if procDict is not None: return procDict[z[0]]
    else:
        proc = Internal.getNodeFromName2(z, 'proc')
        proc = Internal.getValue(proc)
        return proc
        
#==============================================================================
# Retourne le proc de z si z est global
#==============================================================================
def getProcGlobal__(zoneName, t, procDict):
    if procDict is not None: return procDict[zoneName]
    else:
        z = Internal.getNodeFromName2(t, zoneName)
        proc = Internal.getNodeFromName2(z, 'proc')
        proc = Internal.getValue(proc)
        return proc

#==============================================================================
# Calcule le graph
# IN: t: arbre contenant des noeuds 'proc' entierement charge
# IN: type: type de graph voulu
# type='bbox' si intersection de bbox des zones (full/bbox)
# type='bbox2' si intersection de bbox et pas sur la meme base (full/bbox)
# type='bbox3' si intersection de bbox avec t2 (t:full/bbox et t2:full/bbox)
# type='match' si match entre zones (full/skel/load skel/partial[+procDict])
# type='ID' si donnees d'interpolation entre zones (full/skel/load skel/partial+procDict)
# type='IBCD' si donnees IBC entres zones (full/skel/load skel/partial+procDict)
# type='ALLD' si toutes donnees (Interp+IBC) (full/skel/load skel/partial+procDict)
# type='proc' si la zone a un noeud proc different du proc courant (full/skel/load skel/partial)
# intersectionsDict: dictionnaire Python contenant la liste de zones intersec-
# tantes, comme produit par "X.getIntersectingDomains". Attention, si type='bbox3',
# l'utilisateur doit fournir l'arbre t2 et intersectionDict doit decrire les 
# intersections entre t et t2, comme produit par "X.getIntersectingDomains(t,t2)".
# OUT: graph: dictionnaire contenant des informations d'envoie
# des zones entre processeurs
# graph est construit de telle sorte que:
# graph[P0][P1] renvoie la liste des zones du processeur P0
# "connectees" avec au moins une zone du processeur P1
# Ex: addXZones: envoie les zones graph[rank][opp] au proc opp,
# attend ensuite les zones graph[opp][rank] pour tout opp.
#==============================================================================
def computeGraph(t, type='bbox', t2=None, procDict=None, rank=0,
                 intersectionsDict=None):
    """Return the communication graph for different block relation types."""
    zones = Internal.getZones(t)
    graph = {}
    if type == 'bbox': # zone de t de P0 intersectent une zone de t de P1
        if not intersectionsDict:
            try: import Connector.PyTree as X
            except: raise ImportError("computeGraph: requires Connector module.")
            intersectionsDict = X.getIntersectingDomains(t)        
        for z in zones:
            proc = getProcLocal__(z, procDict)
            for z2 in zones:
                if z[0] in intersectionsDict[z2[0]]:
                    popp = getProcGlobal__(z2[0], t, procDict) 
                    updateGraph__(graph, proc, popp, z[0])         
        #import Connector.PyTree as X 
        #for z in zones: 
        #    proc = getProcLocal__(z, procDict) 
        #    doms = X.getBBIntersectingDomains(t, z, tol=1.e-12) 
        #    for d in doms: 
        #        popp = getProcGlobal__(d[0], t, procDict)  
        #        updateGraph__(graph, proc, popp, z[0]) 

    elif type == 'bbox2': # zone de t de P0 qui intersecte une zone de t de P1 mais qui n'est pas sur la meme base
        if not intersectionsDict:
            try: import Connector.PyTree as X
            except: raise ImportError("computeGraph: requires Connector module.")
            intersectionsDict = X.getIntersectingDomains(t)
        for z in zones:
            (p, c) = Internal.getParentOfNode(t, z)
            base = p[0]
            proc = getProcLocal__(z, procDict)
            for z2 in zones:
                if z[0] in intersectionsDict[z2[0]]:
                    (p, c) = Internal.getParentOfNode(t, z2)
                    based = p[0]
                    popp = getProcGlobal__(z2[0], t, procDict) 
                    if (popp != proc and base != based):
                        if proc not in graph: graph[proc] = {popp:[z[0]]}
                        else:
                            if popp not in graph[proc]: graph[proc][popp] = [z[0]]
                            else:
                                if z[0] not in graph[proc][popp]:
                                    graph[proc][popp].append(z[0])
        # import Connector.PyTree as X
        # for z in zones:
        #     (p, c) = Internal.getParentOfNode(t, z)
        #     base = p[0]
        #     proc = getProcLocal__(z, procDict)
        #     doms = X.getBBIntersectingDomains(t, z, tol=1.e-12)
        #     for d in doms:
        #         (p, c) = Internal.getParentOfNode(t, d)
        #         based = p[0]
        #         popp = getProcGlobal__(d[0], t, procDict) 
        #         if (popp != proc and base != based):
        #             if proc not in graph: graph[proc] = {popp:[z[0]]}
        #             else:
        #                 if popp not in graph[proc]: graph[proc][popp] = [z[0]]
        #                 else:
        #                     if z[0] not in graph[proc][popp]:
        #                         graph[proc][popp].append(z[0])

    elif type == 'bbox3': # zone de t sur P0 qui intersecte une zone de t2 sur P1
        #if not t2: raise ValueError("computeGraph: type bbox3 requires a t2.")
        zones2 = Internal.getZones(t2)
        if not intersectionsDict:
            try: import Connector.PyTree as X
            except: raise ImportError("computeGraph: requires Connector module.")
            intersectionsDict = X.getIntersectingDomains(t, t2)
        for z in zones:
            proc = getProcLocal__(z, procDict)
            for z2 in zones2:
                if z2[0] in intersectionsDict[z[0]]:
                    popp = getProcGlobal__(z2[0], t2, None)
                    updateGraph__(graph, proc, popp, z[0])
        #import Connector.PyTree as X
        #for z in zones:
        #    proc = getProcLocal__(z, procDict)
        #    doms = X.getBBIntersectingDomains(t2, z, tol=1.e-12)
        #    for d in doms:
        #        popp = getProcGlobal__(d[0], t2, None)
        #        updateGraph__(graph, proc, popp, z[0])

    elif type == 'ID': # base sur les interpolations data
        for z in zones:
            proc = getProcLocal__(z, procDict)
            subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
            subRegions = []
            for s in subRegions2:
                sname = s[0]
                if sname.split('_')[0] == 'ID': subRegions.append(s)
            for s in subRegions:
                donor = Internal.getValue(s)
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    popp = getProcGlobal__(donor, t, procDict)
                    updateGraph__(graph, proc, popp, z[0])

    elif type == 'IBCD': # base sur les IBC data
        for z in zones:
            proc = getProcLocal__(z, procDict)
            subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
            subRegions = []
            for s in subRegions2:
                sname = s[0]
                if sname.split('_')[0] == 'IBCD': subRegions.append(s)
            for s in subRegions:
                donor = Internal.getValue(s)
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des IBC
                    popp = getProcGlobal__(donor, t, procDict)
                    updateGraph__(graph, proc, popp, z[0])

    elif type == 'ALLD': # base sur les Interpolations+IBC data
        for z in zones:
            proc = getProcLocal__(z, procDict)
            subRegions = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
            for s in subRegions:
                donor = Internal.getValue(s)
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations/IBC
                    popp = getProcGlobal__(donor, t, procDict)
                    updateGraph__(graph, proc, popp, z[0])

    elif type == 'match': # base sur les raccords matchs
        for z in zones:
            proc = getProcLocal__(z, procDict)
            GC = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for c in GC:
                donor = Internal.getValue(c)
                popp = getProcGlobal__(donor, t, procDict)
                updateGraph__(graph, proc, popp, z[0])

    elif type == 'proc':
        for z in zones:
            if not isZoneSkeleton__(z):
                popp = getProcLocal__(z, procDict)
                proc = rank
                updateGraph__(graph, proc, popp, z[0])
    return graph

#==============================================================================
# Retourne le dictionnaire proc['blocName']
# a partir d'un arbre distribue contenant des noeuds proc
#==============================================================================
def getProcDict(t):
    """Return the dictionary proc['zoneName']."""
    proc = {}
    zones = Internal.getZones(t)
    for z in zones:
        nproc = getProc(z)
        proc[z[0]] = nproc
    return proc

#==============================================================================
# getProc in zone (if exists), otherwise return -1
# IN: a: zone node
#==============================================================================
def getProc(a):
    """Return the proc of zone."""
    nproc = Internal.getNodeFromName2(a, 'proc')
    if nproc is not None: nproc = Internal.getValue(nproc)
    else: nproc = -1
    return nproc

#==============================================================================
# setProc in zone
# IN: z: zone node
# IN: nproc: proc to set (int)
#==============================================================================
def setProc(t, nproc):
  """Set the proc to a zone or a set of zones."""
  tp = Internal.copyRef(t)
  _setProc(tp, nproc)
  return tp

def _setProc(t, nproc):
    zones = Internal.getZones(t)
    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#Param')
        if node is not None: param = node
        else:
            param = ['.Solver#Param', None, [], 'UserDefinedData_t']
            z[2].append(param)
        v = numpy.zeros((1,1), numpy.int32); v[0,0] = nproc
        node = Internal.getNodeFromName(param, 'proc')
        if node is not None: node[1] = v
        else:
            a = ['proc', v, [], 'DataArray_t']
            param[2].append(a)
    return None
 
