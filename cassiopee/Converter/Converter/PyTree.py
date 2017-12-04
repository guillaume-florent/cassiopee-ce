import Converter
import Internal
import numpy, fnmatch
__version__ = Converter.__version__

# Variables globales
__ZoneNameServer__ = {}
__BCNameServer__ = {}
__BaseNameServer__ = {}

#==============================================================================
# -- Gestion du nommage unique --
#==============================================================================

# -- getZoneName
# Retourne un nom de zone (unique par rapport au __ZoneNameServer__)
# IN: proposedName: nom propose
# OUT: nouveau nom unique
def getZoneName(proposedName):
  global __ZoneNameServer__
  (name, __ZoneNameServer__) = getUniqueName(proposedName, __ZoneNameServer__)
  return name

# Enregistre les noms de zones de t dans le __ZoneNameServer__
def registerZoneNames(t):
  global __ZoneNameServer__
  nodes = Internal.getZones(t)
  for i in nodes:
    __ZoneNameServer__[i[0]] = 0

# -- getBCName
# Retourne un nom de BC unique, stocke les noms de BC_t, GridConnectivity_t
# IN: proposedName: nom propose
def getBCName(proposedName):
  global __BCNameServer__
  (name, __BCNameServer__) = getUniqueName(proposedName, __BCNameServer__)
  return name

# Enregistre les noms des BCs de t dans le __BCNameServer__
def registerBCNames(t):
  global __BCNameServer__
  nodes = Internal.getNodesFromType(t, 'BC_t')
  for i in nodes: __BCNameServer__[i[0]] = 0
  nodes = Internal.getNodesFromType(t, 'GridConnectivity1to1_t')
  for i in nodes: __BCNameServer__[i[0]] = 0
  nodes = Internal.getNodesFromType(t, 'GridConnectivity_t')
  for i in nodes: __BCNameServer__[i[0]] = 0

# -- getBaseName
# Retourne un nom de base (unique)
# IN: proposedName: nom propose
# OUT: nouveau nom unique
def getBaseName(proposedName):
  global __BaseNameServer__
  (name, __BaseNameServer__) = getUniqueName(proposedName, __BaseNameServer__)
  return name

# Enregistre les noms de base de t dans le __BaseNameServer__
def registerBaseNames(t):
  global __BaseNameServer__
  nodes = Internal.getBases(t)
  for i in nodes:
    __BaseNameServer__[i[0]] = 0

# Enregistre les Zone names, les Base names, les BC names
def registerAllNames(t):
  registerZoneNames(t); registerBaseNames(t); registerBCNames(t)

# Clear all names
def clearAllNames():
  global __ZoneNameServer__
  __ZoneNameServer__ = {}
  global __BCNameServer__
  __BCNameServer__ = {}
  global __BaseNameServer__
  __BaseNameServer__ = {}
  return None

# Retourne proposedName ou proposedName.count
def getUniqueName(proposedName, server):
  namespl = proposedName.rsplit('.', 1)
  if len(namespl) == 2:
    try: c = int(namespl[1]); name = namespl[0]
    except: name = proposedName
  else: name = proposedName
  if not server.has_key(name):
    server[name] = 0
    return (name, server)
  else:
    c = server[name]; ret = 1
    while ret == 1:
      name2 = '%s.%d'%(name,c)
      if not server.has_key(name2): ret = 0
      else: ret = 1
      c += 1
    server[name2] = 0
    server[name] = c
    return (name2, server)

#==============================================================================
# -- pyTree informations --
#==============================================================================

# -- printTree
# print a tree to file or screen (obsolete, use Internal)
def printTree(t, file=None, stdOut=None, editor=False):
  Internal.printTree(t, file, stdOut, editor)

# -- getNPts
# Retourne le nombre de pts dans t
def getNPts(t):
  """Return the number of points in t.
  Usage: getNpts(t)"""
  zones = Internal.getZones(t)
  npts = 0
  for z in zones:
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Structured': npts += dim[1]*dim[2]*dim[3]
    elif dim[0] == 'Unstructured': npts += dim[1]
  return npts

# -- getNCells
# Retourne le nombre de cellules dans t
def getNCells(t):
  """Return the number of cells in t.
  Usage: getNCells(t)"""
  zones = Internal.getZones(t)
  ncells = 0
  for z in zones:
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Structured':
      ni1 = max(dim[1]-1,1); nj1 = max(dim[2]-1,1); nk1 = max(dim[3]-1,1)
      ncells += ni1*nj1*nk1
    elif dim[0] == 'Unstructured': ncells += dim[2]
  return ncells

# -- convertPyTree2ZoneNames
# Return the list of path of zones of a python tree as a list of strings
def convertPyTree2ZoneNames(t):
  """Return the list of zone names of a py tree.
  Usage: convertPyTree2ZoneNames(t)"""
  l = []
  toptree = Internal.isTopTree(t)
  if toptree:
    bases = Internal.getBases(t)
    for b in bases:
      baseName = b[0]
      zones = Internal.getNodesFromType1(b, 'Zone_t')
      for z in zones:
        zoneName = z[0]
        l.append('%s/%s'%(baseName,zoneName))
  return l

# -- getStdNodesFromName. Applique uniquement sur une zone.
# Retourne les noeuds des champs de nom 'name' dans les conteneurs standards
# Accepte les variables 'centers:var', les noms de containers (GridCoordinates)
# Si pas trouve, retourne []
def getStdNodesFromName(z, name):
  loc = '*'
  v = name.split(':')
  if len(v) > 1: var = v[1]; loc = v[0]
  else: var = v[0]
  result = []
  if loc == 'nodes' or loc == '*':
    node = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if node is not None:
      if name == Internal.__GridCoordinates__: result.append(node)
      for j in node[2]:
        if j[0] == var: result.append(j); break
    node = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if node is not None:
      if name == Internal.__FlowSolutionNodes__: result.append(node)
      for j in node[2]:
        if (j[0] == var): result.append(j); break
  if loc == 'centers' or loc == '*':
    node = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if node is not None:
      if name == Internal.__FlowSolutionCenters__: result.append(node)
      for j in node[2]:
        if j[0] == var: result.append(j); break
  return result

# -- getVarNames
# Retourne une liste des noms de variables presents pour chaque zone
# avec leur localisation
# localisation=centers: ou nom du container:
# Ex: pour t contenant 2 zones:
# [ ['CoordinateX','CoordinateY'], ['centers:F'] ]
# Si excludeXYZ=True, les coordonnees ne font pas partie des champs retournes
# Si loc='both', retourne les variables en centres et en noeuds
# Sinon loc='centers' ou loc='nodes'
def getVarNames(t, excludeXYZ=False, loc='both'):
  allvars = []
  nodes = Internal.getZones(t)
  for z in nodes:
    varlist = []
    if not excludeXYZ:
      nodesGC = Internal.getNodesFromType1(z, 'GridCoordinates_t')
      for i in nodesGC:
        nodesVar = Internal.getNodesFromType1(i, 'DataArray_t')
        for j in nodesVar: varlist.append(j[0])

    if loc == 'nodes' or loc == 'both':
      nodesSol = Internal.getNodesFromName1(z, Internal.__FlowSolutionNodes__)
      location = ''
      for i in nodesSol:
         nodesVar = Internal.getNodesFromType1(i, 'DataArray_t')
         for j in nodesVar: varlist.append(location+j[0])

    if loc == 'centers' or loc == 'both':
      nodesSol = Internal.getNodesFromName1(z, Internal.__FlowSolutionCenters__)
      location = 'centers:'
      for i in nodesSol:
        nodesVar = Internal.getNodesFromType1(i, 'DataArray_t')
        for j in nodesVar: varlist.append(location+j[0])
    allvars.append(varlist)
  return allvars

# -- isNamePresent
# isNamePresent dans un arbre t
# Retourne -1: la variable n'est presente dans aucune zone de t
# Retourne 0: la variable est presente dans au moins une zone, mais pas toutes.
# Retourne 1: la variables est presente dans toutes les zones
def isNamePresent(t, varname):
  vars = getVarNames(t)
  if len(vars) == 0: return -1
  one = 0
  for z in vars:
    found = 0
    for v in z:
      if v == varname: found = 1; one = 1; break
    if found == 0:
      if one == 1: return 0
  if one == 0: return -1
  else: return 1

# -- getNobNozOfZone
# IN: a: zone
# IN: t: le top tree de a
# OUT: (nob, noz): no de la base et no de la zone a dans toptree
# If zone can not be found, return (-1, -1)
# Must be faster than getParentOfNode
def getNobNozOfZone(a, t):
   """Return the (nob, noz) of a in t.
    Usage: getNobNozOfZone(a, t)"""
   nob = 0; noz = 0
   for b in t[2]:
     noz = 0
     for z in b[2]:
       if id(z) == id(a): return (nob, noz)
       noz += 1
     nob += 1
   return (-1, -1)

# -- getNobOfBase
# IN: b: base
# IN: t: le top tree de a
# OUT: nob: no de la base dans le top tree
# If base can not be found, return -1
# Must be faster than getParentOfNode
def getNobOfBase(base, t):
  """Return the nob of a base in t.
  Usage: getNobOfBase(base, t)"""
  nob = 0
  for b in t[2]:
    if id(b) == id(base): return nob
    nob += 1
  return -1

# -- GetZoneNames
# Retourne une liste contenant le nom des zones de l'arbre, trie ainsi:
# En premier, les zones structurees
# En second, les zones non structurees
def getZoneNames(t, prefixByBase=True):
    bases = Internal.getBases(t)
    names = []
    if (bases == [] or prefixByBase == False): # Juste zone name
        zones = Internal.getZones(t)
        for z in zones:
            type = Internal.getZoneType(z)
            if type == 1: names.append(z[0])
        for z in zones:
            type = Internal.getZoneType(z)
            if type == 2: names.append(z[0])
    else: # Base/Zone name
        for b in bases:
            baseName = b[0]
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                type = Internal.getZoneType(z)
                if type == 1: names.append(baseName+Internal.SEP1+z[0])
        for b in bases:
            baseName = b[0]
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                type = Internal.getZoneType(z)
                if type == 2: names.append(baseName+Internal.SEP1+z[0])
    return names

# -- getConnectedZones
# IN: a: zone, tree, base, liste de zones
# IN: toptree: le toptree de a
# OUT: la liste des zones connectees a a par:
# - des matchs
# - des near matchs
# Attention: le nom des zones doit etre unique dans l'arbre!
def getConnectedZones(a, topTree, type='all'):
  """Return the list of zones connected to a through match or nearMatch.
  Usage: getConnectedZones(a, topTree, type)"""
  zones = Internal.getZones(a)
  out = set()
  match = 0
  if type == 'all' or type == 'BCMatch': match = 1
  nearMatch = 0
  if type == 'all' or type == 'BCMatch': nearMatch = 1

  if match == 1:
    for z in zones:
      nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
      for n in nodes:
        donor = Internal.getValue(n)
        out.add(donor)

  if nearMatch == 1:
    for z in zones:
      nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
      for n in nodes:
        donor = Internal.getValue(n)
        out.add(donor)

  # Retrouve les zones
  outz = []
  zones = Internal.getZones(topTree)
  for i in out:
    for z in zones:
      if z[0] == i: outz.append(z); break
  return outz

#==============================================================================
# -- Get / Set field values --
#==============================================================================

# -- getValue
# Retourne la valeur d'un champ standard pour un indice donne
def getValue(t, var, ind):
  """Return the values for a point of index ind or (i,j,k)...
  Usage: getValue(t, var, ind)"""
  if isinstance(var, list):
    out = []
    for i in var:
      ret = getValue__(t, i, ind)
      if isinstance(ret, list): out += ret
      else: out.append(ret)
    return out
  else: return getValue__(t, var, ind)

def getValue__(t, var, ind):
  result = []
  zones = Internal.getZones(t)
  if zones == []: raise ValueError("getValue: not a zone node.")
  z = zones[0]

  # localisation
  u = var.split(':')
  if len(u) >= 2: loc = u[0]; var = u[1]
  else: loc = '*'; var = u[0]

  # dim
  dim = Internal.getZoneDim(z); cellDim = dim[4]
  if dim[0] == 'Structured': # output im,jm,km
    ni = dim[1]; nj = dim[2]; ni1 = max(ni-1,1); nj1 = max(nj-1,1)
    ninj = ni*nj; ni1nj1 = ni1*nj1
    if isinstance(ind, tuple):
      if len(ind) == 3:
        im = ind[0]-1; jm = ind[1]-1; km = ind[2]-1
      elif len(ind) == 2:
        im = ind[0]-1; jm = ind[1]-1; km = 0
      elif len(ind) == 1:
        im = ind[0]-1; jm = 0; km = 0
      else:
        raise ValueError("getValue: too much values in tuple.")
    else:
      if (loc == 'nodes' or loc == '*'):
        km = ind / (ninj)
        jm = (ind - km*ninj) / ni
        im = ind - jm*ni - km*ninj
      else:
        km = ind / (ni1nj1)
        jm = (ind - km*ni1nj1) / ni1
        im = ind - jm*ni1 - km*ni1nj1

    # GridCoordinates
    v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__GridCoordinates__)
            and (loc == 'nodes' or loc == '*')
            and i[3] == 'DataArray_t'):
          if cellDim == 3:
            result.append(i[1][im,jm,km])
          elif cellDim == 2: result.append(i[1][im,jm])
          else: result.append(i[1][im])

    # FlowSolutionNodes
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
            and (loc == 'nodes' or loc == '*')
            and i[3] == 'DataArray_t'):
          if cellDim == 3: result.append(i[1][im,jm,km])
          elif cellDim == 2: result.append(i[1][im,jm])
          else: result.append(i[1][im])

    # FlowSolutionCenters
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
            and (loc == 'centers' or loc == '*')
            and i[3] == 'DataArray_t'):
          if cellDim == 3: result.append(i[1][im,jm,km])
          elif cellDim == 2: result.append(i[1][im,jm])
          else: result.append(i[1][im])

  else: # output index
    # GridCoordinates
    v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__GridCoordinates__)
            and (loc == 'nodes' or loc == '*')
            and i[3] == 'DataArray_t'):
          result.append(i[1][ind])

    # FlowSolutionNodes
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
            and (loc == 'nodes' or loc == '*')
            and i[3] == 'DataArray_t'):
          result.append(i[1][ind])

    # FlowSolutionCenters
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
            and (loc == 'centers' or loc == '*')
            and i[3] == 'DataArray_t'):
          result.append(i[1][ind])

  if len(result) == 1: return result[0]
  else: return result

# -- setValue
# Set a value in a zone field pour l'indice donne
# Attention: pas de copie ici, c'est directement le champ qui est modifie
def setValue(t, var, ind, val):
  """Set the values in an array for a point of index ind or (i,j,k)...
  Usage: setValue(t, var, ind, value)"""
  nodes = Internal.getZones(t)
  if (nodes == []): raise ValueError("setValue: not a zone node.")
  z = nodes[0]

  # localisation
  u = var.split(':')
  if (len(u) >= 2): loc = u[0]; var = u[1]
  else: loc = '*'; var = u[0]

  c = 0
  if not isinstance(val, list): val = [val]

  # dim
  dim = Internal.getZoneDim(z); cellDim = dim[4]
  if (dim[0] == 'Structured'): # output im,jm,km
    ni = dim[1]; nj = dim[2]; ni1 = max(ni-1,1); nj1 = max(nj-1,1)
    ninj = ni*nj; ni1nj1 = ni1*nj1
    if isinstance(ind, tuple):
      ll = len(ind)
      if (ll == 3):
        im = ind[0]-1; jm = ind[1]-1; km = ind[2]-1
      elif (ll == 2):
        im = ind[0]-1; jm = ind[1]-1; km = 0
      elif (ll == 1):
        im = ind[0]-1; jm = 0; km = 0
      else:
        raise ValueError("setValue: too much values in tuple.")
    else:
      if (loc == 'nodes' or loc == '*'):
        km = ind / (ninj)
        jm = (ind - km*ninj) / ni
        im = ind - jm*ni - km*ninj
      else:
        km = ind / (ni1nj1)
        jm = (ind - km*ni1nj1) / ni1
        im = ind - jm*ni1 - km*ni1nj1

    # GridCoordinates
    v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__GridCoordinates__)
        and (loc == 'nodes' or loc == '*')):
          if (cellDim == 3):
            i[1][im,jm,km] = val[c]; c = c+1
          elif (cellDim == 2):
            i[1][im,jm] = val[c]; c = c+1
          else: i[1][im] = val[c]; c = c+1

    # FlowSolutionNodes
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
        and (loc == 'nodes' or loc == '*')):
          if (cellDim == 3):
            i[1][im,jm,km] = val[c]; c += 1
          elif (cellDim == 2):
            i[1][im,jm] = val[c]; c += 1
          else: i[1][im] = val[c]; c += 1

    # FlowSolutionCenters
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
        and (loc == 'centers' or loc == '*')):
          if (cellDim == 3):
            i[1][im,jm,km] = val[c]; c += 1
          elif (cellDim == 2):
            i[1][im,jm] = val[c]; c += 1
          else: i[1][im] = val[c]; c += 1

  else:
    # GridCoordinates
    v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__GridCoordinates__)
        and (loc == 'nodes' or loc == '*')):
          i[1][ind] = val[c]; c += 1

    # FlowSolutionNodes
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
        and (loc == 'nodes' or loc == '*')):
          i[1][ind] = val[c]; c += 1

    # FlowSolutionCenters
    v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if v is not None:
      for i in v[2]:
        if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
        and (loc == 'centers' or loc == '*')):
          i[1][ind] = val[c]; c += 1
  return None

#==============================================================================
# -- Create PyTree --
#==============================================================================

# -- newPyTree
def newPyTree(args=[]):
  """Create a new PyTree.
  Usage: newPyTree(['Base',3, z1, ...])"""
  t = Internal.createRootNode()
  t[2].append(Internal.createCGNSVersionNode())
  l = len(args)
  base = None; cellDim = 3
  for i in args:
    if isinstance(i, str): # a baseName
      base = Internal.createBaseNode(i, cellDim); t[2].append(base)
    elif isinstance(i, int): # base dim
      base[1][0] = i
    else:
      if len(i) == 4: # maybe a standard node
        if i[3] == 'Zone_t':
          if base is None:
            base = Internal.createBaseNode(name, cellDim); t[2].append(base)
          base[2].append(i)
        elif i[3] == 'CGNSBase_t':
          base = i; t[2].append(base)
      else: # a list of zones?
        for z in i:
          if isinstance(z, list):
            if len(z) == 4 and z[3] == 'Zone_t':
              if base is None:
                base = Internal.createBaseNode(name, cellDim); t[2].append(base)
              base[2].append(z)
  return t

# -- addBase2PyTree
# IN: cellDim=1,2,3 (dimension des cellules des zones de cette base)
def addBase2PyTree(t, baseName, cellDim=3):
  """Add a base name to a pyTree.
  Usage: addBase2PyTree(t, baseName, cellDim)"""
  if not isinstance(baseName, str):
    raise TypeError("addBase2PyTree: baseName must be a string.")
  a = Internal.copyRef(t)
  if a == []:
    a = Internal.createRootNode()
    a[2].append(Internal.createCGNSVersionNode())
    base = Internal.createBaseNode(baseName, cellDim)
    a[2].append(base)
  else: _addBase2PyTree(a, baseName, cellDim)
  return a

def _addBase2PyTree(a, baseName, cellDim=3):
  bases = Internal.getBases(a)
  found = False
  for i in bases:
    if i[0] == baseName: found = True
  if not found:
    base = Internal.createBaseNode(baseName, cellDim)
    a[2].append(base)
  return None

#==============================================================================
# -- delete certain nodes --
#==============================================================================

# -- deleteFlowSolutions__
# Enleve les noeuds FlowSolutionNodes, ou FlowSolutionCenters ou les 2
# loc='nodes', 'centers', 'both' respectivement
def deleteFlowSolutions__(a, loc='both'):
  b = Internal.copyRef(a)
  _deleteFlowSolutions__(b, loc)
  return b

def _deleteFlowSolutions__(a, loc='both'):
  if loc == 'centers' or loc == 'both':
    centers = Internal.getNodesFromName3(a, Internal.__FlowSolutionCenters__)
    if centers != []: _rmNodes(a, Internal.__FlowSolutionCenters__)
  if loc == 'nodes' or loc == 'both':
    nodes = Internal.getNodesFromName3(a, Internal.__FlowSolutionNodes__)
    if nodes != []: _rmNodes(a, Internal.__FlowSolutionNodes__)
  return None

# -- deleteGridConnectivity__
# Enleve les noeuds GridConnectivity
# type='None' par defaut: enleve le noeud complet
# type='BCMatch': enleve les BCMatch.
#   si kind='self' detruit uniquement les match sur la meme zone
#   si kind='other', n'enleve pas les raccords coincidents sur la meme zone
# type='BCOverlap' : enleve uniquement les BCOverlap
#   si kind='self' detruit les BCOverlap en autoattach (donorZoneName=zoneName)
#   si kind='other' detruit les BCOverlap pas autoattach
# Par defaut, kind='all', detruit alors les 2 (self + other)
def deleteGridConnectivity__(a, type='None', kind='all'):
  if type == 'None':
    cn = Internal.getNodesFromName3(a, 'ZoneGridConnectivity')
    if cn != []: a = rmNodes(a, 'ZoneGridConnectivity')
  elif type == 'BCMatch':
    if kind == 'self': _deleteSelfBCMatch__(a)
    elif kind == 'other': _deleteOtherBCMatch__(a)
    else: a = rmBCOfType(a, 'BCMatch')
  elif type == 'BCOverlap':
    if kind == 'other': _deleteBCOverlapWithDonorZone__(a)
    elif kind == 'self': _deleteBCOverlapWithoutDonorZone__(a)
    else: a = rmBCOfType(a, 'BCOverlap')
  return a

def _deleteGridConnectivity__(a, type='None', kind='all'):
  if type == 'None':
    cn = Internal.getNodesFromName3(a, 'ZoneGridConnectivity')
    if cn != []: _rmNodes(a, 'ZoneGridConnectivity')
  elif type == 'BCMatch':
    if kind == 'self': _deleteSelfBCMatch__(a)
    elif kind == 'other': _deleteOtherBCMatch__(a)
    else: _rmBCOfType(a, 'BCMatch')
  elif type == 'BCOverlap':
    if kind == 'other': _deleteBCOverlapWithDonorZone__(a)
    elif kind == 'self': _deleteBCOverlapWithoutDonorZone__(a)
    else: _rmBCOfType(a, 'BCOverlap')
  return None

# enleve les BCMatch sur une meme zone
def _deleteSelfBCMatch__(a):
  zones = Internal.getZones(a)
  for z in zones:
    zoneName = z[0]
    nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for i in nodes:
      donorName = i[1]
      if donorName == zoneName:
        (parent, d) = Internal.getParentOfNode(z, i)
        del parent[2][d]
  return None

# enleve les BCMatch entre deux zones differentes
def _deleteOtherBCMatch__(a):
  zones = Internal.getZones(a)
  for z in zones:
    zoneName = z[0]
    nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for i in nodes:
      donorName = i[1]
      if (donorName != zoneName):
        (parent, d) = Internal.getParentOfNode(z, i)
        del parent[2][d]
  return None

# enleve les BCOverlap de type autoattach
def _deleteBCOverlapWithoutDonorZone__(a):
  zones = Internal.getZones(a)
  for z in zones:
    zoneName = z[0]
    nodes = Internal.getNodesFromType1(z, 'GridConnectivity_t')
    for i in nodes:
      donorName = Internal.getValue(i)
      r = Internal.getNodeFromType(i, 'GridConnectivityType_t')
      if r is not None:
        val = Internal.getValue(r)
        if (val == 'Overset' and zoneName == donorName):
          (parent, d) = Internal.getParentOfNode(z, i)
          del parent[2][d]
  return None

# enleve les BCOverlap avec domaine attache
def _deleteBCOverlapWithDonorZone__(a):
  zones = Internal.getZones(a)
  for z in zones:
    zoneName = z[0]
    nodes = Internal.getNodesFromType1(z, 'GridConnectivity_t')
    for i in nodes:
      donorName = Internal.getValue(i)
      r = Internal.getNodeFromType(i, 'GridConnectivityType_t')
      if r is not None:
        val = Internal.getValue(r)
        if (val == 'Overset' and zoneName != donorName):
          (parent, d) = Internal.getParentOfNode(z, i)
          del parent[2][d]
  return None

# -- deleteZoneBC
# Enleve les noeuds ZoneBC d'un type donne
# si type != 'None' enleve seulement le type de BC specifie
def deleteZoneBC__(a, type='None'):
  if type == 'None':
    bc = Internal.getNodesFromName3(a, 'ZoneBC')
    if bc != []: a = rmNodes(a, 'ZoneBC')
  else: a = rmBCOfType(a, type)
  return a

def _deleteZoneBC__(a, type='None'):
  if type == 'None':
    bc = Internal.getNodesFromName3(a, 'ZoneBC')
    if bc != []: _rmNodes(a, 'ZoneBC')
  else: _rmBCOfType(a, type)
  return None

# -- deleteAllBCAndSolutions__
# Enleve les noeuds ZoneBC, ZoneGridConnectivity, FlowSolution,
# FlowSolutionCenters
def deleteAllBCAndSolutions__(a):
  a = deleteGridConnectivity__(a)
  a = deleteZoneBC__(a)
  a = deleteFlowSolutions__(a)
  return a

def _deleteAllBCAndSolutions__(a):
  _deleteGridConnectivity__(a)
  _deleteZoneBC__(a)
  _deleteFlowSolutions__(a)
  return None

def deleteChimeraInfo__(a):
  ap = Internal.copyRef(a)
  _deleteChimeraInfo__(ap)
  return ap

def _deleteChimeraInfo__(a):
  Internal._rmNodesByName(a, 'OversetHoles')
  Internal._rmNodesByName(a, 'ID_*')
  return None

# -- delete Nodes specific to solvers --
def _deleteSolverNodes__(a):
  Internal._rmNodesByName(a, ':elsA#Hybrid') # elsA
  Internal._rmNodesByName(a, '.Solver#ownData') # Fast
  return None

# -- deleteEmptyZones
# Supprime les zones ayant un nombre de noeuds ou d'elements
# nul (non structure) ou un ni, nj, nk de 0 (structure)
# IN: t: tree, base, liste de zones
# OUT: isomorphe a l'entree moins les zones ayant aucun element
def deleteEmptyZones(t):
  """Delete zones with null number of points or elements."""
  tp = Internal.copyRef(t)
  _deleteEmptyZones(tp)
  return tp

def _deleteEmptyZones(t):
  zones = Internal.getZones(t)
  for z in zones:
    dim = Internal.getZoneDim(z)
    if (dim[0] == 'Unstructured' and dim[1]*dim[2] == 0) or (dim[0] == 'Structured' and dim[1]*dim[2]*dim[3] == 0):
      (p, c) = Internal.getParentOfNode(t, z)
      if (id(p) == id(t)): del p[c] # this is a list
      else: del p[2][c]
  return None

# -- rmNodes
def rmNodes(z, name):
  """Remove nodes name from z.
  Usage: rmNodes(z, name)"""
  zc = Internal.copyRef(z)
  _rmNodes(zc, name)
  return zc

def _rmNodes(z, name):
  zn = Internal.getZones(z)
  for i in zn:
    if isinstance(name, list):
      for v in name:
        nodes = Internal.getNodesFromName2(i, v)
        for j in nodes:
          (parent, d) = Internal.getParentOfNode(i, j)
          del parent[2][d]
    else:
      nodes = Internal.getNodesFromName2(i, name)
      for j in nodes:
        (parent, d) = Internal.getParentOfNode(i, j)
        del parent[2][d]
  return None

#==============================================================================
# -- File / pyTree conversions --
#==============================================================================

# -- convertFile2PyTree
def convertFile2PyTree(fileName, format=None, nptsCurve=20, nptsLine=2,
                       density=-1., skeletonData=None, dataShape=None):
  """Read a file and return a pyTree containing file data.
  Usage: convertFile2PyTree(fileName, format, options)"""
  if format is None:
    format = Converter.convertExt2Format__(fileName); autoTry = True
  else: autoTry = False
  try: file = open(fileName, 'r')
  except: raise IOError("convertFile2PyTree: file %s not found."%fileName)
  file.close()

  if format == 'bin_cgns' or format == 'unknown':   
    format = Converter.checkFileType(fileName)
    
  if format == 'bin_cgns' or format == 'bin_adf' or format == 'bin_hdf':
    try:
      t = Converter.converter.convertFile2PyTree(fileName, format, skeletonData, dataShape)
      t = Internal.createRootNode(children=t[2])
      Internal._correctPyTree(t, level=10) # force CGNS names
      Internal._correctPyTree(t, level=2) # force unique name
      Internal._correctPyTree(t, level=7) # create familyNames
      registerAllNames(t)
      return t
    except:
      if format == 'bin_cgns' or format == 'bin_adf':
        try:
          t = Converter.converter.convertFile2PyTree(fileName, 'bin_hdf', skeletonData)
          t = Internal.createRootNode(children=t[2])
          Internal._correctPyTree(t, level=10) # force CGNS names
          Internal._correctPyTree(t, level=2) # force unique name
          Internal._correctPyTree(t, level=7) # create familyNames
          registerAllNames(t)
          return t
        except: pass
      else: # adf par defaut
        try:
          t = Converter.converter.convertFile2PyTree(fileName, 'bin_adf', skeletonData)
          t = Internal.createRootNode(children=t[2])
          Internal._correctPyTree(t, level=10) # force CGNS names
          Internal._correctPyTree(t, level=2) # force unique name
          Internal._correctPyTree(t, level=7) # create familyNames
          registerAllNames(t)
          return t
        except: pass
  elif format == 'unknown':
    try:
      t = Converter.converter.convertFile2PyTree(fileName, 'bin_adf', skeletonData)
      t = Internal.createRootNode(children=t[2])
      Internal._correctPyTree(t, level=10) # force CGNS names
      Internal._correctPyTree(t, level=2) # force unique name
      Internal._correctPyTree(t, level=7) # create familyNames
      registerAllNames(t)
      return t
    except: pass
    try:
      t = Converter.converter.convertFile2PyTree(fileName, 'bin_hdf', skeletonData)
      t = Internal.createRootNode(children=t[2])
      Internal._correctPyTree(t, level=10) # force CGNS names
      Internal._correctPyTree(t, level=2) # force unique name
      Internal._correctPyTree(t, level=7) # create familyNames
      registerAllNames(t)
      return t
    except: pass

  if format == 'bin_pickle':
    import cPickle as pickle
    print 'Reading %s (bin_pickle)...'%fileName,
    try:
      file = open(fileName, 'rb')
      a = pickle.load(file)
      file.close()
    except:
      raise TypeError("convertFile2PyTree: file %s can not be read."%fileName)
    else: print 'done.'

    if Internal.isTopTree(a): 
      registerAllNames(a)     
      return a # OK
    ret = Internal.isStdNode(a)
    if ret != -2:
      t, ntype = Internal.node2PyTree(a)
      registerAllNames(t)
      return t # OK
    # sinon, c'est un arrays (traite dans la suite)

  zn = []; bcfaces = []
  if format != 'bin_pickle': # autres formats
    if autoTry: format = None
    a = Converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine,
                                     density, zoneNames=zn, BCFaces=bcfaces)

  t = newPyTree([])
  base1 = False; base2 = False; base3 = False; base = 1; c = 0

  for i in a:
    if len(i) == 5: # Structure
      if (i[3] == 1 and i[4] == 1):
        if not base1:
          t = addBase2PyTree(t, 'Base1', 1); base1 = base; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base1][2].append(z)
      elif (i[4] == 1):
        if not base2:
          t = addBase2PyTree(t, 'Base2', 2); base2 = base; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base2][2].append(z)
      else:
        if not base3:
          t = addBase2PyTree(t, 'Base', 3); base3 = base; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base3][2].append(z)
    else: # non structure
      if (i[3] == 'BAR'):
        if not base1:
          t = addBase2PyTree(t, 'Base1', 1); base1 = base; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base1][2].append(z)
      elif (i[3] == 'TRI' or i[3] == 'QUAD'):
        if not base2:
          t = addBase2PyTree(t, 'Base2', 2); base2 = base; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base2][2].append(z)
      else:
        if not base3:
          t = addBase2PyTree(t, 'Base', 3); base3 = base; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base3][2].append(z)
    z[0] = zn[c]
    c += 1
  Internal._correctPyTree(t, level=10) # force CGNS names
  Internal._correctPyTree(t, level=2) # force unique name
  _addBCFaces(t, bcfaces)
  Internal._correctPyTree(t, level=7) # create familyNames
  registerAllNames(t)
  return t

# -- convertPyTree2File
def convertPyTree2File(t, fileName, format=None, isize=4, rsize=8,
                       endian='big', colormap=0, dataFormat='%.9e ', links=[]):
  """Write a pyTree to a file.
  Usage: convertPyTree2File(t, fileName, format, options)"""
  if t == []: print 'Warning: convertPyTree2File: nothing to write.'; return
  if format is None:
    format = Converter.convertExt2Format__(fileName)
    if format == 'unknown': format = 'bin_cgns'
  if (format == 'bin_cgns' or format == 'bin_adf' or format == 'bin_hdf'):
    tp, ntype = Internal.node2PyTree(t)
    Internal._adaptZoneNamesForSlash(tp)
    Converter.converter.convertPyTree2File(tp[2], fileName, format, links)
  elif format == 'bin_pickle':
    import cPickle as pickle
    file = open(fileName, 'wb')
    print 'Writing '+fileName+'...',
    pickle.dump(t, file, protocol=pickle.HIGHEST_PROTOCOL); file.close()
    print 'done.'
  else:
    a = center2Node(t, Internal.__FlowSolutionCenters__)
    a = getAllFields(a, 'nodes')
    a = Internal.clearList(a)
    zoneNames = getZoneNames(t, prefixByBase=False)
    BCFaces = getBCFaces(t)
    Converter.convertArrays2File(a, fileName, format, isize, rsize, endian,
                                 colormap, dataFormat, zoneNames, BCFaces)


# -- convertFile2PyTree
# def convertFile2PartialPyTree(fileName, comm, format=None, nptsCurve=20, nptsLine=2,
#                                         density=-1., skeletonData=None):
#   """Convert a file to pyTree.
#   Usage: convertFile2PyTree(fileName, format, options)"""
#   if format is None:
#     format = Converter.convertExt2Format__(fileName); autoTry = True
#   else: autoTry = False
#   try: file = open(fileName, 'r')
#   except: raise IOError("convertFile2PyTree: file %s not found."%fileName)

#   t = Converter.converter.convertFile2PartialPyTree(fileName, format, skeletonData, comm)
#   t = Internal.createRootNode(children=t[2])
#   return t

def convertFile2PyTreeFromPath(fileName, Filter, 
                               format=None, nptsCurve=20, nptsLine=2,
                               density=-1., skeletonData=None):
  """Convert a file to pyTree.
  Usage: convertFile2PyTree(fileName, format, options)"""
  if format is None:
    format = Converter.convertExt2Format__(fileName); autoTry = True
  else: autoTry = False
  try: file = open(fileName, 'r')
  except: raise IOError("convertFile2PartialPyTreeFromPath2: file %s not found."%fileName)

  t = Converter.converter.convertFile2PyTreeFromPath(fileName, format, Filter)
  # t = Internal.createRootNode(children=t[2])
  return t

def convertFile2PartialPyTreeFromPath(fileName, Filter, comm=None,
                                       format=None, nptsCurve=20, nptsLine=2,
                                       density=-1., skeletonData=None):
  """Convert a file to pyTree.
  Usage: convertFile2PyTree(fileName, format, options)"""
  if format is None:
    format = Converter.convertExt2Format__(fileName); autoTry = True
  else: autoTry = False
  try: file = open(fileName, 'r')
  except: raise IOError("convertFile2PartialPyTreeFromPath2: file %s not found."%fileName)

  t = Converter.converter.convertFile2PartialPyTree(fileName, format, skeletonData, comm, Filter)
  # t = Internal.createRootNode(children=t[2])
  return t

# -- convertPyTree2File
def convertPyTree2FilePartial(t, fileName, comm, Filter, ParallelHDF=False,
                              format=None, isize=4, rsize=8,
                              endian='big', colormap=0, dataFormat='%.9e '):
  """Convert a pyTree to a file.
  Usage: convertPyTree2File(t, fileName, format, options)"""

  # > GardeFou
  if t == []: print 'Warning: convertPyTree2File: nothing to write.'; return
  format = 'bin_hdf'

  if not ParallelHDF:
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Write Tree Data execpt Data in Filter
    SkeletonTree = Internal.copyRef(t)
    for path in Filter.keys():
      print path
      Node = Internal.getNodeFromPath(SkeletonTree, path)
      Node[1] = None
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Si MPI Mode Off (HDF Not Parralel)
    if(comm.Get_rank() == 0):
      convertPyTree2File(SkeletonTree, fileName, format)
      # > Fill up Dimension
      skeletonData = None
      Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wait for Skeleton write
    comm.barrier()

    # > Set skeletonData to Not None
    skeletonData = []

    # > Cette maniere de faire provoque de la non reproductibilite ...
    # Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > On serialize
    for lock in xrange(comm.Get_size()):
      if (lock == comm.Get_rank()):
        Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
      comm.barrier()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  else:  # > Si MPI Mode Off (HDF Not Parallel)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Write Tree Data execpt Data in Filter
    SkeletonTree = Internal.copyRef(t)
    for path in Filter.keys():
      Node = Internal.getNodeFromPath(SkeletonTree, path)
      Node[1] = None
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if (comm.Get_rank() == 0):
      convertPyTree2File(SkeletonTree, fileName, format)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # > On wait l'ecriture Skelette ...
    comm.barrier()

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Write data in filter in file (With creation of DataSpace )
    skeletonData = None  # Skeleton Data is inecfective (Normaly)
    Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# -- convertPyTree2File
def convertPyTree2FileMPI(t, fileName, comm, SkeletonTree, ParallelHDF=False,
                             format=None, isize=4, rsize=8,
                             endian='big', colormap=0, dataFormat='%.9e '):
  """Convert a pyTree to a file.
  Usage: convertPyTree2File(t, fileName, format, options)"""

  # > GardeFou
  if t == []: print 'Warning: convertPyTree2File: nothing to write.'; return
  format = 'bin_hdf'
    
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > First step : Prepare a dictonary of Filter and a dictionnary of Property 
  #   in order to prepare for each procs the data to write ...
  Filter = dict()
  Proper = dict()
  Elmts  = dict()
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  import CGNS.PAT.cgnsutils as CGU
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  printTree(t)
  
  # 
  
  for Zone in Internal.getZones(t):
    pathsToArray = Internal.getPathsFromType(Zone, 'IndexRange_t')
    print '*'*100
    print pathsToArray
    print '*'*100
    pathsToArray  = CGU.getPathsByTypeSet(Zone, 'IndexRange_t')
    print pathsToArray
    print '*'*100
  
  # > The path who wants to effectivly write is DataArray_t
  pathsToArray = Internal.getPathsFromType(t, 'IndexRange_t')
  pathsToArray += Internal.getPathsFromType(t, 'DataArray_t')
  print pathsToArray
  
  pathsToArray  = CGU.getPathsByTypeSet(t, 'IndexRange_t')
  pathsToArray += CGU.getPathsByTypeSet(t, 'DataArray_t')
  print pathsToArray
  for path in pathsToArray:
     print path
     node = Internal.getNodeFromPath(t, path)
     
     # > The ideo is to make a global Dataspace full for the current proc and void for the other 
     NbE = list(node[1].shape)
     Beg = [0]*len(NbE); Sti = [1]*len(NbE); Blk = [1]*len(NbE)
     DataSpaceMMRY = [Beg, Sti, NbE, Blk]
     
     # > Partial filter (Voluntary not fill completely see after ...)
     Filter[path] = DataSpaceMMRY 
     
     # > You need to setup all label to rebuild the tree 
     Label    = []
     ListPath = path.split("/")[1:]
     topnode = t
     for i,l in enumerate(ListPath):
        topnode = Internal.getNodeFromName1(topnode, l)
        Label.append(topnode[3])

     # > Fill property dictionnary 
     Proper[path] = {'ProcNumber' : comm.Get_rank(), 
                     'DataType'   : node[1].dtype,
                     'DataShape'  : node[1].shape,
                     'Label'      : Label} 
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  pathsToElmts = Internal.getPathsFromType(t, 'Elements_t')
  for p in pathsToElmts:
    node = Internal.getNodeFromPath(t, p)
    Elmts[p] = [node[0], node[1], [], node[3]]
    
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Merge filter among procs ...
  ListFilter = comm.gather(Filter, root=0)
  ListProper = comm.gather(Proper, root=0)
  ListElmts  = comm.gather(Elmts , root=0)

  if(comm.Get_rank() == 0):
    DictFilterAll = dict()
    DictProperAll = dict()
    DictElmtsAll  = dict()
    for Proc in ListFilter:
       for Path in Proc.keys():
         DictFilterAll[Path] = Proc[Path]
    for Proc in ListProper:
      for Path in Proc.keys():
         DictProperAll[Path] = Proc[Path]
    for Proc in ListElmts:
      for Path in Proc.keys():
         DictElmtsAll[Path] = Proc[Path]
  else:
    DictFilterAll = None
    DictProperAll = None
    DictElmtsAll  = None
  Filter   = comm.bcast(DictFilterAll, root=0)
  Proper   = comm.bcast(DictProperAll, root=0)
  Elmts    = comm.bcast(DictElmtsAll , root=0)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Perform sort of the receive dictionnary 
  for path in Filter.keys():
    if(Proper[path]['ProcNumber'] != comm.Get_rank()):
       # > Change the global size 
       TMP    = Filter[path]
       DataSpaceGLOB = TMP[2]
       TMP[2] = [0]*len(TMP[2])
       Filter[path] = TMP+TMP+[DataSpaceGLOB]
    else:
      Filter[path] = Filter[path]+Filter[path]+[Filter[path][2]]
    # print path, Filter[path]
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > A partir de l'arbre on recreer tout les noeuds 
  for path in Filter.keys():
    # Node = Internal.getNodeFromPath(t, path)
    ListPath = path.split("/")[1:]
    topnode = t
    for i,l in enumerate(ListPath):
      if(Internal.getNodeFromName1(topnode, l) is None):
        if(i == len(ListPath)-1):
          shape   = [0]*len(list(Proper[path]['DataShape']))
          lvalue  = numpy.empty(shape, dtype=Proper[path]['DataType'])
          # topnode = Internal.createUniqueChild(topnode, l, 'DataArray_t', value=lvalue)
          topnode = Internal.createUniqueChild(topnode, l, Proper[path]['Label'][i], value=lvalue)
        else:
          topnode = Internal.createUniqueChild(topnode, l, Proper[path]['Label'][i])
      else:
        topnode = Internal.getNodeFromName1(topnode, l)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Make a proper Skeleton tree wit
  SkeletonTree2 = Internal.copyRef(t)
  for path in Filter.keys():
    Node = Internal.getNodeFromPath(SkeletonTree2, path)
    Node[1] = None
    
  if(comm.Get_rank() == 0):
    for path in Elmts.keys():
      ListPath = path.split("/")[1:-1]
      EndPath  = path.split("/")[-1]
      topnode  = SkeletonTree2
      for l in ListPath:
        print l
        topnode = Internal.getNodeFromName1(topnode, l)
      print path
      if(Internal.getNodeFromName1(topnode, EndPath) == None):
        Internal._addChild(topnode, Elmts[path])
      else:
        Node = Internal.getNodeFromName1(topnode, EndPath)
        Node[1] = Elmts[path][1]    
    
  for Zone in Internal.getZones(SkeletonTree2):
     Zone[1] = Internal.getNodeFromName2(SkeletonTree, Zone[0])[1]
    
  SkeletonTree = Internal.merge([SkeletonTree2, SkeletonTree])
  
  for path in Filter.keys():
    Node = Internal.getNodeFromPath(SkeletonTree, path)
    Node[1] = None
  
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if(comm.Get_rank() == 0):
    printTree(SkeletonTree)
    
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Write skeleton 
  if(comm.Get_rank() == 0):
    convertPyTree2File(SkeletonTree, fileName, format)
  comm.barrier()
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Write data in filter in file (With creation of DataSpace )
  skeletonData = None  # Skeleton Data is inecfective (Normaly)
  Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  return None

#==============================================================================
# -- Array / PyTree conversions --
#==============================================================================

# -- convertArray2ZoneNode
# Cree une nouvelle zone a partir d'arrays
# IN: A: liste d'arrays. Contient un array ou 2 arrays (dans ce cas, le
# premier array est en noeuds et le deuxieme en centres).
# OUT: noeud zone
def convertArrays2ZoneNode(zoneName, A):
  """Convert arrays to a zone node.
  Usage: convertArrays2ZoneNode(zoneName, A)"""
  if len(A) == 1:
    z = Internal.createZoneNode(getZoneName(zoneName), A[0], [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
  elif len(A) == 2:
    z = Internal.createZoneNode(getZoneName(zoneName), A[0], A[1],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
  else:
    raise TypeError("convertArrays2ZoneNode: A must be [a] or [node, center] when creating %s."%zoneName)
  return z

# -- getField
# Retourne des arrays correspondant a un nom de DataArray_t dans t
# Retourne une liste d'arrays par zone. Si une zone ne possede pas ce champ,
# retourne [] pour cette zone.
# Ex: getField('CoordinateX', t), retourne le champ CoordinateX pour toutes
# les zones de t
# Ex: getField('centers:Density', z), retourne le champ de Density pour la
# zone z.
def getField(name, t):
  nodes = Internal.getZones(t)
  arrays = []
  loc = 'nodes'
  spl = name.split(':')
  if len(spl) != 1:
    if spl[0] == 'centers': loc = 'centers'
    name = spl[1]

  for z in nodes:
    dim = Internal.getZoneDim(z)
    connects = []
    if dim[0] == 'Structured': np = dim[1]*dim[2]*dim[3]
    else:
      np = dim[1]
      connects = Internal.getElementNodes(z)

    info = z[2]
    a = None
    if loc == 'nodes': # regarde les containeurs GridCoordinates et Nodes
      for i in info:
        if (i[0] == Internal.__GridCoordinates__ or i[0] == Internal.__FlowSolutionNodes__):
          info2 = i[2]
          for j in info2:
            if (j[0] == name and j[3] == 'DataArray_t'):
              r = Internal.convertDataNode2Array(j, dim, connects, 0)
              a = r[1]
              if a is not None: break

    elif loc == 'centers':
      for i in info:
        if i[0] == Internal.__FlowSolutionCenters__:
          info2 = i[2]
          for j in info2:
            if (j[0] == name and j[3] == 'DataArray_t'):
              r = Internal.convertDataNode2Array(j, dim, connects, 1)
              a = r[1]
              if a is not None: break

    if a is not None: arrays.append(a)
    else: arrays.append([])
  return arrays

# -- getFields
# Retourne les arrays a partir d'un nom de conteneur a partir de n'importe quel
# noeud d'un arbre
# Retourne une liste de arrays
# IN: t: n'importe quel noeud de l'arbre
# IN: name: GridCoordinates, FlowSolution, FlowSolution#Centers (conteneur)
# OUT: arrays: solution
# OUT: peut contenir des arrays vides ([])
# Rem: un conteneur ne contient que des champs homogenes en localisation
def getFields(name, t, api=1):
  zones = Internal.getZones(t)
  arrays = []
  if isinstance(name, list): names = name
  else: names = [name]
  for z in zones:
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Structured':
      np = dim[1]*dim[2]*dim[3]
      connects = []
    else:
      np = dim[1]
      connects = Internal.getElementNodes(z)

    info = z[2]; out = []; loc = 0
    for i in info:
      if i[0] in names:
        locf = Internal.getNodeFromType1(i, 'GridLocation_t')
        if (locf is not None and Internal.getValue(locf) == 'CellCenter'): loc = 1
        for j in i[2]:
          if j[3] == 'DataArray_t': out.append(j)

    if out != []: 
      if api==1: array = Internal.convertDataNodes2Array(out, dim, connects, loc)
      else: array = Internal.convertDataNodes2Array2(out, dim, connects, loc)
    else: array = []
    arrays.append(array)
  return arrays

# -- getAllFields
# Retourne les arrays correspondant a une localisation
# si loc='nodes', retourne les champs de
# __GridCoordinates__ + __FlowSolutionNodes__
# si loc='centers', retourne les champs de __FlowSolutionCenters__
# OUT: peut contenir des arrays vides ([])
def getAllFields(t, loc):
  nodes = Internal.getZones(t)
  result = []
  for z in nodes:
    if loc == 'nodes':
        f = getFields([Internal.__GridCoordinates__,Internal.__FlowSolutionNodes__], z)[0]
    else:
        f = getFields(Internal.__FlowSolutionCenters__, z)[0]
    result.append(f)
  return result

# -- setPartialFields
# Set field values in zones at given indices
# IN: t: pyTree a modifier
# IN: arrays: liste des arrays des champs a modifier.
# Chaque array correspond a une zone de t a modifier. Il doit
# contenir le meme nbre de champs que chaque zone de t.
# IN: listIndices: liste des numpys des indices a modifier pour chaque zone
# sous forme de numpys d'entiers.
# NB: un array de la liste arrays peut etre vide (=[])
def setPartialFields(t, arrays, listIndices, loc='nodes'):
  """Set some field values for given indices."""
  tp = Internal.copyRef(t)
  _setPartialFields(tp, arrays, listIndices, loc)
  return tp

def _setPartialFields(t, arrays, listIndices, loc='nodes'):
  if loc == 'nodes': locI = 0
  else: locI = 1
  nodes = Internal.getZones(t)
  nzones = len(nodes)
  for c in xrange(nzones):
    a = arrays[c]; indices = listIndices[c]
    z = nodes[c] # zone
    Converter.converter.setPartialFieldsPT(z, a, indices, locI, 
                                           Internal.__GridCoordinates__, 
                                           Internal.__FlowSolutionNodes__, 
                                           Internal.__FlowSolutionCenters__)
  return None

# --setPartialFields1
# Set field values in zones at given indices
# IN: t: pyTree a modifier
# IN: listFields: liste pour chaque zone de t des numpys des champs.
# Chaque zone  est representee par une liste de numpys (1 par champ).
# IN: listIndices: liste des numpys des indices a modifier pour chaque zone
# sous forme de numpys d'entiers.
# NB: un numpy de la liste des numpys peut etre vide (=[])
def setPartialFields1(t, listFields, listIndices, loc='nodes'):
  tp = Internal.copyRef(t)
  _setPartialFields1(tp, listFields, listIndices, loc)
  return tp

def _setPartialFields1(t, listFields, listIndices, loc='nodes'):
  if loc == 'nodes': locI = 0
  else: locI = 1
  zones = Internal.getZones(t)
  nzones = len(zones)
  for c in xrange(nzones):
    Converter.converter._setPartialFields(zones[c], listFields[c], listIndices[c], locI, 
                                          Internal.__GridCoordinates__, 
                                          Internal.__FlowSolutionNodes__, 
                                          Internal.__FlowSolutionCenters__)
  return None

# -- setFields
# Set fields a partir de n'importe quel noeud de l'arbre et
# d'une liste d'arrays dans les conteneurs standards:
# si loc=nodes, conteneurs=__GridCoordinates__, __FlowSolutionNodes__
# si loc=centers, conteneurs=__FlowSolutionCenters__
# IN: arrays: liste d'array qui seront mis dans les noeuds zones
# IN: doit contenir des arrays vides si certains champs sont vides.
# IN: t: noeud designant une base, un arbre ou une zone
# IN: loc='nodes' ou 'centers': localisation des arrays
# IN: writeDim: True: force l'ecriture de la dimension/connect du noeud zone
def setFields(arrays, t, loc, writeDim=True):
  # Recherche de tous les noeuds zones a partir de t
  nodes = Internal.getZones(t)

  # Verification de la coherence
  if not isinstance(arrays[0],list):
    raise TypeError("setFields: arrays should be a list of array.")
  if len(nodes) != len(arrays):
    raise ValueError("setFields: more zones in tree than in arrays.")

  for c in xrange(len(arrays)):
    a = arrays[c]
    z = nodes[c] # zone
    info = z[2]

    if (writeDim == True and loc == 'nodes' and a != []):
      d = Internal.array2PyTreeDim(a)
      cellDim = d.shape[0]
    else:
      cellDim = z[1].shape[0]

    # Remplace les noeuds contenant les variables
    if a == []: vars = []
    else:
      vars = a[0].split(",")
      # un array * ne peut pas etre mis en nodes
      if (loc == 'nodes' and len(a) == 4):
        elt = a[3]
        if (elt[len(elt)-1] == '*'):
          print 'Warning: setFields: %s array is not set.'%elt
          vars = []
    p = 0
    for v in vars:
      renamed = 0 # si le nom est change en nom CGNS = 1
      if Internal.name2CGNS.has_key(v): renamed = 1; variable = Internal.name2CGNS[v]
      else: variable = v
      if ((variable == 'CoordinateX' or variable == 'CoordinateY'
          or variable == 'CoordinateZ') and loc == 'nodes'):
        coordNode = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if coordNode is None:
          info = [Internal.__GridCoordinates__, None, [], 'GridCoordinates_t']
          z[2].append(info)
        else: info = coordNode
        l = Internal.getNodesFromName(info, variable)
        if l != []:
            l[0][1] = Internal.createDataNode(variable, a, p, cellDim)[1]
        else:
          node = Internal.createDataNode(variable, a, p, cellDim)
          info[2].append(node)
          if renamed == 1: Internal._rmNodesByName(info, v)

      else: # FlowSolution
        if loc == 'nodes':
          flowNode = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        else:
          flowNode = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if flowNode is None:
          if loc == 'nodes':
            info = [Internal.__FlowSolutionNodes__, None, [], 'FlowSolution_t']
          else:
            info = [Internal.__FlowSolutionCenters__, None, [], 'FlowSolution_t']
            v = numpy.fromstring('CellCenter', 'c')
            info[2].append(['GridLocation', v, [], 'GridLocation_t'])
          z[2].append(info)
        else:
          info = flowNode
        l = Internal.getNodesFromName(info, variable)
        if loc == 'nodes':
          node = Internal.createDataNode(variable, a, p, cellDim)
        else:
          node = Internal.createDataNode(variable, a, p, cellDim)
        if l != []: l[0][1] = node[1]
        else:
          info[2].append(node)
          if renamed == 1: Internal._rmNodesByName(info, v)

      p += 1

    # update les dimensions si necessaire
    if (writeDim == True and loc == 'nodes' and vars != []):
      z[1] = Internal.array2PyTreeDim(a)
      if len(a) == 5: # structure
        typeNodes = Internal.getNodesFromType1(z, 'ZoneType_t')
        val = Internal.getValue(typeNodes[0])
        if (val != 'Structured'):
          # Supprimer GridElementsNode
          GENodes = Internal.getElementNodes(z)
          for n in GENodes:
            p, r = Internal.getParentOfNode(z, n)
            del p[2][r]
          typeNodes[0][1] = numpy.fromstring('Structured', 'c')

      elif len(a) == 4: # non structure
        typeNodes = Internal.getNodesFromType1(z, 'ZoneType_t')
        val = Internal.getValue(typeNodes[0])
        if val == 'Structured':
          v = numpy.fromstring('Unstructured', 'c')
          typeNodes[0][1] = v
          Internal.setElementConnectivity(z, a)
        else:
          Internal.setElementConnectivity(z, a)
  return t

# -- getNumpyArrays
# Retourne une reference sur le tableau numpy correspondant a un champ
# IN: t: noeud zone, base, tree
# IN: name: nom du champ a extraire (pas de conteneur)
# OUT: liste des numpy arrays (references)
def getNumpyArrays(t, name):
  zones = Internal.getZones(t)
  result = []
  for z in zones:
    sol = Internal.getNodeFromName2(z, name)
    if sol is not None: result.append(sol[1])
    else: result.append([])
  return result

# -- setNumpyArrays
# Remplace une reference sur le tableau numpy correspondant a un champ
# IN: t: noeud zone, base, tree
# IN: name: nom du champ a extraire
# IN: liste des numpy arrays (references)
# OUT: t modifie
def setNumpyArrays(t, name, arrays):
  zones = Internal.getZones(t)
  if (len(arrays) != len(zones)):
    print 'Error: setNumpyArrays: not enough arrays.'; return
  i = 0
  for z in zones:
    sol = Internal.getNodeFromName2(z, name)
    if sol is not None: sol[1] = arrays[i]
    i += 1
  return

# -- ownNumpyArrays
# For numpys allocated with KCore empty (OWNDATA=False), reown it
def _ownNumpyArrays(t):
  n = t[1]
  if (isinstance(n, numpy.ndarray) and n.flags['OWNDATA'] == False):
    b = numpy.copy(n); t[1] = b
  for i in t[2]: _ownNumpyArrays(i)
  return None

# -- convertPyTree2Array
# Convert a python tree node value to an array if possible
def convertPyTree2Array(path, tree):
    """Convert a Python tree node to an array.
    Usage: convertPyTree2Array(path, tree)"""
    p = path.split('/')
    zone = p[0]+'/'+p[1]
    z = Internal.getNodeFromPath(tree, zone)
    if z is None:
      raise ValueError("convertPyTree2Array: zone %s not found."%zone)

    dim = Internal.getZoneDim(z)
    connects = Internal.getElementNodes(z)

    a = Internal.getNodeFromPath(tree, path)

    if a is None:
      raise ValueError("convertPyTree2Array: path %s not found."%path)

    if a[3] == 'DataArray_t':
       array = Internal.convertDataNode2Array(a, dim, connects)[1]
       return array

    elif a[3] == 'IndexRange_t':
      b = numpy.zeros((1,6))
      b.flat = a[1].flat
      return ['index', b, 6, 1, 1]

    elif a[3] == 'GridCoordinates_t':
      info = a[2]; out = []
      for i in info:
        if i[3] == 'DataArray_t':
          a = Internal.convertDataNode2Array(i, dim, connects, 0)[1]
          if a is not None: out.append(a)
      if out != []:
        array = Converter.addVars(out); return array
      else: return None
    elif a[3] == 'FlowSolution_t':
      loc = 0
      locf = Internal.getNodeFromType1(a, 'GridLocation_t')
      if (locf is not None and Internal.getValue(locf) == 'CellCenter'): loc = 1
      info = a[2]; out = []
      for i in info:
        if i[3] == 'DataArray_t':
          a = Internal.convertDataNode2Array(i, dim, connects, loc)[1]
          if a is not None: out.append(a)
      if (out != []):
        array = Converter.addVars(out); return array
      else: return None

    elif a[3] == 'Zone_t':
      if dim[0] == 'Structured':
        np = dim[1]*dim[2]*dim[3]
      else:
        np = dim[1]
      info = a[2]; out = []; out2 = []
      for i in info:
        if (i[3] == 'GridCoordinates_t'):
          info2 = i[2]
          for j in info2:
            if (j[3] == 'DataArray_t'):
              a = Internal.convertDataNode2Array(j, dim, connects, 0)[1]
              if a is not None:
                if (a[1].shape[1] == np): out.append(a)
                else: out2.append(a)

        if i[3] == 'FlowSolution_t':
          loc = 0
          locf = Internal.getNodeFromType1(i, 'GridLocation_t')
          if (locf is not None and Internal.getValue(locf) == 'CellCenter'): loc = 1
          info2 = i[2]
          for j in info2:
            if j[3] == 'DataArray_t':
              a = Internal.convertDataNode2Array(j, dim, connects, loc)[1]
              if a is not None:
                if (a[1].shape[1] == np): out.append(a)
                else: out2.append(a)

      if (out != [] and out2 != []):
        array = Converter.addVars(out)
        array2 = Converter.addVars(out2)
        print "Warning: convertPyTree2Array: only node field are in array."
        return array # return array, array2
      elif out != []:
        array = Converter.addVars(out); return array
      elif out2 != []:
        array2 = Converter.addVars(out2); return array2
      else: return None

    else:
      raise ValueError("convertPyTree2Array: in given path, no data node was found.")

# -- makeAllNumpy
# Remplace tous les noeuds 'valeur' (int ou float) par des noeuds numpy de
# une valeur
def makeAllNumpy(t):
  """Make a pyTree all numpy compliant.
  Usage: makeAllNumpy(a)"""
  a = Internal.copyRef(t)
  _makeAllNumpy(a)
  return a

def _makeAllNumpy(a):
  if Internal.isTopTree(a):
    for c in a[1:]: makeAllNumpy__(c)
  else: makeAllNumpy__(a)
  return None

def makeAllNumpy__(node):
  for c in node[2]:
    if isinstance(c[1], str):
      c[1] = numpy.array(c[1], dtype='c')
    if isinstance(c[1], int):
      c[1] = numpy.array(c[1], dtype=numpy.int32)
    if isinstance(c[1], float):
      c[1] = numpy.array(c[1], dtype=numpy.float64)
    makeAllNumpy__(c)

# -- makeNumpyStringString
# Remplace tous les noeuds numpy strings (tableau de lettres) par
# des noeuds string
def makeNumpyStringString(t):
  """Make a pyTree all numpy string to strings.
  Usage: makeNumpyStringString(a)"""
  a = Internal.copyRef(t)
  if Internal.isTopTree(a):
    for c in a[1:]: makeNumpyStringString__(c)
  else: makeNumpyStringString__(a)
  return a

def makeNumpyStringString__(node):
  for c in node[2]:
    c[1] = getValue(c)
    makeNumpyStringString__(c)

#==============================================================================
# -- Traitements generiques --
#==============================================================================

# -- TZGC
# Traitement agissant sur les coords, effectue par zones.
# Prend le champ GC (noeuds), applique F, remet le resultat dans un
# conteneur suivant locout.
def TZGC(t, locout, F, *args):
  tp = Internal.copyRef(t)
  _TZGC(tp, locout, F, *args)
  return tp

def _TZGC(t, locout, F, *args):
  zones = Internal.getZones(t)
  for z in zones:
    coord = getFields(Internal.__GridCoordinates__, z)[0]
    if coord != []:
      coord = F(coord, *args)
      setFields([coord], z, locout)
  return None

def TZGC2(t, locout, F, *args):
  tp = Internal.copyRef(t)
  zones = Internal.getZones(t)
  for z in zones:
    coord = getFields(Internal.__GridCoordinates__, z, api=2)[0]
    if coord != []:
      coord = F(coord, *args) # copie
      setFields([coord], z, locout)
  return tp

def _TZGC2(t, locout, F, *args):
  zones = Internal.getZones(t)
  for z in zones:
    coord = getFields(Internal.__GridCoordinates__, z, api=2)[0]
    if coord != []:
      F(coord, *args) # in place array2
  return None

# -- TZGF
# Traitement agissant sur le conteneur fieldName, effectue par zones.
# Prend les champs definis dans le conteneur fieldName, applique F, met
# le resultat dans un conteneur suivant locout.
def TZGF(t, locout, fieldName, F, *args):
  tp = Internal.copyRef(t)
  _TZGF(tp, locout, fieldName, F, *args)
  return tp

def _TZGF(t, locout, fieldName, F, *args):
  zones = Internal.getZones(t)
  for z in zones:
    fa = getFields(fieldName, z)[0]
    if fa != []:
      if F is not None: fa = F(fa, *args)
      setFields([fa], z, locout)
  return None

# -- TZA
# Traitement requerant les coord + tous les champs, effectue par zones.
# Prend les champs situes dans locin, effectue le traitement F sur ces champs,
# Remet le resultat dans locout.
# IN: locin: nodes, centers, both
# IN: locout: nodes, centers, both
# sont permis nodes / nodes, nodes / centers, centers / nodes,
# centers / centers, both / both (dans ce cas => nodes/nodes + centers/centers)
# IN: F: fonction de traitement pour les noeuds (possible None)
# IN: Fc: fonction de traitement pour les centres (possible None)
# IN: args: argument de la fonction F
# Si locin=both, args doit contenir les arguments de F pour les noeuds
# suivis des arguments de Fc pour les centres. Les champs en centres ne
# contiennent pas les coordonnees. La fonction F est appelee une
# fois pour les noeuds, Fc une fois pour les centres.
def TZA(t, locin, locout, F, Fc, *args):
  tp = Internal.copyRef(t)
  _TZA(tp, locin, locout, F, Fc, *args)
  return tp

def _TZA(t, locin, locout, F, Fc, *args):
  zones = Internal.getZones(t)
  for z in zones:
    if locin == 'nodes':
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
      if (fc != [] and fa != []):
        f = Converter.addVars([fc, fa])
        fp = F(f, *args)
        setFields([fp], z, locout)
      elif fa != []:
        fp = F(fa, *args)
        setFields([fp], z, locout)
      elif fc != []:
        fp = F(fc, *args)
        setFields([fp], z, locout)
    elif locin == 'centers':
      fa = getFields(Internal.__FlowSolutionCenters__, z)[0]
      if fa != []:
        fp = Fc(fa, *args)
        setFields([fp], z, locout)
    else: # both
      # Dans ce cas, on suppose que F ne change pas la localisation
      l = len(args);
      args1 = args[0:l/2]; args2 = args[l/2:]
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
      fb = getFields(Internal.__FlowSolutionCenters__, z)[0]
      if (fc != [] and fa != []):
        f = Converter.addVars([fc, fa])
        fp = F(f, *args1)
        setFields([fp], z, 'nodes')
      elif (fa != []):
        fp = F(fa, *args1)
        setFields([fp], z, 'nodes')
      elif (fc != []):
        fp = F(fc, *args1)
        setFields([fp], z, 'nodes')
      if (fb != []):
        if Fc is not None: fb = Fc(fb, *args2)
        setFields([fb], z, 'centers')
  return None

# -- TZAGC
# Traitement effectue pour tous les champs + coord. memes pour les centres.
# Dans ce cas, on reconstruit un maillage en centres qui sera associe
# au champ en centres.
# IN: t: arbre a traiter
# IN: locin: nodes, centers, both
# IN: locout: nodes, centers, both
# IN: F: fonction a appliquer pour les noeuds
# IN: Fc: fonction a appliquer pour les centres
# IN: args: arguments de F pour les noeuds + arguments de Fc pour les centres
# dans le cas both.
# La fonction F est appelee une fois pour les noeuds, Fc une fois pour les
# centres.
def TZAGC(t, locin, locout, F, Fc, *args):
  tp = Internal.copyRef(t)
  _TZAGC(tp, locin, locout, F, Fc, *args)
  return tp

def _TZAGC(t, locin, locout, F, Fc, *args):
  zones = Internal.getZones(t)
  for z in zones:
    if (locin == 'nodes'):
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
      if fc != [] and fa != []:
        f = Converter.addVars([fc, fa])
        fp = F(f, *args)
        setFields([fp], z, locout)
      elif fa != []:
        fp = F(fa, *args)
        setFields([fp], z, locout)
      elif fc != []:
        fp = Fc(fc, *args)
        setFields([fp], z, locout)
    elif locin == 'centers':
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionCenters__, z)[0]

      if (fc != [] and fa != []):
        fc2 = Converter.node2Center(fc)
        f = Converter.addVars([fc2, fa])
        fp = Fc(f, *args)
        st = fp[0].split(',')
        vars = []
        for i in st:
          if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            vars.append(i)
        if vars != []:
          fp = Converter.extractVars(fp, vars)
          setFields([fp], z, locout)
      elif fc != []:
        fc2 = Converter.node2Center(fc)
        fp = Fc(fc2, *args)
        st = fp[0].split(',')
        vars = []
        for i in st:
          if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            vars.append(i)
        if vars != []:
          fp = Converter.extractVars(fp, vars)
          setFields([fp], z, locout)
      elif fa != []:
        fp = Fc(fa, *args)
        setFields([fp], z, locout)
    else: # both
      l = len(args)
      args1 = args[0:l/2]; args2 = args[l/2:]
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
      fb = getFields(Internal.__FlowSolutionCenters__, z)[0]
      
      if (fc != [] and fa != []):
        f = Converter.addVars([fc, fa])
        fp = F(f, *args1)
        setFields([fp], z, 'nodes')
      elif fa != []:
        fp = F(fa, *args1)
        setFields([fp], z, 'nodes')
      elif fc != []:
        fp = F(fc, *args1)
        setFields([fp], z, 'nodes')

      if (fc != [] and fb != []):
        fc2 = Converter.node2Center(fc)
        f = Converter.addVars([fc2, fb])
        fp = Fc(f, *args2)
        st = fp[0].split(',')
        vars = []
        for i in st:
          if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            vars.append(i)
        if vars != []:
          fp = Converter.extractVars(fp, vars)
          setFields([fp], z, 'centers')
      elif fb != []:
        fp = Fc(fb, *args2)
        setFields([fp], z, 'centers')
  return None

# -- TZANC
# Traitement effectue pour tous les champs + coord. memes pour les centres.
# Dans ce cas, on passe le champ aux centres en noeuds, on applique F,
# et on repasse le champ en centres
# IN: t: arbre a traiter
# IN: locin: nodes, centers, both
# IN: locout: nodes, centers, both
# IN: F: fonction a appliquer pour les noeuds
# IN: Fc: fonction a appliquer pour les centres
# IN: args: arguments de F pour les noeuds + arguments de Fc pour les centres
# dans le cas both.
# La fonction F est appelee une fois pour les noeuds, Fc une fois pour les
# centres.
def TZANC(t, locin, locout, F, Fc, *args):
  tp = Internal.copyRef(t)
  _TZANC(tp, locin, locout, F, Fc, *args)
  return tp

def _TZANC(t, locin, locout, F, Fc, *args):
  zones = Internal.getZones(t)
  l = len(args)
  args1 = args[0:l/2]; args2 = args[l/2:]
  for z in zones:
    if (locin == 'nodes'):
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
      if (fc != [] and fa != []):
        f = Converter.addVars([fc, fa])
        fp = F(f, *args)
        setFields([fp], z, locout)
      elif fa != []:
        fp = F(fa, *args)
        setFields([fp], z, locout)
      elif fc != []:
        fp = F(fc, *args)
        setFields([fp], z, locout)
    elif (locin == 'centers'):
      zp = Internal.copyRef(z)
      _deleteFlowSolutions__(zp, 'nodes')
      zp = center2Node(zp, Internal.__FlowSolutionCenters__)
      fa = getFields(Internal.__FlowSolutionNodes__, zp)[0]
      if fa != []:
        fa = Fc(fa, *args)
        if (locout == 'nodes'): setFields([fa], z, 'nodes')
        else:
          zp = node2Center(zp, Internal.__FlowSolutionNodes__)
          fa = getFields(Internal.__FlowSolutionCenters__, zp)[0]
          setFields([fa], z, 'centers')

    else: # both
      l = len(args)
      args1 = args[0:l/2]; args2 = args[l/2:]
      fc = getFields(Internal.__GridCoordinates__, z)[0]
      fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
      zp = Internal.copyRef(z)
      _deleteFlowSolutions__(zp, 'nodes')
      zp = center2Node(zp, Internal.__FlowSolutionCenters__)
      fb = getFields(Internal.__FlowSolutionNodes__, zp)[0]
      if fc != [] and fa != []:
        f = Converter.addVars([fc, fa])
        fp = F(f, *args1)
        setFields([fp], z, 'nodes')
      elif fa != []:
        fp = F(fa, *args1)
        setFields([fp], z, 'nodes')
      elif fc != []:
        fp = F(fc, *args1)
        setFields([fp], z, 'nodes')

      if fb != []:
        if Fc is not None: fb = Fc(fb, *args2)
        setFields([fb], zp, 'nodes')
        zp = node2Center(zp, Internal.__FlowSolutionNodes__)
        fa = getFields(Internal.__FlowSolutionCenters__, zp)[0]
        setFields([fa], z, 'centers')
  return None

# -- TLAGC
# Traitement effectue PAR LOT pour tous les champs + coord. memes pour
# les centres.
# Dans ce cas, on reconstruit un maillage en centres qui sera associe
# au champ en centres.
# IN: t: arbre a traiter
# IN: F: fonction a appliquer
# IN: args: arguments de F pour les noeuds + arguments de F pour les centres
# OUT: le resultat du traitement par lot (entier, reels, arrays)
# Si sortie=arrays de Converter
#  1- c'est a l'appelant de retourner une zone/liste de zones/arbre
#  2- traitement des centres: l'appelant doit enlever les coordonnees des centres
def TLAGC(t, F, *args):
  tp = Internal.copyRef(t)
  zones = Internal.getZones(tp)
  l = len(args)
  args1 = args[0:l/2]; args2 = args[l/2:]
  allfc = []; allfa = []; allfb = []
  nzones = len(zones)
  if nzones == 0: return tp

  for z in zones:
    fc = getFields(Internal.__GridCoordinates__, z)[0]
    fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
    fb = getFields(Internal.__FlowSolutionCenters__, z)[0]
    allfc.append(fc); allfa.append(fa); allfb.append(fb)

  # Application par lot des noeuds
  allf = []
  for i in xrange(nzones):
    if (allfc[i] != [] and allfa[i] != []):
      allf.append(Converter.addVars([allfc[i], allfa[i]]))
    elif (allfa[i] != []): allf.append(allfa[i])
    elif (allfc[i] != []): allf.append(allfc[i])
  if allf != []: allf = F(allf, *args1)
  allfa = allf # ref

  # Application par lot des centres
  allf = []
  allfc = Converter.node2Center(allfc)
  for i in xrange(nzones):
    if (allfc[i] != [] and allfb[i] != []):
      allf.append(Converter.addVars([allfc[i],allfb[i]]))
    elif (allfb[i] != []): allf.append(allfb[i])

  if (allf != []): allf = F(allf, *args2)
  allfb = allf # ref

  if allfa != [] and allfb != []: return [allfa,allfb]
  elif allfa != []: return [allfa]
  elif allfb != []: return [allfb]
  else: return []

#==============================================================================
# -- Fields / Vars management --
#==============================================================================

# -- addVars: ajoute des variables dans un pyTree
def addVars(t, vars):
  """Add variables to a pyTree.
  Usage: a = addVars(t, vars)"""
  tp = Internal.copyRef(t)
  _addVars(tp, vars)
  return tp

def _addVars(t, vars):
  zones = Internal.getZones(t)
  if isinstance(vars, list):
    for var in vars:
      loc = 'nodes'; v = var.split(':')
      if len(v) > 1: var = v[1]; loc = v[0]
      for z in zones:
        variable = Internal.getCGNSName(var)
        found = 0
        if loc == 'centers':
          node = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
          if node is not None:
            nodeofname = Internal.getNodeFromName1(node, variable)
            if nodeofname is not None: found = 1
        if loc == 'nodes':
          node = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
          if node is not None:
            nodeofname = Internal.getNodeFromName1(node, variable)
            if nodeofname is not None: found = 1
          node = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
          if node is not None:
            nodeofname = Internal.getNodeFromName1(node, variable)
            if nodeofname is not None: found = 1
        if found == 0:
          dim = Internal.getZoneDim(z)
          _addAVar__(z, dim, '%s:%s'%(loc,variable))
  else:
    loc = 'nodes'; v = vars.split(':')
    if len(v) > 1: vars = v[1]; loc = v[0]
    variable = Internal.getCGNSName(vars)
    for z in zones:
      found = 0
      if loc == 'centers':
        node = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if node is not None:
          nodeofname = Internal.getNodeFromName1(node, variable)
          if nodeofname is not None: found = 1
      if loc == 'nodes':
        node = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if node is not None:
          nodeofname = Internal.getNodeFromName1(node, variable)
          if nodeofname is not None: found = 1
        node = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if node is not None:
          nodeofname = Internal.getNodeFromName1(node, variable)
          if nodeofname is not None: found = 1
      if found == 0:
        dim = Internal.getZoneDim(z)
        _addAVar__(z, dim, '%s:%s'%(loc,variable))
  return None

def _addAVar__(z, dim, var):
  vref = getVarNames(z)
  for i in vref:
    if i == var: return # variable deja existante
  loc = 'nodes'
  v = var.split(':')
  if len(v) > 1: var = v[1]; loc = v[0]
  else: var = v[0]

  if (dim[0] == 'Structured'):
    if (loc == 'nodes'):
      a = Converter.array(var, dim[1], dim[2], dim[3])
    if (loc == 'centers'):
      a = Converter.array(var, max(dim[1]-1,1), max(dim[2]-1,1), max(dim[3]-1,1))
  else:
    if dim[3] != 'NGON':
      if loc == 'nodes': a = Converter.array(var, dim[1], dim[2], dim[3])
      elif loc == 'centers': a = Converter.array(var, dim[1], dim[2], dim[3]+'*')
    else: # NGON
      if loc == 'nodes':
        a = getAllFields(z, 'nodes')[0]
        a = Converter.addVars(a, var)
      else:
        a = getAllFields(z, 'centers')[0]
        if a == []:
          a = getFields(Internal.__GridCoordinates__,z)[0]
          a = Converter.node2Center(a)
        a = Converter.addVars(a, var)
        a = Converter.extractVars(a, [var])
  setFields([a], z, loc, False)
  return None

# -- initVars: initialise une variable
def initVars(t, varNameString, val1=[], val2=[]):
  """Init variables defined by varNameString.
  Usage: a = initVars(array, varNameString, val)
  or
  Usage: a = initVars(array, varNameString, F, [strings])"""
  tp = Internal.copyRef(t)
  _initVars(tp, varNameString, val1, val2)
  return tp

def _initVars(t, varNameString, val1=[], val2=[]):
  # Ajoute la variable si necessaire
  s = varNameString.split('=')
  if len(s) == 1: # pas formule
    loc = 'nodes'
    v = varNameString.split(':')
    if len(v) > 1: var = v[1]; loc = v[0]
    else: var = v[0]
    _addVars(t, varNameString)
    c = 0
    for i in val2:
        v = i.split(':')
        if len(v) > 1: val2[c] = v[1]
        else: val2[c] = v[0]
        c += 1
    _TZAGC(t, loc, loc, Converter.initVars, Converter.initVars,
           var, val1, val2)
  else: # formule
    loc = 'nodes'
    v = s[0].split(':')
    if len(v) > 1: loc = v[0]; loc = loc.replace('{', '')
    _TZAGC(t, loc, loc, Converter.initVars, Converter.initVars,
           varNameString, val1, val2)
  return None

# Merge BCDataSets
def _mergeBCDataSets(t):
  zones = Internal.getZones(t)
  if zones == []: zones = [t] # must be a BC node
  for z in zones:
    bcs = Internal.getNodesFromType2(z, 'BC_t')
    for b in bcs:
      Internal._mergeBCDataSets__(z,b)
  return None


def _nullifyBCDataSetVectors(t, bndType, loc='FaceCenter', 
                             vectors=[['VelocityX','VelocityY','VelocityZ']]):
  locI = 0
  if loc == 'FaceCenter': locI = 1; FSCont = Internal.__FlowSolutionCenters__
  else: locI = 0; FSCont = Internal.__FlowSolutionNodes__

  zones = Internal.getZones(t)
  families = getFamilyBCNamesOfType(t, bndType)
  for z in zones:
    dimZ = Internal.getZoneDim(z)
    niZ = dimZ[1]; njZ = dimZ[2]; nkZ = dimZ[3]
    allbcs = Internal.getNodesFromType2(z, 'BC_t')
    bcs = Internal.getNodesFromValue(allbcs, bndType)
    bcs += getFamilyBCs(allbcs,families)    
    for bc in bcs:
      PR = Internal.getNodeFromName1(bc, 'PointRange')
      PL = Internal.getNodeFromName1(bc, 'PointList')
      np = 0
      if PR is not None: 
        win =  Internal.range2Window(PR[1])
        imin = win[0]; imax = win[1]
        jmin = win[2]; jmax = win[3]
        kmin = win[4]; kmax = win[5]
        if locI == 0: 
          di = max(1,imax-imin+1)
          dj = max(1,jmax-jmin+1)
          dk = max(1,kmax-kmin+1)
        else:
          di = max(1,imax-imin)
          dj = max(1,jmax-jmin)
          dk = max(1,kmax-kmin)
        np = di*dj*dk
      elif PL is not None:
        np = PL[1].size
      else:
        raise(ValueError,"nullifyVectorAtBCDataSet: no PointRange/PointList in BC.")
      
      datas = Internal.getBCDataSet(z,bc)
      if datas == []: # create the BCDataSet
        d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                  gridLocation=loc, parent=bc)
        d = Internal.newBCData('BCDirichlet', parent=d)
        cont, noc = Internal.getParentOfNode(z, d)
        for vect in vectors:
            vxname = vect[0]; fxInt = numpy.zeros((np),numpy.float64)
            vyname = vect[1]; fyInt = numpy.zeros((np),numpy.float64)
            vzname = vect[2]; fzInt = numpy.zeros((np),numpy.float64)
            if PR is not None:
              Converter.converter.nullifyVectorAtBCFaceStruct(z, fxInt, fyInt, fzInt,
                                                              imin, imax, jmin, jmax, kmin, kmax, locI, 
                                                              vxname, vyname, vzname,  
                                                              Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)

            elif PL is not None: 
              print "nullifyVectorAtBCDataSet: not implemented for PointList."
            Internal._createUniqueChild(d, vxname, 'DataArray_t', value=fxInt)
            Internal._createUniqueChild(d, vyname, 'DataArray_t', value=fyInt)
            Internal._createUniqueChild(d, vzname, 'DataArray_t', value=fzInt)

      else: # BCDataSet exists: add missing variables
        d = Internal.getNodeFromType2(bc,'BCData_t')
        for vect in vectors:
          vxname = vect[0]; fxInt = numpy.zeros((np),numpy.float64)
          vyname = vect[1]; fyInt = numpy.zeros((np),numpy.float64)
          vzname = vect[2]; fzInt = numpy.zeros((np),numpy.float64)
          Converter.converter.nullifyVectorAtBCFaceStruct(z, fxInt, fyInt, fzInt, 
                                                          imin, imax, jmin, jmax, kmin, kmax, locI, 
                                                          vxname, vyname, vzname,
                                                          Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)
          Internal._createUniqueChild(d, vxname, 'DataArray_t', value=fxInt)
          Internal._createUniqueChild(d, vyname, 'DataArray_t', value=fyInt)
          Internal._createUniqueChild(d, vzname, 'DataArray_t', value=fzInt)
  return None

# Create a BCDataSet for a bndType by Oth-order extrapolation from interior
# loc='FaceCenter' or 'Vertex'
# update = True: update the BCDataSet by extrapolation if field already exists
def _createBCDataSetOfType(t, bndType, loc='FaceCenter', update=True, vectors=[]):
  locI = 0
  if loc == 'FaceCenter': locI = 1; FSCont = Internal.__FlowSolutionCenters__
  else: locI = 0; FSCont = Internal.__FlowSolutionNodes__

  zones = Internal.getZones(t)
  families = getFamilyBCNamesOfType(t,bndType)
  for z in zones:
    dimZ = Internal.getZoneDim(z)
    niZ = dimZ[1]; njZ = dimZ[2]; nkZ = dimZ[3]
    allbcs = Internal.getNodesFromType2(z, 'BC_t')
    bcs = Internal.getNodesFromValue(allbcs, bndType)
    bcs += getFamilyBCs(allbcs,families)
    FSNode = Internal.getNodeFromName1(z,FSCont)
    varnames=[]
    for fs in FSNode[2]: 
      if fs[3]=='DataArray_t': varnames.append(fs[0])
    for bc in bcs:
      PR = Internal.getNodeFromName1(bc, 'PointRange')
      PL = Internal.getNodeFromName1(bc, 'PointList')
      np = 0
      if PR is not None: 
        win =  Internal.range2Window(PR[1])
        imin = win[0]; imax = win[1]
        jmin = win[2]; jmax = win[3]
        kmin = win[4]; kmax = win[5]
        if locI == 0: 
          di = max(1,imax-imin+1)
          dj = max(1,jmax-jmin+1)
          dk = max(1,kmax-kmin+1)
        else:
          di = max(1,imax-imin)
          dj = max(1,jmax-jmin)
          dk = max(1,kmax-kmin)
        np = di*dj*dk
      elif PL is not None:
        np = PL[1].size
      else:
        raise(ValueError,"createBCDataSetOfType: no PointRange/PointList in BC.")
      
      datas = Internal.getBCDataSet(z,bc)
      if datas == []: # create the BCDataSet
        d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                  gridLocation=loc, parent=bc)
        d = Internal.newBCData('BCDirichlet', parent=d)
        cont, noc = Internal.getParentOfNode(z, d)
        for fs in FSNode[2]: 
          if fs[3]=='DataArray_t': 
            varname = fs[0]
            fInt = numpy.zeros((np),numpy.float64)
            if PR is not None:
              Converter.converter.extrapInterior2BCFaceStruct(z, fInt, imin, imax, jmin, jmax, kmin, kmax, locI, varname, 
                                                              Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)

            elif PL is not None: 
              print "createBCDataSetOfType: not implemented for PointList."
            Internal._createUniqueChild(d, varname, 'DataArray_t', value=fInt)

      else: # BCDataSet exists: add missing variables
        d = Internal.getNodeFromType2(bc,'BCData_t')
        dataSetNames = []
        for dataSet in d[2]: dataSetNames.append(dataSet[0])
        for varname in varnames:
          if varname not in dataSetNames or update is True:
            fInt = numpy.zeros((np),numpy.float64)
            Converter.converter.extrapInterior2BCFaceStruct(z, fInt, imin, imax, jmin, jmax, kmin, kmax, locI, varname, 
                                                            Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)
            Internal._createUniqueChild(d, varname, 'DataArray_t', value=fInt)                      
  if bndType == 'BCWallInviscid':
    _nullifyBCDataSetVectors(t, bndType, loc=loc, vectors=[['VelocityX','VelocityY','VelocityZ']])
  elif bndType == 'BCSymmetryPlane':
    if vectors == []: 
      _nullifyBCDataSetVectors(t, bndType, loc=loc, vectors=[['VelocityX','VelocityY','VelocityZ']])
    else:
      _nullifyBCDataSetVectors(t, bndType, loc=loc, vectors=[['VelocityX','VelocityY','VelocityZ']]+vectors)
  return None

# Apply init to all bcdataset
def _initBCDataSet(t, varNameString, val1=[], val2=[]):
    zones = Internal.getZones(t)
    if zones == []: zones = [t] # must be a BC node
    for z in zones:
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            datas = Internal.getBCDataSet(z, b)
            fields = []; connects = []
            for d in datas: # build array
              np = d[1].size; ne = np-1
              dim = ['Unstructured', np, ne, 'NODE', 1]
              f = Internal.convertDataNode2Array(d, dim, connects)
              fields.append(f[1])
            if fields != []:
              fields = Converter.addVars(fields)
            else:
              np1 = Internal.getNodeFromName1(b, 'PointList')
              np2 = Internal.getNodeFromName1(b, 'PointRange')
              if np2 is not None:
                win = Internal.range2Window(np2[1])
                imin = win[0]; imax = win[1]
                jmin = win[2]; jmax = win[3]
                kmin = win[4]; kmax = win[5]
                np = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
              elif np1 is not None:
                np = np1[1].size
              else: raise ValueError, 'initBCDataSet: no PointRange or PointList in BC.'
              fields = Converter.array('empty',np,1,1)
            fn = Converter.initVars(fields, varNameString, val1, val2)
            nofld = fn[1].shape[0]-1
            varName = varNameString.split('=')[0]
            varName = varName.replace('{', '')
            varName = varName.replace('}', '')
            varName = varName.replace('centers:', '')
            varName = varName.replace('nodes:', '')
            f = Converter.extractVars(fn, [varName])
            fieldFaceNode = Internal.createDataNode(varName, f, 0, cellDim=1)
            if datas != []:
              cont, c = Internal.getParentOfNode(z, datas[0])
              Internal._createUniqueChild(cont, varName, 'DataArray_t', value=fieldFaceNode[1])
            else:
              d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                        gridLocation='FaceCenter', parent=b)
              d = Internal.newBCData('BCData', parent=d)
              Internal._createUniqueChild(d, varName, 'DataArray_t', value=fieldFaceNode[1])
    return None

# -- randomizeVar: ajoute du bruit sur une variable
def randomizeVar(t, var, deltaMin, deltaMax):
  """Randomize a field defined by var within a range [a-deltamin,a+deltamax]
  Usage: a = randomizeVar(a, var, deltaMin, deltaMax)"""
  tp = Internal.copyRef(t)
  _randomizeVar(tp, var, deltaMin, deltaMax)
  return tp

def _randomizeVar(t, var, deltaMin, deltaMax):
  loc = 'nodes'
  varname = var
  spl = var.split(':')
  if len(spl) != 1:
    if spl[0] == 'centers': loc = 'centers'
    varname = spl[1]

  _TZA(t, loc, loc, Converter.randomizeVar, Converter.randomizeVar, varname, deltaMin, deltaMax)
  return None

# -- fillMissingVariables: remplit les variables manquantes pour que toutes
# les zones aient les memes variables
def fillMissingVariables(t):
  """Fill FlowSolution nodes with variables, such that all the zones have
  the same variables."""
  tp = Internal.copyRef(t)
  _fillMissingVariables(tp)
  return tp

def _fillMissingVariables(t):
  # scan des variables
  varsn = []; varsc = []
  nodes = Internal.getZones(t)
  for z in nodes:
    vars = getVarNames(z, excludeXYZ=True, loc='nodes')[0]
    for var in vars:
      found = 0
      for v in varsn:
        if v == var: found = 1
      if found == 0: varsn.append(var)

  for z in nodes:
    vars = getVarNames(z, excludeXYZ=True, loc='centers')[0]
    for var in vars:
      found = 0
      for v in varsc:
        if v == var: found = 1
      if found == 0: varsc.append(var)

  # add vars
  _addVars(t, varsn+varsc)

  # Reorder vars for all zones
  nodes = Internal.getZones(t)

  varx = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
  for z in nodes:
    # Tri les coordonnees
    cont = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if cont is not None:
      children = cont[2]
      childrenNames = []
      for c in children: childrenNames.append(c[0])
      new = []
      for name in varx:
        try:
          s1 = childrenNames.index(name)
          new.append(children[s1])
        except: pass
      cont[2] = new

  for z in nodes:
    # Tri les variables en noeuds
    cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if cont is not None:
      children = cont[2]
      childrenNames = []
      for c in children: childrenNames.append(c[0])
      new = []
      for name in varsn:
        try:
          s1 = childrenNames.index(name)
          new.append(children[s1])
        except: pass
      pos = 0
      for name in childrenNames:
        try:
          s1 = varsn.index(name)
        except: new.append(children[pos])
        pos += 1
      cont[2] = new

  for z in nodes:
    # tri les variables en centres
    cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if cont is not None:
      children = cont[2]
      childrenNames = []
      for c in children: childrenNames.append('centers:'+c[0])
      new = []
      for name in varsc:
        try:
          s1 = childrenNames.index(name)
          new.append(children[s1])
        except: pass
      pos = 0
      for name in childrenNames:
        try:
          s1 = varsc.index(name)
        except: new.append(children[pos])
        pos += 1
      cont[2] = new

    for z in nodes:
      # reordonne les containers
      c = 0; i0 = -1; i1 = -1; i2 = -1
      for j in z[2]:
        if j[0] == Internal.__GridCoordinates__: i0 = c
        if j[0] == Internal.__FlowSolutionNodes__: i1 = c
        if j[0] == Internal.__FlowSolutionCenters__: i2 = c
        c += 1
      if i0 > i1 and i1 != -1: tmp = z[2][i1]; z[2][i1] = z[2][i0]; z[2][i0] = tmp; tmp = i1; i1 = i0; i0 = tmp
      if i0 > i2 and i2 != -1: tmp = z[2][i2]; z[2][i2] = z[2][i0]; z[2][i0] = tmp; tmp = i2; i2 = i0; i0 = tmp
      if i1 > i2 and i2 != -1: tmp = z[2][i2]; z[2][i2] = z[2][i1]; z[2][i1] = tmp; tmp = i2; i2 = i1; i1 = tmp
  return None

# -- cpVars
def cpVars(t1, var1, t2, var2):
  """Copy field variables."""
  tp2 = Internal.copyRef(t2)
  zones1 = Internal.getZones(t1)
  zones2 = Internal.getZones(tp2)
  if len(zones1) != len(zones2):
    raise TypeError("cpVars: zones in t1 and t2 must be coherents.")
  l = 0
  for z1 in zones1:
    z2 = zones2[l]
    _cpVars__(z1, var1, z2, var2)
    l += 1
  return tp2

# -- cpVars: in place in t2
def _cpVars(t1, var1, t2, var2):
  zones1 = Internal.getZones(t1)
  zones2 = Internal.getZones(t2)
  l = 0
  for z1 in zones1:
    z2 = zones2[l]
    _cpVars__(z1, var1, z2, var2)
    l += 1
  return None

# internal cpVars for zones - check supprimes car utilisation avec le solver
def _cpVars__(z1, var1, z2, var2):
  """Copy variables var1 from z1 to variables var2 of z2.
  Usage: cpVars(z1, var1, z2, var2)"""
  nodes1 = getStdNodesFromName(z1, var1)
  nodes2 = getStdNodesFromName(z2, var2)

  # Switch nodes
  if (nodes2 != []): # var2 existe deja dans z2
    nodes2[0][1] = nodes1[0][1]
    return None
  # Create nodes
  s2 = var2.split(':')
  if (len(s2) == 1 or s2[0] == 'nodes'): loc2 = 'nodes'
  else: loc2 = 'centers'
  (r, c) = Internal.getParentOfNode(z1, nodes1[0])
  if (r[0] == Internal.__GridCoordinates__):
    place2 = Internal.getNodeFromName1(z2, Internal.__GridCoordinates__)
  elif (loc2 == 'nodes'):
    place2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionNodes__)
  elif (loc2 == 'centers'):
    place2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)

  if place2 is None: # Le containeur n'existe pas; create it
    if r[0] == Internal.__GridCoordinates__:
      z2[2].append([Internal.__GridCoordinates__, None, [], 'GridCoordinates_t'])
    elif (loc2 == 'nodes'):
      z2[2].append([Internal.__FlowSolutionNodes__, None, [], 'FlowSolution_t'])
    elif (loc2 == 'centers'):
      z2[2].append([Internal.__FlowSolutionCenters__, None, [], 'FlowSolution_t'])
      f = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)
      Internal._createUniqueChild(f, 'GridLocation', 'GridLocation_t', 'CellCenter')
    h = z2[2][len(z2[2])-1]
  else:
    h = place2
  var2s = var2.split(':')
  if (len(var2s) > 1): var2 = var2s[1]
  h[2].append([var2, nodes1[0][1], nodes1[0][2], nodes1[0][3]])
  return None

# -- extractVars
def extractVars(t, vars):
  """Only keep the given variables in t.
  Usage: extractVars(z, var)"""
  tc = Internal.copyRef(t)
  _extractVars(tc, vars)
  return tc

def _extractVars(t, vars):
  if isinstance(vars, str): vars = [vars]
  zones = Internal.getZones(t)
  for z in zones:
    varNames = getVarNames(z)[0]
    for v in varNames:
      if v not in vars: _rmVars(z, v)
  return None

# -- rmVars
def rmVars(z, var):
  """Remove variables var from t.
  Usage: rmVars(t, var)"""
  zc = Internal.copyRef(z)
  _rmVars(zc, var)
  return zc

def _rmVars(z, var):
  zn = Internal.getZones(z)
  if zn == []: return None
  for i in zn:
    if isinstance(var, list):
      for v in var:
        nodes = getStdNodesFromName(i, v)
        if nodes != []:
          (parent, d) = Internal.getParentOfNode(i, nodes[0])
          del parent[2][d]
    else:
      nodes = getStdNodesFromName(i, var)
      if nodes != []:
        (parent, d) = Internal.getParentOfNode(i, nodes[0])
        del parent[2][d]
  return None

# -- normalize: normalise un jeu de variables
def normalize(t, vars):
  """Normalize the field defined by vars in the tree.
  Usage: normalize(t, vars)"""
  tp = Internal.copyRef(t)
  _normalize(tp, vars)
  return tp

def _normalize(t, vars):
  loc = ''
  vars2 = []
  for v in vars:
    s = v.split(':')
    if len(s) == 2 and s[0] == 'centers': #centres
      if loc == '': loc = s[0]
      elif loc != s[0]: raise ValueError("normalize: vector components must have the same location (centers or nodes).")
      vars2.append(s[1])
    elif (len(s) == 1 or (len(s) == 2 and s[0] == 'nodes')): #noeuds
      if loc == '': loc = 'nodes'
      elif loc == 'centers': raise ValueError("normalize: vector components must have the same location (centers or nodes).")
      if len(s) == 2: vars2.append(s[1])
      else: vars2.append(s[0])
    else:
      raise ValueError("normalize: invalid vector component.")
  _TZA(t, loc, loc, Converter.normalize, Converter.normalize, vars2)
  return None

# -- magnitude: calcul la norme d'un jeu de variables
def magnitude(t, vars):
  """Compute the magnitude of the fields defined by vars in the tree.
  Usage: magnitude(t, vars)"""
  tp = Internal.copyRef(t)
  _magnitude(tp, vars)
  return tp

def _magnitude(t, vars):
  loc = ''
  vars2 = []
  for v in vars:
    s = v.split(':')
    if len(s) == 2: #centres
      if loc == '': loc = s[0]
      elif loc != s[0]: raise ValueError("magnitude: vector components must have the same location (centers or nodes).")
      vars2.append(s[1])
    elif len(s) == 1: #noeuds
      if loc == '': loc = 'nodes'
      elif loc == 'centers': raise ValueError("magnitude: vector components must have the same location (centers or nodes).")
    else:
      raise ValueError("magnitude: invalid vector component.")

  if loc == 'nodes':
    _TZA(t, loc, loc, Converter.magnitude, Converter.magnitude, vars)
  else:
    _TZA(t, loc, loc, Converter.magnitude, Converter.magnitude, vars2)
  return None

# -- normL0
def normL0(t, var):
  """Get the L0 norm of the field defined by varName in t.
  If celln exists in the array, the norm for blanked points is not computed.
  Usage: normL0(t, varName)"""
  A = getField(var, t)
  v = var.split(':')
  if len(v) > 1: var = v[1]
  return Converter.normL0(A, var)

# -- normL2
def normL2(t, var):
  """Get the L2 norm of the field defined by varName in t.
  If celln exists in the array, the norm for blanked points is not computed.
  Usage: normL0(t, varName)"""
  A = getField(var, t)
  v = var.split(':')
  if len(v) > 1: var = v[1]
  return Converter.normL2(A, var)

# -- getArgMin
def getArgMin(t, var):
  """Get value where the variable defined by varName is minimum.
  Usage: getArgMin(t, var)"""
  A = getField(var, t)
  v = var.split(':')
  if len(v) > 1: var = v[1]
  return Converter.getArgMin(A, var)

# -- getArgMax
def getArgMax(t, var):
  """Get value where the variable defined by varName is maximum.
  Usage: getArgMin(t, var)"""
  A = getField(var, t)
  v = var.split(':')
  if len(v) > 1: var = v[1]
  return Converter.getArgMax(A, var)

# -- getMinValue
def getMinValue(t, varName):
  """Get the minimum value of variable defined by var.
  Usage: getMinValue(t, var)"""
  if not isinstance(varName, list): varNames = [varName]
  else: varNames = varName
  out = []
  if varNames[0] == Internal.__GridCoordinates__: varNames = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
  for v in varNames:
    A = getField(v, t)
    va = v.split(':')
    if len(va) > 1: v = va[1]
    out.append(Converter.getMinValue(A, v))
  if len(out) == 1: return out[0]
  else: return out

# -- getMaxValue
def getMaxValue(t, varName):
  """Get the maximum value of variable defined by var.
  Usage: getMaxValue(t, var)"""
  if not isinstance(varName, list): varNames = [varName]
  else: varNames = varName
  out = []
  if varNames[0] == Internal.__GridCoordinates__: varNames = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
  for v in varNames:
    A = getField(v, t)
    va = v.split(':')
    if len(va) > 1: v = va[1]
    out.append(Converter.getMaxValue(A, v))
  if len(out) == 1: return out[0]
  else: return out

# -- getMeanValue
def getMeanValue(t, var):
  """Get the mean value of variable defined by var.
  Usage: getMeanValue(t, var)"""
  A = getField(var, t)
  v = var.split(':')
  if len(v) > 1: var = v[1]
  return Converter.getMeanValue(A, var)

# -- getMeanRangeValue
def getMeanRangeValue(t, var, rmin, rmax):
  """Get the mean value of variable defined by var in the sorted range between rmin and rmax.
  Usage: getMeanRangeValue(t, var, rmin, rmax)"""
  A = getField(var, t)
  v = var.split(':')
  if len(v) > 1: var = v[1]
  return Converter.getMeanRangeValue(A, var, rmin, rmax)

#==============================================================================
# -- Topology conversions --
#==============================================================================

# -- Convert a BAR array without branches, closed into an i-array
def convertBAR2Struct(t):
  """Convert a BAR array without branches, closed into an i-array.
  Usage: convertBAR2Struct(t)"""
  return TZA(t, 'both', 'both', Converter.convertBAR2Struct, None)

def _convertBAR2Struct(t):
  _TZA(t, 'both', 'both', Converter.convertBAR2Struct, None)
  return None

# -- convertArray2Tetra
# split='simple': pas d'ajout de points
# split='withBarycenters': ajout de points aux centres des elements et des faces
def convertArray2Tetra(t, split='simple'):
  """Convert a zone to an unstructured zone.
  Unstructured array can be triangular in 2D and tetraedrical in 3D.
  Usage: convertArray2Tetra(t, split)"""
  tp = Internal.copyRef(t)
  _convertArray2Tetra(tp, split)
  return tp

def _convertArray2Tetra(t, split='simple'):
  _deleteZoneBC__(t)
  _deleteGridConnectivity__(t)
  if split == 'simple':
    _TZANC(t, 'both', 'both', Converter.convertArray2Tetra,
           Converter.convertArray2Tetra, split, split)
  else:
    nodes = Internal.getZones(t)
    for z in nodes:
      fieldn = getAllFields(z, 'nodes')[0]
      fieldc = getAllFields(z, 'centers')[0]
      if fieldc == []:
        res = Converter.convertArray2Tetra1__(fieldn,split='withBarycenters')
        z = setFields([res], z, 'nodes', writeDim=True)
      else:
        res = Converter.convertArray2Tetra1__(fieldn, fieldc, split='withBarycenters')        
        z = setFields([res[0]], z, 'nodes', writeDim=True)
        z = setFields([res[1]], z, 'centers', writeDim=False)
  return None

# -- convertArray2Hexa
def convertArray2Hexa(t):
  """Convert a structured zone to an unstructured quad/hexa zone.
  Convert an unstructured zone to a quad/hexa zone. If the original
  zone is a TRI,TETRA or PENTA zone, return a QUAD/HEXA/HEXA zone with
  degenerated edges.
  Usage: convertArray2Hexa(t)"""
  tp = Internal.copyRef(t)
  _convertArray2Hexa(tp)
  return tp

def _convertArray2Hexa(t):
  _deleteZoneBC__(t)
  _deleteGridConnectivity__(t)
  _TZA(t, 'both', 'both', Converter.convertArray2Hexa, None)
  return None

# -- convertArray2NGon
def convertArray2NGon(t, recoverBC=True):
  """Convert a zone to a NGON zone.
  Usage: convertArray2NGon(t)"""
  tp = Internal.copyRef(t)
  _convertArray2NGon(tp,recoverBC=recoverBC)
  return tp

def _convertArray2NGon(t,recoverBC=True):
  zones = Internal.getZones(t)
  if recoverBC:
    gbcs = [] 
    for z in zones:
      dims = Internal.getZoneDim(z)
      if dims[0] == 'Unstructured' and dims[3] == 'NGON': gbcs.append([])
      else: gbcs.append(getBCs(z))
  else: _deleteZoneBC__(t)
  _deleteGridConnectivity__(t)
  _TZA(t, 'both', 'both', Converter.convertArray2NGon, None)

  # Recover BCs for NGon
  if recoverBC:
    zones = Internal.getZones(t); c = 0
    for z in zones:
      if (gbcs[c] != [] and len(gbcs[c][0]) > 0): _recoverBCs(z, gbcs[c])
      c += 1

  return None

# -- convertArray2Node
def convertArray2Node(t):
  """Convert an array to a node array.
  Usage: convertArray2Node(t)"""
  tp = Internal.copyRef(t)
  _convertArray2Node(tp)
  return tp

def _convertArray2Node(t):
  _deleteFlowSolutions__(t, 'centers')
  _deleteZoneBC__(t)
  _deleteGridConnectivity__(t)
  #_TZA(t, 'nodes', 'nodes', Converter.convertArray2Node, None)
  zones = Internal.getZones(t)
  for z in zones:
    n = Internal.getNodeFromName1(z, 'ZoneType')
    Internal.setValue(n, 'Unstructured')
    x = Internal.getNodeFromName2(z, 'CoordinateX')
    z[1] = numpy.empty((1,3), numpy.int32, order='Fortran')
    z[1][0,0] = x[1].size; z[1][0,1] = 0; z[1][0,2] = 0
    Internal._rmNodesByType(z, 'Elements_t')
    n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
    Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
    Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
  return None

# -- convertTri2Quad
def convertTri2Quad(z, alpha=30.):
  """Convert a TRI zone to a QUAD zone.
  Usage: convertTri2Quad(z, angle)"""
  a = getAllFields(z, 'nodes')
  b, c = Converter.convertTri2Quad(a, alpha)
  z1 = convertArrays2ZoneNode(getZoneName('quad'), [b])
  z2 = convertArrays2ZoneNode(getZoneName('tri'), [c])
  return z1, z2

# -- conformizeNGon
def conformizeNGon(a, tol=1.e-6):
    """Conformize topologically a NGON zone.
    Usage: conformizeNGon(a, tol)"""
    return TZGC(a, 'nodes', Converter.conformizeNGon, tol)

def _conformizeNGon(a, tol=1.e-6):
    _TZGC(a, 'nodes', Converter.conformizeNGon, tol)
    return None

#=============================================================================
# -- Create BC(s) to a zone node --
#=============================================================================

# -- addBC2Zone
def addBC2Zone(zone, bndName, bndType, range=[],
               zoneDonor=[], rangeDonor=[], trirac=[1,2,3],
               rotationCenter=[], rotationAngle=[], translation=[],
               faceList=[], elementList=[], elementRange=[], data=None,
               subzone=None, faceListDonor=None, elementListDonor=None,
               elementRangeDonor=None, tol=1.e-12):
  """Add a BC to a zone node.
  Usage: addBC2Zone(zone, bndName, bndType, range)"""
  zonep = Internal.copyRef(zone)
  _addBC2Zone(zonep, bndName, bndType, range,
              zoneDonor, rangeDonor, trirac,
              rotationCenter, rotationAngle, translation,
              faceList, elementList, elementRange, data, subzone,
              faceListDonor, elementListDonor, elementRangeDonor, tol)
  return zonep

def _addBC2Zone(a, bndName, bndType, range=[],
                zoneDonor=[], rangeDonor=[], trirac=[1,2,3],
                rotationCenter=[], rotationAngle=[], translation=[],
                faceList=[], elementList=[], elementRange=[], data=None,
                subzone=None, faceListDonor=None, elementListDonor=None,
                elementRangeDonor=None, tol=1.e-12):
  bndName = getBCName(bndName)
  zones = Internal.getZones(a)
  for z in zones:
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured':
      eltType = dims[3]
      if eltType == 'NGON':
        if (faceList == [] and subzone is None):
          raise TypeError("addBC2Zone: unstructured grids requires a faceList or a subzone.")
        _addBC2NGonZone__(z, bndName, bndType, faceList, data, subzone,
                          zoneDonor, faceListDonor,
                          rotationCenter, rotationAngle, translation, tol)
      else: # basic elements
        if (elementList == [] and elementRange == [] and subzone is None and faceList == []):
          raise TypeError("addBC2Zone: unstructured grids requires a elementList, a elementRange or a subzone.")
        _addBC2UnstructZone__(z, bndName, bndType, elementList, elementRange,
                              faceList, data, subzone,
                              zoneDonor, elementListDonor, elementRangeDonor,
                              faceListDonor, rotationCenter, rotationAngle,
                              translation)
    else: # structured
      _addBC2StructZone__(z, bndName, bndType, range, faceList,
                          zoneDonor, rangeDonor, faceListDonor, trirac,
                          rotationCenter, rotationAngle, translation, data)
  return None

# -- Ajout Infos peridiques pour les BC match periodiques
def _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation):
  if (rotationCenter != [] or rotationAngle != [] or translation != []):
    if rotationCenter == []: rotCenter = (0,0,0)
    else: rotCenter = rotationCenter
    if rotationAngle == []: rotAngle = (0,0,0)
    else: rotAngle = rotationAngle
    if translation == []: trans = (0,0,0)
    else: trans = translation

    info[2].append(['GridConnectivityProperty', None, [], 'GridConnectivityProperty_t'])
    info = info[2][len(info[2])-1]
    info[2].append(['Periodic', None, [], 'Periodic_t'])
    info = info[2][0]
    if len(rotCenter) != 3:
      raise ValueError("addBC2Zone: rotationCenter must be (Cx, Cy, Cz).")
    if len(rotAngle) != 3:
      raise ValueError("addBC2Zone: rotationAngle must be (alpha, beta, gamma).")
    if len(trans) != 3:
      raise ValueError("addBC2Zone: translation must be (tx, ty, tz).")
    v = numpy.zeros((3), numpy.float64)
    v[0] = rotCenter[0]; v[1] = rotCenter[1]; v[2] = rotCenter[2]
    info[2].append(['RotationCenter', v, [], 'DataArray_t'])
    v = numpy.zeros((3), numpy.float64)
    v[0] = rotAngle[0]; v[1] = rotAngle[1]; v[2] = rotAngle[2]
    info[2].append(['RotationAngle', v, [], 'DataArray_t'])
    v = numpy.zeros((3), numpy.float64)
    v[0] = trans[0]; v[1] = trans[1]; v[2] = trans[2]
    info[2].append(['Translation', v, [], 'DataArray_t'])
  return None


# typeZone=0 : struct, 1 : NGON, 2: unstructured BE
def _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=0, faceList=[], elementList=[], elementRange=[], pointRange=[], zoneDonor=[]):
  zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
  if zoneGC == []:
    z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
    zoneGC = z[2][len(z[2])-1]
  else:
    zoneGC = zoneGC[0]
  # Cree le noeud de GC
  if zoneDonor == []:
    # autoattach
    v = numpy.fromstring(z[0], 'c')
    zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])
  else:
    # donors donnes
    st = ""
    for i in zoneDonor:
      if isinstance(i, str):
        if st == "": st = i
        else: st = st+","+i
      else:
        if st == "": st = i[0]
        else: st = st+","+i[0]
    v = numpy.fromstring(st, 'c')
    zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])
  info = zoneGC[2][len(zoneGC[2])-1]

  if typeZone == 0: # STRUCTURED
    if pointRange==[]: raise ValueError("_addFamilyOfStageGC__: pointRange is empty.")
    r = Internal.window2Range(pointRange)
    info[2].append(['PointRange', r, [], 'IndexRange_t'])
  
  elif typeZone == 1: # NGON
    if faceList==[]: raise ValueError("_addFamilyOfStageGC__: faceList is empty.")

    if isinstance(faceList, numpy.ndarray): r = faceList
    else: r = numpy.array(faceList, dtype=numpy.int32)
    r = r.reshape((1,r.size), order='Fortran')
    info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

  elif typeZone == 2: #UNS BE
    if (elementList != []):
      if isinstance(elementList, numpy.ndarray): r = elementList
      else: r = numpy.array(elementList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__ELEMENTLIST__, r, [], 'IndexArray_t'])
    elif (elementRange != []):
      r = numpy.empty((1,2), numpy.int32, order='Fortran')
      r[0,0] = elementRange[0]
      r[0,1] = elementRange[1]
      info[2].append([Internal.__ELEMENTRANGE__, r, [], 'IndexRange_t'])
    elif (faceList != []):
      v = numpy.fromstring('FaceCenter', 'c')
      info[2].append(['GridLocation', v, [], 'GridLocation_t'])
      if isinstance(faceList, numpy.ndarray): r = faceList
      else: r = numpy.array(faceList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
    else:
      raise ValueError("_addFamilyOfStageGC__: elementList, elementRange and faceList are all empty.")

  else:
    raise ValueError("_addFamilyOfStageGC__: typeZone not valid.")

  v = numpy.fromstring('Abutting', 'c')
  info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])
  v = numpy.fromstring(bndType2, 'c')
  info[2].append(['FamilyStage', v, [], 'FamilyName_t'])
  return None

# -- addBC2Zone pour les grilles structurees
def _addBC2StructZone__(z, bndName, bndType, range=[], faceList=[],
                        zoneDonor=[], rangeDonor=[], faceListDonor=[], trirac=[1,2,3],
                        rotationCenter=[], rotationAngle=[],
                        translation=[], data=None):
  # Type de condition aux limites definies par une famille
  # par ex: FamilySpecified:OVERLAP
  s = bndType.split(':')
  bndType1 = s[0]
  if len(s) > 1: bndType2 = s[1]
  else: bndType2 = ''

  typeR = -1 # 0  : PR, 1 PL
  if range != []: typeR=0 
  elif faceList != []: typeR=1
  else:
    raise ValueError("addBC2Zone: match connectivity requires a range or face list.")

  # Range defini par une chaine ou non
  if typeR==0:
    if isinstance(range, str):      
      range = convertStringRange2Range__(range, z)#fenetre complete
    if isinstance(rangeDonor, str):
      if isinstance(zoneDonor, str):
        raise ValueError("addBC2Zone: donor range must be explicitly specified.")
      else:
        rangeDonor = convertStringRange2Range__(rangeDonor, zoneDonor)
  
  if bndType1 == 'BCMatch':
    if typeR==0 and rangeDonor == [] and trirac==[]: 
      raise ValueError("addBC2Zone: match connectivity requires a donor point range and a trirac.")
    elif typeR==1 and faceListDonor == []: 
      raise ValueError("addBC2Zone: match connectivity requires a donor face list.")

    if (zoneDonor == []):
      raise ValueError("addBC2Zone: match connectivity requires a donor zone.")
    # Cree le noeud zoneGridConnectivity si besoin
    zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    if zoneGC == []:
      z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
      zoneGC = z[2][len(z[2])-1]
    else:
      zoneGC = zoneGC[0]
    # Cree le noeud de GC1-1
    if isinstance(zoneDonor, str): v = numpy.fromstring(zoneDonor, 'c')
    else: v = numpy.fromstring(zoneDonor[0], 'c')
    zoneGC[2].append([bndName, v, [], 'GridConnectivity1to1_t']); l = len(zoneGC[2])
    info = zoneGC[2][l-1]

    if typeR==0:
      r = Internal.window2Range(range)
      info[2].append([Internal.__FACERANGE__, r, [], 'IndexRange_t'])
      o = Internal.window2Range(rangeDonor)
      info[2].append([Internal.__FACERANGE__+'Donor', o, [], 'IndexRange_t'])
      size = len(trirac)
      o = numpy.zeros((size), numpy.int32)
      for i in xrange(size): o[i] = trirac[i]
      info[2].append(['Transform', o, [], '\"int[IndexDimension]\"'])

    elif typeR==1:
      v = numpy.fromstring('FaceCenter', 'c')
      info[2].append(['GridLocation', v, [], 'GridLocation_t'])
      if isinstance(faceList, numpy.ndarray): r = faceList
      else: r = numpy.array(faceList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
      if isinstance(faceListDonor, numpy.ndarray): r = faceListDonor
      else: r = numpy.array(faceListDonor, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])
      v = numpy.fromstring('Abutting1to1', 'c')
      info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])

    # Ajout pour les BC match periodiques
    _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation)

  elif bndType1 == 'BCNearMatch':
    if typeR != 0: 
      raise ValueError("addBC2Zone: nearmatch connectivity requires a point range.")

    if (rangeDonor == [] or zoneDonor == [] or trirac == []):
      raise ValueError("addBC2Zone: nearmatch connectivity requires a donor and transform.")

    # Cree le noeud zoneGridConnectivity si besoin
    zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    if zoneGC == []:
      z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
      zoneGC = z[2][len(z[2])-1]
    else:
      zoneGC = zoneGC[0]
    if isinstance(zoneDonor, str):
      v = numpy.fromstring(zoneDonor, 'c')
    else: v = numpy.fromstring(zoneDonor[0], 'c')
    zoneGC[2].append([bndName, v, [], 'GridConnectivity_t']); l = len(zoneGC[2])

    info = zoneGC[2][l-1]
    r = Internal.window2Range(range)
    info[2].append(['PointRange', r, [], 'IndexRange_t'])
    v = numpy.fromstring('Abutting', 'c')
    info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])
    c = numpy.ones((3,1), numpy.int32)
    info[2].append(['PointListDonor', c, [], 'IndexArray_t'])

    # Nearmatch: rajoute les noeuds PointRangeDonor et Transform en UserDefinedData
    info[2].append(['UserDefinedData', None, [], 'UserDefinedData_t'])
    l = len(info[2])-1; info = info[2][l]
    dd = Internal.window2Range(rangeDonor)
    info[2].append(['PointRangeDonor', dd, [], 'DataArray_t'])
    size = len(trirac); o = numpy.zeros((size), numpy.int32)
    for i in xrange(size): o[i] = trirac[i]
    info[2].append(['Transform', o, [], 'DataArray_t'])
    ratio = getNMRatio__(range, rangeDonor, trirac)
    size = len(ratio); ro = numpy.zeros((size), numpy.float64)
    for i in xrange(size): ro[i] = ratio[i]
    info[2].append(['NMRatio', ro, [], 'DataArray_t'])

  elif bndType1 == 'BCOverlap':
    # Cree le noeud zoneGridConnectivity si besoin
    zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    if zoneGC == []:
      z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
      zoneGC = z[2][len(z[2])-1]
    else:
      zoneGC = zoneGC[0]
    # Cree le noeud de GC
    if zoneDonor == []:
      # autoattach
      v = numpy.fromstring(z[0], 'c')
      zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])
    else:
      # donors donnes
      st = ""
      for i in zoneDonor:
        if isinstance(i, str):
          if st == "": st = i
          else: st = st+","+i
        else:
          if st == "": st = i[0]
          else: st = st+","+i[0]
      v = numpy.fromstring(st, 'c')
      zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])
    info = zoneGC[2][len(zoneGC[2])-1]
    d = numpy.ones((3,1), numpy.int32)
    c = numpy.ones((3,1), numpy.float64)
    r = Internal.window2Range(range)
    info[2].append(['PointRange', r, [], 'IndexRange_t'])
    v = numpy.fromstring('Overset', 'c')
    info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])
    info[2].append(['CellListDonor', d, [], 'IndexArray_t'])
    info[2].append(['InterpolantsDonor', c, [], 'DataArray_t'])
    if rangeDonor == 'doubly_defined':
      info[2].append(['UserDefinedData', None, [], 'UserDefinedData_t'])
      l = len(info[2])-1; info = info[2][l]
      dd = numpy.ones((1), numpy.int32)
      info[2].append(['doubly_defined', dd, [], 'DataArray_t'])

  elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
    _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=0, pointRange=range, zoneDonor=zoneDonor)

  else: # classical BC
    # Cree le noeud zoneBC si besoin
    zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
    if zoneBC == []:
      z[2].append(['ZoneBC', None, [], 'ZoneBC_t'])
      zoneBC = z[2][len(z[2])-1]
    else: zoneBC = zoneBC[0]
    # Cree le noeud de BC
    v = numpy.fromstring(bndType1, 'c')
    zoneBC[2].append([bndName, v, [], 'BC_t']); l = len(zoneBC[2])
    info = zoneBC[2][l-1]
    r = Internal.window2Range(range)
    info[2].append(['PointRange', r, [], 'IndexRange_t'])

    # Ajoute la famille si necessaire
    if bndType1 == 'FamilySpecified':
      v = numpy.fromstring(bndType2, 'c')
      info[2].append(['FamilyName', v, [], 'FamilyName_t'])

    # Ajoute les Data si necessaire (donnees Dirichlet)
    if data is not None:
      node1 = Internal.createNode('State', 'DataArray_t', value=data)
      node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
      node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
      info[2].append(node3)
  return None

# -- addBC2Zone pour les zones NGon
def _addBC2NGonZone__(z, bndName, bndType, faceList, data, subzone,
                      zoneDonor, faceListDonor,
                      rotationCenter, rotationAngle, translation, tol):
  # Type de condition aux limites definies par une famille
  s = bndType.split(':')
  bndType1 = s[0]
  if (len(s) > 1): bndType2 = s[1]
  else: bndType2 = ''

  # si subzone: on cree le faceList par identification
  if subzone is not None:
    hook = createHook(z, 'faceCenters')
    faceList = identifyElements(hook, subzone, tol)
    freeHook(hook)

  if (bndType1 == 'BCMatch'):
    if (zoneDonor == [] or (faceListDonor is None and subzone is None)):
      raise ValueError("addBC2Zone: NGON match connectivity requires a donor and a faceListDonor or a subzone.")
    # cree le faceListDonor si subzone fourni
    if (subzone is not None and faceListDonor is None):
      hook = createHook(zoneDonor, 'faceCenters')
      faceListDonor = identifyElements(hook, subzone, tol)
      freeHook(hook)

    # Cree le noeud zoneGridConnectivity si besoin
    zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    if (zoneGC == []):
      z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
      zoneGC = z[2][len(z[2])-1]
    else:
      zoneGC = zoneGC[0]

    if isinstance(zoneDonor, str): v = numpy.fromstring(zoneDonor, 'c')
    else: v = numpy.fromstring(zoneDonor[0], 'c')
    zoneGC[2].append([bndName, v, [], 'GridConnectivity_t']); l = len(zoneGC[2])
    info = zoneGC[2][l-1]
    v = numpy.fromstring('FaceCenter', 'c')
    info[2].append(['GridLocation', v, [], 'GridLocation_t'])
    if isinstance(faceList, numpy.ndarray): r = faceList
    else: r = numpy.array(faceList, dtype=numpy.int32)
    r = r.reshape((1,r.size), order='Fortran')
    info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
    if isinstance(faceListDonor, numpy.ndarray): r = faceListDonor
    else: r = numpy.array(faceListDonor, dtype=numpy.int32)
    r = r.reshape((1,r.size), order='Fortran')
    info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])
    v = numpy.fromstring('Abutting1to1', 'c')
    info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])
    # Ajout pour les BC match periodiques
    _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation)

  elif bndType1 == 'BCNearMatch':
    print 'addBC2Zone: BCNearMatch not valid for NGON zones.'

  # elif bndType1 == 'BCOverlap':
  #   # Cree le noeud zoneGridConnectivity si besoin
  #   zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
  #   if (zoneGC == []):
  #     z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
  #     zoneGC = z[2][len(z[2])-1]
  #   else:
  #     zoneGC = zoneGC[0]
  #   # Cree le noeud de GC
  #   if (zoneDonor == []):
  #     # autoattach
  #     v = numpy.fromstring(z[0], 'c')
  #     zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])
  #   else:
  #     # donors donnes
  #     st = ""
  #     for i in zoneDonor:
  #       if isinstance(i, str):
  #         if (st == ""): st = i
  #         else: st = st+","+i
  #       else:
  #         if (st == ""): st = i[0]
  #         else: st = st+","+i[0]
  #     v = numpy.fromstring(st, 'c')
  #     zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])

    # l = len(zoneGC[2])
    # info = zoneGC[2][l-1]
    # v = numpy.fromstring('FaceCenter', 'c')
    # info[2].append(['GridLocation', v, [], 'GridLocation_t'])
    # if isinstance(faceList, numpy.ndarray): r = faceList
    # else: r = numpy.array(faceList, dtype=numpy.int32)
    # info[2].append(['PointList', r, [], 'IndexArray_t'])
    # v = numpy.fromstring('Overset', 'c')
    # info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])

    # d = numpy.ones((3,1), numpy.int32)
    # c = numpy.ones((3,1), numpy.float64)
    # info[2].append(['PointListDonor', d, [], 'IndexArray_t'])
    # info[2].append(['InterpolantsDonor', c, [], 'DataArray_t'])


  elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
    _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=1, faceList=faceList, zoneDonor=zoneDonor)

  else: # BC classiques
    # Cree le noeud zoneBC si besoin
    zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
    if zoneBC == []:
      z[2].append(['ZoneBC', None, [], 'ZoneBC_t'])
      zoneBC = z[2][len(z[2])-1]
    else: zoneBC = zoneBC[0]

    # Cree le noeud de BC
    v = numpy.fromstring(bndType1, 'c')
    zoneBC[2].append([bndName, v, [], 'BC_t']); l = len(zoneBC[2])
    info = zoneBC[2][l-1]
    v = numpy.fromstring('FaceCenter', 'c')
    info[2].append(['GridLocation', v, [], 'GridLocation_t'])
    if isinstance(faceList, numpy.ndarray): r = faceList
    else: r = numpy.array(faceList, dtype=numpy.int32)
    r = r.reshape((1,r.size), order='Fortran')
    info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

    # Ajoute la famille si necessaire
    if bndType1 == 'FamilySpecified':
      v = numpy.fromstring(bndType2, 'c')
      info[2].append(['FamilyName', v, [], 'FamilyName_t'])

    # Ajoute les Data si necessaire (donnees Dirichlet)
    if data is not None:
      node1 = Internal.createNode('State', 'DataArray_t', value=data)
      node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
      node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
      info[2].append(node3)
  return None

# -- addBC2Zone pour les zones non structurees a elements basiques
def _addBC2UnstructZone__(z, bndName, bndType, elementList, elementRange,
                          faceList, data, subzone,
                          zoneDonor, elementListDonor, elementRangeDonor, faceListDonor,
                          rotationCenter, rotationAngle, translation):
  s = bndType.split(':')
  bndType1 = s[0]
  if len(s) > 1: bndType2 = s[1]
  else: bndType2 = ''

  # si subzone: on cree la connectivite BC, on passe en elementRange
  if subzone is not None:
    bcn = Internal.getNodeFromName1(z, subzone[0])
    if bcn is None:
      _mergeConnectivity(z, subzone, boundary=1)
    bcn = Internal.getNodeFromName1(z, subzone[0])
    bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
    elementRange = [bcnr[1][0], bcnr[1][1]]

  if (bndType1 == 'BCMatch'):
    if (zoneDonor == [] or
        faceListDonor is None and subzone is None and elementListDonor is None and elementRangeDonor is None):
      raise ValueError("addBC2Zone: unstructured match connectivity requires a donor zone and a faceListDonor or a subzone or an elementRangeDonor or an elementListDonor.")
    # si subzone fournie: cree le elementRangeDonor
    if subzone is not None:
      bcn = Internal.getNodeFromName1(zoneDonor, subzone[0])
      if bcn is None:
        _mergeConnectivity(zoneDonor, subzone, boundary=1)
      bcn = Internal.getNodeFromName1(zoneDonor, subzone[0])
      bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
      elementRangeDonor = [bcnr[1][0], bcnr[1][1]]

    # Cree le noeud zoneGridConnectivity si besoin
    zoneGC = Internal.createUniqueChild(z, 'ZoneGridConnectivity',
                                        'ZoneGridConnectivity_t')

    if isinstance(zoneDonor, str): v = numpy.fromstring(zoneDonor, 'c')
    else: v = numpy.fromstring(zoneDonor[0], 'c')
    zoneGC[2].append([bndName, v, [], 'GridConnectivity_t']); l = len(zoneGC[2])
    info = zoneGC[2][l-1]

    if (elementList != []):
      if isinstance(elementList, numpy.ndarray): r = elementList
      else: r = numpy.array(elementList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__ELEMENTLIST__, r, [], 'IndexArray_t'])
    elif (elementRange != []):
      r = numpy.empty((1,2), numpy.int32, order='Fortran')
      r[0,0] = elementRange[0]
      r[0,1] = elementRange[1]
      info[2].append([Internal.__ELEMENTRANGE__, r, [], 'IndexRange_t'])
    elif (faceList != []):
      v = numpy.fromstring('FaceCenter', 'c')
      info[2].append(['GridLocation', v, [], 'GridLocation_t'])
      if isinstance(faceList, numpy.ndarray): r = faceList
      else: r = numpy.array(faceList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

    if elementListDonor is not None:
      if isinstance(elementListDonor, numpy.ndarray):
        r = elementListDonor
      else: r = numpy.array(elementListDonor, dtype=numpy.int32)
      r = r.reshape((1,r.size))
      info[2].append([Internal.__ELEMENTLIST__+'Donor', r, [], 'IndexArray_t'])
    elif elementRangeDonor is not None:
      r = numpy.empty((1,2), numpy.int32, order='Fortran')
      r[0,0] = elementRangeDonor[0]
      r[0,1] = elementRangeDonor[1]
      info[2].append([Internal.__ELEMENTRANGE__+'Donor', r, [], 'IndexRange_t'])
    elif faceListDonor is not None:
      if isinstance(faceListDonor, numpy.ndarray): r = faceList
      else: r = numpy.array(faceListDonor, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])

    v = numpy.fromstring('Abutting1to1', 'c')
    info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])
    # Ajout pour les BC match periodiques
    _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation)

  elif (bndType1 == 'BCOverlap'):
    print 'addBC2Zone: BCOverlap not valid for unstructured zones.'

  elif (bndType1 == 'BCNearMatch'):
    print 'addBC2Zone: BCNearMatch not valid for unstructured zones.'

  elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
    _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=2, elementRange=elementRange,
                          elementList=elementList, faceList=faceList, zoneDonor=zoneDonor)

  else: # BC classique
    # Cree le noeud zoneBC si besoin
    zoneBC = Internal.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')

    # Cree le noeud de BC
    info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType1)
    if elementList != []:
      if isinstance(elementList, numpy.ndarray): r = elementList
      else: r = numpy.array(elementList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      Internal.createChild(info, INTERNAL.__ELEMENTLIST__, 'IndexArray_t', value=r)
    elif elementRange != []:
      n = numpy.empty((1,2), numpy.int32, order='Fortran')
      n[0,0] = elementRange[0]
      n[0,1] = elementRange[1]
      Internal.createUniqueChild(info, Internal.__ELEMENTRANGE__,
                                 'IndexRange_t', value=n)
    elif (faceList != []):
      Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                                 value='FaceCenter')
      if isinstance(faceList, numpy.ndarray): r = faceList
      else: r = numpy.array(faceList, dtype=numpy.int32)
      r = r.reshape((1,r.size), order='Fortran')
      info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

    # Ajoute la famille si necessaire
    if bndType1 == 'FamilySpecified':
      v = numpy.fromstring(bndType2, 'c')
      info[2].append(['FamilyName', v, [], 'FamilyName_t'])

    # Ajoute les Data si necessaire (donnees Dirichlet)
    if data is not None:
      node1 = Internal.createNode('State', 'DataArray_t', value=data)
      node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
      node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
      info[2].append(node3)
  return None

# -- converti une chaine en range (z Structure seulement)
# IN: 'imin', 'jmin', ...
# OUT: [imin,imax,jmin,jmax,kmin,kmax]
def convertStringRange2Range__(range, z):
  if range == 'doubly_defined': return range
  dim = Internal.getZoneDim(z)
  typeGrid = dim[0]
  if typeGrid == 'Unstructured':
    raise TypeError("addBC2Zone: range invalid for unstructured grids.")
  if range == 'imin': range = [1,1,1,dim[2],1,dim[3]]
  elif range == 'imax': range = [dim[1],dim[1],1,dim[2],1,dim[3]]
  elif range == 'jmin': range = [1,dim[1],1,1,1,dim[3]]
  elif range == 'jmax': range = [1,dim[1],dim[2],dim[2],1,dim[3]]
  elif range == 'kmin': range = [1,dim[1],1,dim[2],1,1]
  else: range = [1,dim[1],1,dim[2],dim[3],dim[3]]
  return range

# -- converti un BCRange en tableau de BCFaces
def convertBCRange2BCFace__(ni, nj, nk, range0):
  win = Internal.range2Window(range0)
  imin = win[0]; imax = win[1]
  jmin = win[2]; jmax = win[3]
  kmin = win[4]; kmax = win[5]

  dir = 0
  di = imax-imin; dj = jmax-jmin; dk = kmax-kmin
  if di == 0:
    di = 1
    if ni != 1: dir = 1
  if dj == 0:
    dj = 1
    if nj != 1: dir = 2
  if dk == 0:
    dk = 1
    if nk != 1: dir = 3

  if dir == 1: nfaces = dj*dk
  elif dir == 2: nfaces = di*dk
  else: nfaces = di*dj
  ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
  ninti = ni*nj1*nk1; nintj = ni1*nj*nk1
  bcfaces = numpy.empty((nfaces),numpy.int32)
  ninj = ni*nj
  nof = 0
  if nk == 1: kmax += 1
  if nj == 1: jmax += 1
  if ni == 1: jmax += 1
  if dir == 1:
    for k in xrange(kmin,kmax):
      for j in xrange(jmin,jmax):
        indf = (imin-1) + (j-1)*ni+(k-1)*ni*nj1
        bcfaces[nof] = indf+1; nof+=1
  elif dir == 2:
    for k in xrange(kmin,kmax):
      for i in xrange(imin,imax):
        indf = (i-1) + (jmin-1)*ni1+(k-1)*ni1*nj
        bcfaces[nof] = indf+1+ninti; nof+=1
  else:
    for j in xrange(jmin,jmax):
      for i in xrange(imin,imax):
        indf = (i-1) + (j-1)*ni1+(kmin-1)*ni1*nj1
        bcfaces[nof] = indf+1+ninti+nintj; nof+=1
  return bcfaces

# -- Computes the nearmatch ratio between two near-matching grids (structured)
def getNMRatio__(win, winopp, trirac):
  i1 = win[0]; j1 = win[2]; k1 = win[4]
  i2 = win[1]; j2 = win[3]; k2 = win[5]
  i1o = i1; j1o = j1; k1o = k1; i2o = i2; j2o = j2; k2o = k2
  i1opp = winopp[0]; j1opp = winopp[2]; k1opp = winopp[4]
  i2opp = winopp[1]; j2opp = winopp[3]; k2opp = winopp[5]
  oi = trirac[0];
  if len(trirac) > 1: oj = trirac[1]
  else: oj = 2
  if len(trirac) > 2: ok = trirac[2]
  else: ok = 3
  if oi == 2 or oi == -2: j1o = i1; j2o = i2
  elif oi == 3 or oi == -3: k1o = i1; k2o = i2

  if oj == 1 or oj == -1: i1o = j1; i2o = j2
  elif oj == 3 or oj == -3: k1o = j1; k2o = j2

  if ok == 1 or ok == -1: i1o = k1; i2o = k2
  elif ok == 2 or ok == -2: j1o = k1; j2o = k2

  diopp = max(1,i2opp-i1opp); dio = max(1,i2o-i1o)
  djopp = max(1,j2opp-j1opp); djo = max(1,j2o-j1o)
  dkopp = max(1,k2opp-k1opp); dko = max(1,k2o-k1o)

  ri = float(diopp)/dio
  rj = float(djopp)/djo
  rk = float(dkopp)/dko
  res = [1.,1.,1.]

  if   (oi == 2 or oi ==-2): res[0] = rj
  elif (oi == 3 or oi ==-3): res[0] = rk
  elif (oi == 1 or oi ==-1): res[0] = ri

  if   (oj == 1 or oj==-1): res[1] = ri
  elif (oj == 3 or oj==-3): res[1] = rk
  elif (oj == 2 or oj==-2): res[1] = rj

  if   (ok == 1 or ok ==-1): res[2] = ri
  elif (ok == 2 or ok ==-2): res[2] = rj
  elif (ok == 3 or ok ==-3): res[2] = rk
  return res

# -- recoverBCs
# Identifie des subzones comme BC Faces
# IN: liste de geometries de BC, liste de leur nom, liste de leur type
# OUT: a modifie avec des BCs ajoutees en BCFaces
def _recoverBCs(a, (BCs, BCNames, BCTypes), tol=1.e-11):
  try:import Post.PyTree as P
  except: raise ImportError("_recoverBCs: requires Post module.")
  _deleteZoneBC__(a)
  zones = Internal.getZones(a)
  for z in zones:
    indicesF = []
    f = P.exteriorFaces(z, indices=indicesF)
    indicesF = indicesF[0]
    hook = createHook(f, 'elementCenters')

    for c in xrange(len(BCs)):
      b = BCs[c]
      if b == []:
        raise ValueError("_recoverBCs: boundary is probably ill-defined.")
      # Break BC connectivity si necessaire
      elts = Internal.getElementNodes(b)
      size = 0
      for e in elts:
        erange = Internal.getNodeFromName1(e, 'ElementRange')[1]
        size += erange[1]-erange[0]+1
      n = len(elts)
      if n == 1 :
        ids = identifyElements(hook, b[0], tol)
      else:
        bb = breakConnectivity(b)
        ids = numpy.array([],dtype=numpy.int32)
        for bc in bb:
          ids = numpy.concatenate([ids, identifyElements(hook, bc, tol)])
      # Cree les BCs
      ids = ids[ids > -1]
      sizebc = ids.size
      if sizebc > 0:
        id2 = numpy.empty(sizebc, numpy.int32)
        id2[:] = indicesF[ids[:]-1]
        _addBC2Zone(z, BCNames[c], BCTypes[c], faceList=id2)
    freeHook(hook)

  return None

# -- pushBC
# Recopie les BCs de z1 sur z2 (zones a zones)
# par identification geometrique. Seules les BCs correspondant a des faces
# identiques geometriquement seront reportees.
# IN: z1: any type (STRUCT,BE,NGON)
# IN: z2 tout sauf STRUCT
# IN: type='F': pour tout type: output en face list
#     type='BCC': seulement si z2 est BE: output en BCC (BC connect)
# IN: overwriteBC: efface les BC de t2 si existent
def pushBC(t1, t2, type='F', tol=1.e-12, overwriteBC=True):
  """Put BCs of t1 to t2.
  Usage: t2 = pushBC(t1, t2)"""
  zones1 = Internal.getZones(t1)
  t2p = Internal.copyRef(t2)
  zones2 = Internal.getZones(t2p)
  nz = min(len(zones1),len(zones2))
  for c in xrange(nz):
    z1 = zones1[c]; zp = zones2[c]
    BCs = Internal.getNodesFromType2(z1, 'BC_t')
    if len(BCs) != 0:
      if overwriteBC: _deleteZoneBC__(zp)
      dims = Internal.getZoneDim(zp)
      outBCC = 0
      if (dims[0] == 'Unstructured' and dims[3] != 'NGON' and type != 'F'):
        outBCC = 1
      if outBCC == 0: KDT = createHook(zp, function='faceCenters')

      for b in BCs:
        name = b[0]
        BCType = Internal.getValue(b)
        if BCType == 'FamilySpecified':
          familyName = Internal.getNodeFromType1(b, 'FamilyName_t')
          if familyName is not None: BCType = 'FamilySpecified:%s'%Internal.getValue(familyName)
        ext = extractBCOfName(z1, name)
        if outBCC == 0:
          n = identifyElements(KDT, ext, tol)
          n = n[n>0]
          if n.size > 0: _addBC2Zone(zp, name, BCType, faceList=n)
        else: # output BCC
          try:
            import Transform.PyTree as T
            ext = T.breakElements(ext[0])
          except: pass
          valid = []
          for e in ext:
            d = Internal.getZoneDim(e)
            if d[0] == 'Structured': valid.append(convertArray2Hexa(e))
            elif (d[0] == 'Unstructured' and d[3] != 'NGON'): valid.append(e)
          if len(valid) == 1:
            _addBC2Zone(zp, name, BCType, subzone=valid[0])
          elif len(valid) > 1:
            f = 0
            for v in valid: _addBC2Zone(zp, name+str(f), BCType,
                                        subzone=v); f += 1
      if outBCC == 0: freeHook(KDT)
  return t2p

def identifyBC(t, infos, tol=1.e-12):
  try: import Post.PyTree as P
  except: raise ImportError("identifyBC: requires Post.PyTree module.")
  if infos == []: return t
  allWins = []
  for info in infos: allWins.append(Converter.node2Center(info[3]))
  # Creation d'un hook global a partir de toutes les fenetres
  indirWins = [] # permet d'identifier la fenetre a laquelle se rapporte un pt du hook
  globalHook, indirWins = Converter.createGlobalHook(allWins, 'nodes',
                                                     extended=1)
  # Identify and gather...
  tpp,typen = Internal.node2PyTree(t)
  for nob in xrange(len(tpp[2])):
    if tpp[2][nob][3] == 'CGNSBase_t':
      for noz in xrange(len(tpp[2][nob][2])):
        if tpp[2][nob][2][noz][3] == 'Zone_t':
          z = tpp[2][nob][2][noz]
          dimZ = Internal.getZoneDim(z)
          if dimZ[0] == 'Structured':
            niZ = dimZ[1]; njZ = dimZ[2]; nkZ = dimZ[3]
            faces = P.exteriorFacesStructured(z)
            dirf = 0
            for face in faces:
              dirf += 1
              dimW = Internal.getZoneDim(face)
              niw = dimW[1]; njw = dimW[2]
              if niw == 1 and njw == 1: pass
              else:
                face = node2Center(face)
                dimW = Internal.getZoneDim(face)
                niw = dimW[1]; njw = dimW[2]
                idHook = identifyNodes(globalHook, face, tol)
                if max(idHook) > -1:
                  ranges,nowins = Internal.gatherInStructPatch2D__(idHook,indirWins,niw,njw,dirf,niZ,njZ,nkZ)
                  if ranges != []:
                    nor = 0
                    for r in ranges:
                      noinfo = nowins[nor]
                      info = infos[noinfo]
                      bcname = info[0]
                      bctype = info[1]
                      bcdata = info[2]
                      win = info[3]
                      famName = ''
                      if bctype.split(':')[0] == 'FamilySpecified': famName=bctype.split(':')[1]

                      # Passage en noeuds des indices
                      istart=r[0]; iend=r[1]
                      jstart=r[2]; jend=r[3]
                      kstart=r[4]; kend=r[5]
                      if istart != iend:
                        if iend>1: iend=min(niZ,iend+1)
                      if jstart != jend:
                        if jend>1: jend=min(njZ,jend+1)
                      if kstart != kend:
                        if kend>1: kend=min(nkZ,kend+1)
                      else: # cas 2D
                        if nkZ == 2 and kend == 1: kend = 2
                      nor += 1

                      # addBC
                      z = addBC2Zone(z,bcname,bctype,[istart,iend,jstart,jend,kstart,kend],data=bcdata)
                      tpp[2][nob][2][noz] = z
  Converter.freeHook(globalHook)
  tp = Internal.pyTree2Node(tpp, typen)
  return tp

# -- tagDefinedBC
# tag des points definis par une CL:
# tag=0 si pas de BC
#    =1 point frontiere
#    =2 point interieur
def tagDefinedBC(t):
  a = Internal.copyRef(t)
  # si a = arbre
  toptree = Internal.isTopTree(a)
  if toptree:
    bases = Internal.getBases(a)
    for b in bases:
      zones = Internal.getNodesFromType1(b, 'Zone_t')
      for z in zones:
        dim = Internal.getZoneDim(z)[4]
        p = Internal.getParentOfNode(b, z)
        if dim == 3: b[2][p[1]] = tagDefinedBCForZone3D__(z)
        else: b[2][p[1]] = tagDefinedBCForZone2D__(z)
  else:
    # si a = base
    base = Internal.getBases(a)
    if base != []:
      zones = Internal.getZones(a)
      for z in zones:
        dim = Internal.getZoneDim(z)[4]
        p = Internal.getParentOfNode(a, z)
        if dim == 3: a[2][p[1]] = tagDefinedBCForZone3D__(z)
        else: a[2][p[1]] = tagDefinedBCForZone2D__(z)
    else:
      # liste zones ou zone ?
      stdNode = Internal.isStdNode(a)
      if stdNode == 0 : # liste de zones
        zones = Internal.getZones(a)
        nzones = len(zones)
        for noz in xrange(nzones):
          z = zones[noz]
          dim = Internal.getZoneDim(z)[4]
          if dim == 3: a[noz] = tagDefinedBCForZone3D__(z)
          else: zones[noz] = tagDefinedBCForZone2D__(z)
      else:
        dim = Internal.getZoneDim(a)[4]
        if dim == 3: a = tagDefinedBCForZone3D__(a)
        else: a = tagDefinedBCForZone2D__(a)

  return a

def tagDefinedBCForZone2D__(z):
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured': return z
    ni = dims[1]; nj = dims[2]; nk = dims[3]
    tag = Converter.array('definedBC',ni,nj,1)
    taga = tag[1][0]; wins = []
    # BC classique
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    for bc in bnds:
        range = Internal.getNodeFromName1(bc, 'PointRange')
        r = range[1]
        wins.append(Internal.range2Window(r))
    # BCNearMatch
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for bc in bnds:
      type = Internal.getNodeFromName1(bc, 'GridConnectivityType')
      if type is not None:
        val = Internal.getValue(type)
        if val == 'Abutting':
          range = Internal.getNodeFromName1(bc, 'PointRange')
          r = range[1]
          wins.append(Internal.range2Window(r))

    # BC match
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for bc in bnds:
        range = Internal.getNodeFromName1(bc, 'PointRange')
        r = range[1]
        wins.append(Internal.range2Window(r))
    # BC Overlap
    bnds = Internal.getNodesFromType1(z, 'GridConnectivity_t')
    for bc in bnds:
        r = Internal.getNodesFromName1(bc, 'GridConnectivityType')
        if r is not None:
            val = Internal.getValue(r)
            if val == 'Overset':
                dd = 0 # doubly defined
                userDef = Internal.getNodesFromName(bc, 'UserDefinedData')
                if userDef != []:
                    if len(userDef[0]) == 4:
                        info = userDef[0][2][0]
                        if info[0] == 'doubly_defined': dd = 1
                if dd == 0:
                    range = Internal.getNodeFromName(bc, 'PointRange')
                    wins.append(Internal.range2Window(range[1]))

    tag = Converter.converter.tagDefinedBC(tag, wins, 2)
    z = setFields([tag], z, 'nodes')
    return z

def tagDefinedBCForZone3D__(z):
  dims = Internal.getZoneDim(z)
  if dims[0] == 'Unstructured': return z
  ni = dims[1]; nj = dims[2]; nk = dims[3]

  tag = Converter.array('definedBC', ni, nj, nk)
  wins = []
  # BC classique
  bnds = Internal.getNodesFromType2(z, 'BC_t')
  for bc in bnds:
      range = Internal.getNodeFromName1(bc, 'PointRange')
      r = range[1]
      wins.append(Internal.range2Window(r))
  # BCNearMatch
  bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
  for bc in bnds:
    type = Internal.getNodeFromName1(bc, 'GridConnectivityType')
    if type is not None:
      val = Internal.getValue(type)
      if val == 'Abutting':
        range = Internal.getNodeFromName1(bc, 'PointRange')
        r = range[1]
        wins.append(Internal.range2Window(r))
  # BC match
  bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
  for bc in bnds:
      range = Internal.getNodeFromName1(bc, 'PointRange')
      r = range[1]
      wins.append(Internal.range2Window(r))
  # BC Overlap
  bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
  for bc in bnds:
      r = Internal.getNodeFromName1(bc, 'GridConnectivityType')
      if r is not None:
          val = Internal.getValue(r)
          if val == 'Overset':
              dd = 0 # doubly defined
              userDef = Internal.getNodesFromName(bc, 'UserDefinedData')
              if userDef != []:
                  if len(userDef[0]) == 4:
                      info = userDef[0][2][0]
                      if info[0] == 'doubly_defined': dd = 1
              if dd == 0:
                  range = Internal.getNodeFromName1(bc, 'PointRange')
                  wins.append(Internal.range2Window(range[1]))
  tag = Converter.converter.tagDefinedBC(tag, wins, 3)
  z = setFields([tag], z, 'nodes')
  return z

# -- updateDefinedBCForWins__
# update the definedBC field for a list of windows wins defined by BC
# wins are [i1,i2,j1,j2,k1,k2]
def updateDefinedBCForWins__(z, wins, dim=3):
  tag = getField('definedBC', z)[0]
  for now in xrange(len(wins)):
    if len(wins[now]) < 6:
      for i in xrange(len(wins[now]),7): wins[now].append(1)
  tag = Converter.converter.tagDefinedBC(tag, wins, dim)
  z = setFields([tag], z, 'nodes')
  return z

#=============================================================================
# -- Remove BCs of type or name from a tree or a zone --
#=============================================================================

# -- rmBCOfType
def rmBCOfType(t, bndType):
  """Remove BCs of given type.
  Usage: rmBCOfType(t, bndType)"""
  tp = Internal.copyRef(t)
  _rmBCOfType(tp, bndType)
  return tp

def _rmBCOfType(t, bndType):
  names = bndType.split(':')
  if len(names) == 2: # Family
    _rmBCOfName(t, bndType)

  zones = Internal.getZones(t)
  if bndType == 'BCMatch':
    for z in zones:
      nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
      for i in nodes:
        (parent, d) = Internal.getParentOfNode(z, i)
        del parent[2][d]
      nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
      for i in nodes:
        r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
        if r is not None:
          val = Internal.getValue(r)
          if val == 'Abutting1to1':
            (parent, d) = Internal.getParentOfNode(z, i)
            del parent[2][d]
  elif bndType == 'BCNearMatch' or bndType == 'BCStage':
    for z in zones:
      nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
      for i in nodes:
        r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
        if r is not None:
          val = Internal.getValue(r)
          if val == 'Abutting':
            (parent, d) = Internal.getParentOfNode(z, i)
            del parent[2][d]
  elif bndType == 'BCOverlap':
    for z in zones:
      nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
      for i in nodes:
        r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
        if r is not None:
          val = Internal.getValue(r)
          if val == 'Overset':
            (parent, d) = Internal.getParentOfNode(z, i)
            del parent[2][d]
  else:
    families = getFamilyBCNamesOfType(t, bndType)
    if bndType == 'BCWall':
      familiesI = getFamilyBCNamesOfType(t, 'BCWallInviscid')
      familiesV = getFamilyBCNamesOfType(t, 'BCWallViscous*')
    for z in zones:
      nodes = Internal.getNodesFromValue(z, bndType)
      nodes += getFamilyBCs(z, families)
      if bndType == 'BCWall':
        nodes += Internal.getNodesFromValue(z, 'BCWallInviscid')
        nodes += Internal.getNodesFromValue(z, 'BCWallViscous*')
        nodes += getFamilyBCs(z, familiesI)
        nodes += getFamilyBCs(z, familiesV)
      for i in nodes:
        (parent, d) = Internal.getParentOfNode(z, i)
        del parent[2][d]
  return None

# -- rmBCOfName
def rmBCOfName(t, bndName):
  """Remove BCs of given name.
  Usage: rmBCOfName(t, bndName)"""
  tp = Internal.copyRef(t)
  _rmBCOfName(tp, bndName)
  return tp

def _rmBCOfName(t, bndName):
  names = bndName.split(':')
  if len(names) == 1: # real bnd name
    zones = Internal.getZones(t)
    for z in zones:
      nodes = Internal.getNodesFromName(z, bndName)
      for i in nodes:
        if Internal.getType(i) in ['BC_t', 'GridConnectivity1to1_t', 'GridConnectivity_t']:
          (parent, d) = Internal.getParentOfNode(z, i)
          del parent[2][d]
  else: # family specified BC
    zones = Internal.getZones(t)
    for z in zones:
      nodes = getFamilyBCs(z, names[1])
      for i in nodes:
        (parent, d) = Internal.getParentOfNode(z, i)
        del parent[2][d]
  return None

#==============================================================================
# -- Get empty BCs from a zone --
#==============================================================================

# -- getEmptyBC
# Return the list of empty BCs:
# return range (structured) or face list (unstructured) of undefined BCs
# for any zone in any basis
# if t=tree: returns [winsBase1, winsBase2,...],
# with winsBase=[winzone1Base1, winzone2Base1,...]
# if all the BCs are defined for a zone, [] is returned for the zone
# IN: dim: dimension of the pb
# IN: splitFactor: used only for unstructured grids, split the windows
# following split angle.
def getEmptyBC(a, dim=3, splitFactor=181.):
  """Return the range or facelist of unset boundary conditions."""
  # si a = arbre
  toptree = Internal.isTopTree(a)
  if toptree:
    winst = []
    bases = Internal.getBases(a)
    for b in bases:
      winsb = []
      zones = Internal.getNodesFromType1(b, 'Zone_t')
      for z in zones:
        winsz = getEmptyBCForZone__(z, dim, splitFactor)
        winsb.append(winsz)
      winst.append(winsb)
    return winst
  else:
    # si a = base
    base = Internal.getBases(a)
    if base != []:
      winsb = []
      zones = Internal.getNodesFromType1(base, 'Zone_t')
      for z in zones:
        winsz = getEmptyBCForZone__(z, dim, splitFactor)
        winsb.append(winsz)
      return winsb
    else:
      # liste zones ou zone ?
      stdNode = Internal.isStdNode(a)
      if stdNode == 0: # liste de zones
        wins = []
        for z in a[0:]:
          winsz = getEmptyBCForZone__(z, dim, splitFactor)
          wins.append(winsz)
        return wins
      else:
        return getEmptyBCForZone__(a, dim, splitFactor)
  return []

# -- Detect empty boundary conditions for zones
def getEmptyBCForZone__(z, pbDim, splitFactor):
  dims = Internal.getZoneDim(z)
  if dims[0] == 'Structured': 
    return getEmptyBCForStructZone__(z, dims, pbDim, splitFactor)
  elif dims[3] == 'NGON': 
    return getEmptyBCForNGonZone__(z, dims, pbDim, splitFactor)
  else: return getEmptyBCForBEZone__(z, dims, pbDim, splitFactor)
  
# -- Detect empty boundary conditions for struct zones
# Renvoie une liste de ranges des BC empty
def getEmptyBCForStructZone__(z, dims, pbDim, splitFactor):
  ni = dims[1]; nj = dims[2]; nk = dims[3]
  ni1 = ni-1; nj1 = nj-1; nk1 = nk-1
  nwins = []; wins = []

  # BC classique
  bnds = Internal.getNodesFromType2(z, 'BC_t')
  for bc in bnds:
    r = Internal.getNodeFromName1(bc, 'PointRange')
    if r is not None: wins.append(Internal.range2Window(r[1]))
  # BC match
  bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
  for bc in bnds:
    r = Internal.getNodeFromName1(bc, 'PointRange')
    if r is not None: wins.append(Internal.range2Window(r[1]))

  # BC Overlap/NearMatch/NoMatch
  bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
  for bc in bnds:
    r = Internal.getNodeFromName1(bc, 'PointRange')
    if r is not None: wins.append(Internal.range2Window(r[1]))

  # Parcours des faces
  directions = []
  if ni != 1: directions += [1,2]
  if nj != 1: directions += [3,4]
  if nk != 1 and pbDim == 3: directions += [5,6]
  for dir in directions:
    nwins = Converter.converter.detectEmptyBC(wins, ni, nj, nk, dir, nwins)
  return nwins

# -- Detect empty boundary conditions for unstruct zones
# Renvoie une liste des indices des faces des BC empty
def getEmptyBCForBEZone__(z, dims, pbDim, splitFactor):
  try: import Transform.PyTree as T; import Post.PyTree as P
  except: raise ImportError("getEmptyBC: requires Transform and Post modules.")
  indicesF = []
  f = P.exteriorFaces(z, indices=indicesF)
  indicesF = indicesF[0]
  _initVars(f, 'centers:__tag__', 1)

  # BC classique
  bnds = Internal.getNodesFromType2(z, 'BC_t')
  # BC Match
  bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
  # BC Overlap/NearMatch/NoMatch
  bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')

  zp = Internal.copyRef(z)
  _deleteZoneBC__(zp); _deleteFlowSolutions__(zp)

  defined = [] # BC deja definies
  for bc in bnds:
    flist = Internal.getNodeFromName1(bc, Internal.__FACELIST__)
    if flist is not None: defined.append(T.subzone(zp, flist[1], type='faces'))
    erange = Internal.getNodeFromName1(bc, Internal.__ELEMENTRANGE__)
    if erange is not None:
      r = erange[1]
      defined.append(selectOneConnectivity(zp, range=[r[0,0],r[0,1]]))

  hook = createHook(f, 'elementCenters')
  if defined != []:
    tag = Internal.getNodeFromName2(f, '__tag__')[1]
    for d in defined:
      d = convertArray2NGon(d)
      id0 = identifyElements(hook, d)
      tag[id0[:]-1] = 0

  sel = P.selectCells2(f, 'centers:__tag__')
  if splitFactor >= 180.: sel = T.splitConnexity(sel)
  else: sel = T.splitSharpEdges(sel, alphaRef=splitFactor)

  id0 = []
  for s in sel:
    id1 = identifyElements(hook, s)
    #mask = (id1[:] >= 0) # enleve les elements non identifies
    id2 = numpy.empty(id1.size, numpy.int32)
    id2[:] = indicesF[id1[:]-1]
    id0.append(id2)
  freeHook(hook)
  return id0

def getEmptyBCForNGonZone__(z, dims, pbDim, splitFactor):
  try: import Transform.PyTree as T; import Post.PyTree as P
  except: raise ImportError("getEmptyBC: requires Transform and Post modules.")
  indicesF = []
  f = P.exteriorFaces(z, indices=indicesF)
  indicesF = indicesF[0]
  # BC classique
  bnds = Internal.getNodesFromType2(z, 'BC_t')
  # BC Match
  bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
  # BC Overlap/NearMatch/NoMatch
  bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')

  indicesBC = []
  for b in bnds:
    f = Internal.getNodeFromName1(b, 'PointList')
    indicesBC.append(f[1])

  undefBC = False
  if indicesBC != []:
    indicesBC = numpy.concatenate(indicesBC, axis=1)
    nfacesExt = indicesF.shape[0]
    nfacesDef = indicesBC.shape[1]
    if nfacesExt<nfacesDef:
      print 'Warning: zone %s: number of faces defined by BCs is greater than the number of external faces. Try to reduce the matching tolerance.'%(z[0])
    elif nfacesExt>nfacesDef:
      indicesE = Converter.converter.diffIndex(indicesF, indicesBC)
      undefBC = True
  else: 
    undefBC = True
    indicesE = indicesF

  if undefBC:
    ze = T.subzone(z, indicesE, type='faces')
    zea = getFields('GridCoordinates', ze)[0]
    id0 = T.transform.splitSharpEdgesList(zea, indicesE, splitFactor)
    return id0
  else: return []

#==============================================================================
# IN: z: zone subzone (structure)
# IN: w: window [imin,imax, jmin, jmax, kmin, kmax]
# Reorder unz subzone pour garder les normales exterieures
#==============================================================================
def _reorderSubzone__(z, w, T):
  dim = Internal.getZoneDim(z)
  ni = dim[1]; nj = dim[2]; nk = dim[3]
  imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
  if imin == imax and imin == ni: T._reorder(z, (1,-2,3))
  elif jmin == jmax and jmin == 1: T._reorder(z, (-1,2,3))
  elif kmin == kmax and kmin == nk: T.reorder(z, (-1,2,3))
  return None

#==============================================================================
# -- Extract all BC of type or name --
#==============================================================================

# Get the geometrical BC and append it to res
# IN: i: BC_t node or GridConnectivity_t node
# IN: z: zone owning BC
# IN: T: transform module
def getBC__(i, z, T, res):
  # IndexRange
  r = Internal.getNodeFromType1(i, 'IndexRange_t')
  if r is not None and r[1].shape[0] > 1: # structure - suppose range in nodes
    wrange = r[1]
    w = Internal.range2Window(wrange)
    #zp = subzoneWithReorder__(z, w)
    imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
    zp = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
    zp[0] = z[0]+Internal.SEP1+i[0]
    # Get BCDataSet if any
    datas = Internal.getBCDataSet(z, i)
    if datas != []:
      f = Internal.createUniqueChild(zp, Internal.__FlowSolutionCenters__,
                                     'FlowSolution_t')
      Internal.newGridLocation(value='CellCenter', parent=f)
      for d in datas: Internal.createUniqueChild(f, d[0], d[3], value=d[1])
    _reorderSubzone__(zp, w, T) # normales ext
    _deleteZoneBC__(zp)
    _deleteGridConnectivity__(zp)
    res.append(zp)
  elif r is not None and r[1].shape[0] == 1: # suppose BE + BCC
    r = r[1]
    # zp = selectOneConnectivity(z, range=[r[0,0],r[0,1]])
    zp = selectConnectivity(z, range=[r[0,0],r[0,1]])
    zp[0] = z[0]+Internal.SEP1+i[0]
    _deleteZoneBC__(zp)
    _deleteGridConnectivity__(zp)
    _deleteSolverNodes__(zp)
    res.append(zp)
  # IndexArray
  if r is None: r = Internal.getNodeFromName(i, Internal.__FACELIST__)
  else: r = None
  if r is not None:
    loc = Internal.getNodeFromName1(i, 'GridLocation')
    if loc is not None:
      val = Internal.getValue(loc)
      if val == 'FaceCenter': # Face list (BE ou NGON)
        faceList = r[1]
        rf = Internal.getElementRange(z, type='NGON')
        if rf is not None and rf[0] != 1: # decalage possible du NGON
          faceList2 = numpy.copy(faceList)
          faceList2[:] = faceList[:]-rf[0]+1
          zp = T.subzone(z, faceList2, type='faces')
        else: zp = T.subzone(z, faceList, type='faces')
        zp[0] = z[0]+Internal.SEP1+i[0]
        _deleteZoneBC__(zp)
        _deleteGridConnectivity__(zp)
        _deleteSolverNodes__(zp)
        # Get BCDataSet if any
        datas = Internal.getBCDataSet(z, i)
        if datas != []:
          f = Internal.createUniqueChild(zp, Internal.__FlowSolutionCenters__,
                                         'FlowSolution_t')
          Internal.newGridLocation(value='CellCenter', parent=f)
          for d in datas: Internal.createUniqueChild(f, d[0], d[3], value=d[1])
        res.append(zp)
    else: # suppose FaceList
      faceList = r[1]
      rf = Internal.getElementRange(z, type='NGON')
      if (rf is not None and rf[0] != 1):
        faceList2 = numpy.copy(faceList)
        faceList2[:] = faceList[:]-rf[0]+1
        zp = T.subzone(z, faceList2, type='faces')
      else: zp = T.subzone(z, faceList, type='faces')
      zp[0] = z[0]+Internal.SEP1+i[0]
      _deleteZoneBC__(zp)
      _deleteGridConnectivity__(zp)
      _deleteSolverNodes__(zp)
      # Get BCDataSet if any
      datas = Internal.getBCDataSet(z, i)
      if datas != []:
        f = Internal.createUniqueChild(zp, __FlowSolutionCenters__,
                                       'FlowSolution_t')
        Internal.newGridLocation(value='CellCenter', parent=f)
        for d in datas: Internal.createUniqueChild(f, d[0], d[3], value=d[1])
      res.append(zp)

# -- extractBCOfType
# Extract all BC of a given type
# Recognised bndType: classic CGNS + BCMatch + BCNearMatch + BCOverlap
# topTree: utile si t n'est pas un topTree, permet de trouver les familyBCs
def extractBCOfType(t, bndType, topTree=None):
    """Extract the grid coordinates of given BC type as zones."""
    try: import Transform.PyTree as T
    except:
        raise ImportError("extractBCOfType: requires Transform.PyTree module.")

    names = bndType.split(':')
    if len(names) == 2: # Family
      return extractBCOfName(t, bndType)

    zones = Internal.getZones(t)
    res = []
    if bndType == 'BCMatch':
        for z in zones:
            # Cherche GridConnectivity1to1_t
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for i in nodes: getBC__(i, z, T, res)
            # Cherche GridConnectivity_t + Abutting1to1
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
              type = Internal.getNodeFromName1(i, 'GridConnectivityType')
              if type is not None:
                val = Internal.getValue(type)
                if val == 'Abutting1to1': getBC__(i, z, T, res)
    elif bndType == 'BCNearMatch' or bndType == 'BCStage':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
              type = Internal.getNodeFromName1(i, 'GridConnectivityType')
              if type is not None:
                val = Internal.getValue(type)
                if val == 'Abutting': getBC__(i, z, T, res)
    elif bndType == 'BCOverlap':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                r = Internal.getNodeFromName1(i, 'GridConnectivityType')
                if r is not None:
                  val = Internal.getValue(r)
                  if val == 'Overset': getBC__(i, z, T, res)
    else: # BC physiques
        if topTree is None: topTree = t
        families = getFamilyBCNamesOfType(topTree, bndType)
        if bndType == 'BCWall':
          families1 = getFamilyBCNamesOfType(topTree, 'BCWallInviscid')
          families2 = getFamilyBCNamesOfType(topTree, 'BCWallViscous*')
        for z in zones:
            nodes = Internal.getNodesFromValue(z, bndType)
            nodes += getFamilyBCs(z, families)
            if bndType == 'BCWall':
              nodes += Internal.getNodesFromValue(z, 'BCWallInviscid')
              nodes += Internal.getNodesFromValue(z, 'BCWallViscous*')
              nodes += getFamilyBCs(z, families1)
              nodes += getFamilyBCs(z, families2)              
            for i in nodes: getBC__(i, z, T, res)
    return res

# -- extractBCOfName
def extractBCOfName(t, bndName):
  """Extract the grid coordinates of given BC name as zones.
  Usage: extractBCOfName(t, bndName)"""
  try: import Transform.PyTree as T
  except:
    raise ImportError("extractBCOfName: requires Transform.PyTree module.")
  names = bndName.split(':')
  res = []
  if len(names) == 1: # real bnd name
    for z in Internal.getZones(t):
      # Pas de niveau 2 car pas de wild card autorisee dans ce cas
      nodes = Internal.getNodesFromName(z, bndName)
      for i in nodes:
        if Internal.getType(i) in ['BC_t', 'GridConnectivity1to1_t', 'GridConnectivity_t']:
          getBC__(i, z, T, res)
  else: # family specified BC
    for z in Internal.getZones(t):
      nodes = getFamilyBCs(z, names[1])
      for i in nodes: getBC__(i, z, T, res)
  return res

# -- getBCs
# Retourne la geometrie de toutes les BCs
# IN: t: une zone, une liste de zones, une base ou un arbre
# OUT: BCs: liste des geometries de bcs
# OUT: BCNames: liste des noms des BCs
# OUT: BCTypes: liste des types des BCs
def getBCs(t):
  BCs = []; BCNames = []; BCTypes = []
  for z in Internal.getZones(t):
    nodes = Internal.getNodesFromType2(z, 'BC_t')
    for n in nodes:
      name = n[0]; typeBC = Internal.getValue(n)
      if typeBC == 'FamilySpecified':
        fname = Internal.getNodeFromType1(n,'FamilyName_t')
        fname = Internal.getValue(fname)
        typeBC = 'FamilySpecified:%s'%fname
      zBC = extractBCOfName(z, name)
      if zBC == []:
        name2 = name.split(':')
        if len(name2)>1 and name2[0] == 'FamilySpecified':
          raise ValueError("getBCs: BC of name FamilySpecified:* is not valid.")
      BCs.append(zBC); BCNames.append(name); BCTypes.append(typeBC)

    # Raccords Stage* definis par des familles
    nodes = Internal.getNodesFromType2(z,"ZoneGridConnectivity_t")
    for n in nodes:
      for gc in Internal.getNodesFromType1(n,"GridConnectivity_t"):
        name = Internal.getName(gc)
        fname = Internal.getNodeFromType1(gc,'FamilyName_t') 
        if fname is not None:
          fname = Internal.getValue(fname)
          zBC = extractBCOfName(z,name)
          typeGC = 'FamilySpecified:%s'%fname
          BCs.append(zBC); BCNames.append(name); BCTypes.append(typeGC)
  
  return (BCs, BCNames, BCTypes)

# -- extractBCDataStruct__
# Extract BC characteristics of structured zones [Name,BCType,Data,coords]
# Name is the name of the BC
# Type is the type of BC
# Data is the data attached to some BCs, [] if not used
# coords: coordinates of 9 particular points defining the zone
# BCMatch and BCNearMatch are not taken into account...
def extractBCDataStruct__(z):
  import Transform as T
  infos = []
  dims = Internal.getZoneDim(z)
  ni = dims[1]; nj = dims[2]; nk = dims[3]

  # BCOverlap
  nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
  for i in nodes:
    bcname = i[0]
    bctype = Internal.getValue(i)
    bcdata = Internal.getNodesFromType1(i, 'BCDataSet_t')
    if bcdata != []: bcdata = bcdata[0]
    else: bcdata = None
    fname = Internal.getNodeFromType1(i, 'FamilyName_t')
    if fname is not None: fname = Internal.getValue(fname)

    r = Internal.getNodeFromName1(i, 'GridConnectivityType')
    if r is not None:
      val = Internal.getValue(r)
      if val == 'Overset':
        pr = Internal.getNodeFromName1(i, 'PointRange')
        if pr is not None:
          if fname is not None: bctype = 'FamilySpecified:'+fname
          range0 = pr[1]
          w = Internal.range2Window(range0)
          imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
          coords = getFields(Internal.__GridCoordinates__,z)[0]
          coords = T.subzone(coords, (imin,jmin,kmin), (imax,jmax,kmax))
          info = [bcname, 'BCOverlap', bcdata, coords]
          infos.append(info)

  # All classical BC
  nodes = Internal.getNodesFromType2(z, 'BC_t')
  for i in nodes:
    bcname = i[0]
    bctype = Internal.getValue(i)
    bcdata = Internal.getNodesFromType1(i, 'BCDataSet_t')
    if bcdata != []: bcdata = bcdata[0]
    else: bcdata = None
    fname = Internal.getNodeFromType1(i, 'FamilyName_t')
    if fname is not None: fname = Internal.getValue(fname)
    pr = Internal.getNodeFromName1(i, 'PointRange')
    if pr is not None:
      if fname is not None: bctype = 'FamilySpecified:'+fname
      range0 = pr[1]
      w = Internal.range2Window(range0)
      imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
      coords = getFields(Internal.__GridCoordinates__,z)[0]
      coords = T.subzone(coords, (imin,jmin,kmin), (imax,jmax,kmax))
      info = [bcname, bctype, bcdata, coords]
      infos.append(info)
  return infos

def extractBCInfo(t):
  infos = []
  for z in Internal.getZones(t):
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Structured': infos += extractBCDataStruct__(z)
  return infos

#=============================================================================
# -- Fill empty BCs avec une CL donnee --
#=============================================================================
def fillEmptyBCWith(t, bndName, bndType, dim=3):
  """Fill empty BCs with given type."""
  a = Internal.copyRef(t)
  _fillEmptyBCWith(a, bndName, bndType, dim)
  return a

def _fillEmptyBCWith(t, bndName, bndType, dim=3):
  for z in Internal.getZones(t):
    c = 1
    wins = getEmptyBC(z, dim)
    for w in wins:
      ztype = Internal.getZoneType(z)
      if ztype == 1: # structured
        _addBC2Zone(z, bndName+str(c), bndType, w); c += 1
      else:
        _addBC2Zone(z, bndName+str(c), bndType, faceList=w); c += 1
  return None

#==============================================================================
# -- Family management --
#==============================================================================

# -- tagWithFamily (familyZone)
def tagWithFamily(z, familyName, add=False):
  """Tag a zone node or a BC node with a familyName.
  Usage: tagWithFamily(z, familyName, add)"""
  a = Internal.copyRef(z)
  _tagWithFamily(a, familyName, add)
  return a

def _tagWithFamily__(a, familyName, add=False):
  if a[3] != 'Zone_t' and a[3] != 'BC_t':
    print 'Warning: tagWithFamily: must be used on a Zone_t or BC_t node.'
  if not add:
    Internal._createUniqueChild(a, 'FamilyName', 'FamilyName_t', value=familyName)
  else:
    if Internal.getNodeFromType1(a, 'FamilyName_t') is None:
      Internal._createChild(a, 'FamilyName', 'FamilyName_t', value=familyName)
    else:
      Internal._createChild(a, 'AdditionalFamilyName', 'AdditionalFamilyName_t', value=familyName)
  return None

def _tagWithFamily(a, familyName, add=False):
  zones = Internal.getZones(a)
  if len(zones) > 0: # family of zones
     for z in zones: _tagWithFamily__(z, familyName, add)
  else:
     nodes = Internal.getNodesFromType(a, 'BC_t')
     for n in nodes: _tagWithFamily__(n, familyName, add)

# -- getFamilyZones (wildcard possible on familyName)
def getFamilyZones(t, familyName):
  """Return all zones that have this familyName.
  Usage: getFamilyZones(t, familyName)"""
  out = []
  if isinstance(familyName, str): families = [familyName]
  else: families = familyName

  for z in Internal.getZones(t):
    res = Internal.getNodesFromType1(z, 'FamilyName_t')
    for i in res:
      val = Internal.getValue(i)
      for f in families:
        if ('*' in f)|('?' in f)|('[' in f):
          if fnmatch.fnmatch(val, f): out.append(z)
        else:
          if val == f: out.append(z)
  return out

# -- getFamilyBCs (wildcards possible on familyName)
def getFamilyBCs(t, familyName):
  """Return all BC nodes that have this familyName.
  Usage: getFamilyBCs(t, familyName)"""
  out = []
  if isinstance(familyName, str): families = [familyName]
  else: families = familyName
  for z in Internal.getZones(t):
    nodes = Internal.getNodesFromType2(z, 'BC_t')
    nodes += Internal.getNodesFromType2(z, 'GridConnectivity_t')
    nodes += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for n in nodes:
      res = Internal.getNodesFromType1(n, 'FamilyName_t')
      for i in res:
        val = Internal.getValue(i)
        val = val.strip()
        for f in families:
          if ('*' in f)|('?' in f)|('[' in f):
            if fnmatch.fnmatch(val, f): out.append(n)
          else:
            if val == f: out.append(n)
  return out

# -- getFamilyBCNamesOfType (wildcards possible on bndType)
def getFamilyBCNamesOfType(t, bndType=None):
  """Return the family BC names of a given type.
  Usage: names = getFamilyBCNamesOfType(t, 'BCWall')"""
  out = set()
  nodes = Internal.getNodesFromType2(t, 'Family_t')
  if bndType is None:
    for n in nodes:
      p = Internal.getNodeFromType1(n, 'FamilyBC_t')
      if p is not None: out.add(n[0])
  else:
    for n in nodes:
      p = Internal.getNodeFromType1(n, 'FamilyBC_t')
      if p is not None:
        bctype = Internal.getValue(p)
        if bctype == bndType: out.add(n[0])
  return list(out)

# -- getFamilyBCNamesDict
def getFamilyBCNamesDict(t):
  """Return the dictionary of familyBCs. d['CARTER'] = "BCWall"
  Usage: dict = getFamilyBCNamesDict(t)"""
  d = {}
  nodes = Internal.getNodesFromType2(t, 'Family_t')
  for n in nodes:
    btype = Internal.getNodeFromType1(n, 'FamilyBC_t')
    if btype is not None: 
      d[n[0]] = Internal.getValue(btype)
  return d
  
# -- getFamilyZoneNames
def getFamilyZoneNames(t):
  """Return the family zone names in a tree or a base."""
  out = set()
  nodes = Internal.getNodesFromType2(t, 'Family_t')
  for n in nodes:
    ret = Internal.getNodeFromType(n, 'FamilyBC_t')
    if ret is None: out.add(n[0])
  return list(out)

# -- Add a Family node to a base node
def addFamily2Base(base, familyName, bndType=None):
  """Add a family node to a base node.
  Usage: addFamily2Base(base, familyName, bndType)"""
  a = Internal.copyRef(base)
  _addFamily2Base(a, familyName, bndType)
  return a

def _addFamily2Base(a, familyName, bndType=None):
  bases = Internal.getBases(a)
  for b in bases:
    res = Internal.getNodeFromName1(b, familyName)
    child = Internal.createUniqueChild(b, familyName, 'Family_t', None)

    if bndType == 'BCOverlap': # Cas special
      Internal.createUniqueChild(child, 'FamilyBC', 'FamilyBC_t', 'UserDefined')
      Internal.createUniqueChild(child, '.Solver#Overlap', 'UserDefinedData_t', None)
    elif bndType == 'Mask': # Cas special
      Internal.createUniqueChild(child, 'FamilyBC', 'FamilyBC_t', bndType)
      mask = Internal.createUniqueChild(child, '.Solver#Mask', 'DataArray_t', None)
      Internal.createUniqueChild(mask, 'type', 'DataArray_t', 'xray')
      Internal.createUniqueChild(mask, 'dim1', 'DataArray_t', 1000)
      Internal.createUniqueChild(mask, 'dim2', 'DataArray_t', 1000)
      Internal.createUniqueChild(mask, 'delta', 'DataArray_t', 0.1)
    elif bndType is not None:
      Internal.createUniqueChild(child, 'FamilyBC', 'FamilyBC_t', bndType)
  return None

#==============================================================================
# -- change of localisation --
#==============================================================================

# -- node2Center
def node2Center(t, var=''):
    """Convert zone defined at nodes to a zone defined at centers.
    Usage: node2Center(t, varname)"""
    a = Internal.copyRef(t)
    if var == '':
      if Internal.isStdNode(a) == 0: la = a
      else: la = [a] # noeud standard => liste de noeuds standards
      for i in xrange(len(la)):
        zones = Internal.getZones(la[i])
        for z in zones:
            (p, pos) = Internal.getParentOfNode(la[i], z)
            fieldc = getFields(Internal.__FlowSolutionCenters__, z)
            _deleteFlowSolutions__(z, 'centers')
            _deleteZoneBC__(z)
            _deleteGridConnectivity__(z)
            dim = Internal.getZoneDim(z)[0]
            if (dim == 'Unstructured'):
              f = fieldc[0]
              if f != []: f[3] = f[3].replace('*', '')
              try: import Transform
              except: pass
              else:
                z = TZA(z, 'nodes', 'nodes', Transform.dual, 0, 0)
            else:
              z = TZA(z, 'nodes', 'nodes', Converter.node2Center, None)
            setFields(fieldc, z, 'nodes', writeDim=False)
            if p is None: la[i] = z
            else: p[2][pos] = z
      if Internal.isStdNode(a) == 0: a = la
      else: a = la[0] # liste de noeuds standards => noeud standard
      return a

    elif var == Internal.__GridCoordinates__:
      fieldc = []
      fieldn = getFields(Internal.__GridCoordinates__, a)
      for i in fieldn:
        if i != []:
          b = Converter.node2Center(i)
          fieldc.append(b)
        else:
          fieldc.append([])
      setFields(fieldc, a, 'centers')
      return a
    elif (var == Internal.__FlowSolutionNodes__):
      fieldc = []
      fieldn = getFields(Internal.__FlowSolutionNodes__, a)
      for i in fieldn:
        if i != []:
          b = Converter.node2Center(i)
          fieldc.append(b)
        else:
          fieldc.append([])
      setFields(fieldc, a, 'centers')
      return a
    else:
      if isinstance(var, list):
        for v in var: node2Center__(a, v)
      else: node2Center__(a, var)
      return a

def node2Center__(a, var):
  var, loc = Internal.fixVarName(var)
  if loc == 1: return a
  fieldn = getField(var, a)
  fieldc = Converter.node2Center(fieldn)
  setFields(fieldc, a, 'centers')
  return a

# Adapte les arrays pour les cas nk=1 => passe en nk=2 si la zone est nk=2
def _patchArrayForCenter2NodeNK1__(fields, a):
  zones = Internal.getZones(a)
  c = 0
  for f in fields:
    if f != []:
      z = zones[c]
      dim = Internal.getZoneDim(z)
      ni = dim[1]; nj = dim[2]; nk = dim[3]
      fp = f[1]
      nfld = fp.shape[0]; s = fp.shape[1]
      if (dim[0] == 'Structured' and nk == 2 and s == ni*nj):
        b = numpy.empty( (nfld, ni*nj*2) )
        b[:,0:ni*nj] = fp[:,0:ni*nj]
        b[:,ni*nj:2*ni*nj] = fp[:,0:ni*nj]
        f[1] = b; f[4] = 2
    c += 1

# -- center2Node
# Convert a zone defining centers to nodes or convert a field in a
# base/tree/zone located at centers to nodes
def center2Node(t, var=None, cellNType=0):
    """Converts array defined at centers to an array defined at nodes.
    Usage: center2Node(t, var, cellNType)"""

# Preparation pour le topTree
#     istoptree = Internal.isTopTree(t)
#     if (istoptree == False and topTree != []):
#       zones = C.getConnectedZones(t, topTree=topTree)
#       tp = C.newPyTree(['Base'])
#       zonest = Internal.getZones(t)
#       tp[2][1][2] += zonest; tp[2][1][2] += zones
#       tp = C.center2Node(tp, Internal.__FlowSolutionCenters__)
#       zone = tp[2][1][2][0]

    ghost = Internal.getNodeFromName(t, 'ZoneRind')
    if var is None:
      # solution en centres
      res = Internal.getNodesFromName3(t, Internal.__FlowSolutionCenters__)
      fieldsc = []
      if res != []:
        if Internal.getNodesFromType1(res[0], 'DataArray_t') != []:
          if ghost is None:
            a = Internal.addGhostCells(t, t, 1, modified=[Internal.__FlowSolutionCenters__])
          else: a = Internal.copyRef(t)
          fieldc = getFields(Internal.__FlowSolutionCenters__, a)
          fieldn = []
          listVar = []
          for i in fieldc:
            if i != []:
              b = Converter.center2Node(i, cellNType); fieldn.append(b)
              for va in b[0].split(','):
                if va not in listVar: listVar.append(va)
            else: fieldn.append([])
          _patchArrayForCenter2NodeNK1__(fieldn, a)
          setFields(fieldn, a, 'nodes', writeDim=False)
          if ghost is None:
            a = Internal.rmGhostCells(a, a, 1,
                                      modified=[listVar, Internal.__FlowSolutionCenters__])
          fieldsc = getFields(Internal.__FlowSolutionNodes__, a)

      # destruction
      t = deleteFlowSolutions__(t, 'centers')
      t = TZA(t, 'nodes', 'nodes', Converter.center2Node, cellNType)
      if fieldsc != []: setFields(fieldsc, t, 'centers', writeDim=False)
      return t

    elif var == Internal.__FlowSolutionCenters__:
      res = Internal.getNodesFromName3(t, Internal.__FlowSolutionCenters__)
      if res == []: return t
      if Internal.getNodesFromType1(res[0], 'DataArray_t') == []: return t
      if ghost is None:
        a = Internal.addGhostCells(t, t, 1,
                                   modified=[Internal.__FlowSolutionCenters__])
      else: a = Internal.copyRef(t)
      fieldc = getFields(Internal.__FlowSolutionCenters__, a)
      fieldn = []
      listVar = []
      for i in fieldc:
        if i != []:
          b = Converter.center2Node(i, cellNType); fieldn.append(b)
          for va in b[0].split(','):
            if va not in listVar: listVar.append(va)
        else: fieldn.append([])
      _patchArrayForCenter2NodeNK1__(fieldn, a)
      setFields(fieldn, a, 'nodes', writeDim=False)
      if ghost is None:
        a = Internal.rmGhostCells(a, a, 1,
                                  modified=[listVar,Internal.__FlowSolutionCenters__])
      return a
    else:
      if isinstance(var, list): vars = var
      else: vars = [var]
      if ghost is None:
        a = Internal.addGhostCells(t, t, 1, modified=vars)
      else: a = Internal.copyRef(t)
      for v in vars: center2Node__(a, v, cellNType)
      # if there are center vars in the list, add equivalent node vars because
      # they have been created by center2Node
      ghost2 = Internal.getNodeFromName(a, 'ZoneRind')
      if (ghost is None and ghost2 is not None):
        var2 = vars[:]
        for v in vars:
          variable = v.split(':')
          if (len(variable) == 2 and variable[0] == 'centers'):
            var2.append(variable[1])
        a = Internal.rmGhostCells(a, a, 1, modified=var2)
      return a

def center2Node__(a, var, cellNType):
  fieldc = getField(var, a)
  fieldn = []
  for i in fieldc:
    if i != []: i = Converter.center2Node(i, cellNType)
    fieldn.append(i)
  _patchArrayForCenter2NodeNK1__(fieldn, a)
  setFields(fieldn, a, 'nodes', writeDim=False)
  return a

# -- node2ExtCenters
# Convert a zone to an extended center zone
# If FlowSolutionCenters exist, they are also converted to extended centers
def node2ExtCenter(t, var=''):
  """Convert zones defined on nodes to zones defined on extended centers.
  Usage: node2ExtCenter(a)"""
  a = Internal.copyRef(t)
  if var == '': # on prend ts les champs
    if Internal.isStdNode(a) == 0: la = a
    else: la = [a] # noeud standard => liste de noeuds standards
    for i in xrange(len(la)):
      zones = Internal.getZones(la[i])
      for z in zones:
        (p, pos) = Internal.getParentOfNode(la[i], z)
        fieldn = getAllFields(z,loc='nodes')
        fieldc = getFields(Internal.__FlowSolutionCenters__, z)
        _deleteFlowSolutions__(z, 'centers')
        _deleteZoneBC__(z)
        _deleteGridConnectivity__(z)
        dim = Internal.getZoneDim(z)[0]
        if dim == 'Structured':
          fieldn2 = Converter.node2ExtCenter(fieldn)
          if fieldc != [[]]:
            fieldc2 = Converter.center2ExtCenter(fieldc)
            fieldn2 = Converter.addVars([fieldn2,fieldc2])
          setFields(fieldn2, z, 'nodes', True)
          if p is None: la[i] = z
          else: p[2][pos] = z
    if Internal.isStdNode(a) == 0: a = la
    else: a = la[0] # liste de noeuds standards => noeud standard
    return a

  elif var == Internal.__GridCoordinates__:
    fieldn = getFields(Internal.__GridCoordinates__, a)
    fielde = Converter.node2ExtCenter(fieldn)
    setFields(fielde, a, 'nodes', writeDim=True)
    return a

  elif var == Internal.__FlowSolutionNodes__:
    fieldn = getFields(Internal.__FlowSolutionNodes__, a)
    fielde = Converter.node2ExtCenter(fieldn)
    setFields(fielde, a, 'nodes', writeDim=True)
    return a
  else:
    raise ValueError("node2ExtCenter: only for all fields, coordinates or flow solution located at nodes.")

#==============================================================================
# diff 2 pyTrees
#==============================================================================
def diffArrays(A, B):
  t1 = Internal.copyRef(A); t2 = Internal.copyRef(B)
  zones1 = Internal.getZones(t1)
  zones2 = Internal.getZones(t2)
  nz = len(zones1)
  if nz != len(zones2):
    raise ValueError("diffArrays: different number of zones (A=%d; B=%d)."%(nz,len(zones2)))
  for no in xrange(nz):
    # noeuds
    A1 = getAllFields(zones1[no], 'nodes'); A1 = Internal.clearList(A1)
    A2 = getAllFields(zones2[no], 'nodes'); A2 = Internal.clearList(A2)
    # elimination des solutions aux noeuds
    node = Internal.getNodesFromName1(zones1[no], Internal.__FlowSolutionNodes__)
    if node != []:
      (parent, d) = Internal.getParentOfNode(t1, node[0])
      if parent is not None: del parent[2][d]

    if (A1 != [] and A2 != []):
      diff = Converter.diffArrays(A1, A2)
      setFields(diff, zones1[no], 'nodes')

    # centres
    A1 = getAllFields(zones1[no], 'centers'); A1 = Internal.clearList(A1)
    A2 = getAllFields(zones2[no], 'centers'); A2 = Internal.clearList(A2)
    node = Internal.getNodesFromName1(zones1[no], Internal.__FlowSolutionCenters__)
    if node != []:
      (parent, d) = Internal.getParentOfNode(t1, node[0])
      if parent is not None: del parent[2][d]

    if (A1 != [] and A2 != []):
      diff = Converter.diffArrays(A1, A2)
      setFields(diff, zones1[no], 'centers')
  t1 = rmNodes(t1, Internal.__GridCoordinates__)
  return t1

#==============================================================================
# - add specific nodes -
#==============================================================================

# -- addState
# Add a single state/value or a full reference state
def addState(t, state=None, value=None, adim='adim1',
             MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8,
             UInf=None, TInf=None, PInf=None, RoInf=None, LInf=None,
             Mus=None, MutSMuInf=0.2, TurbLevelInf=1.e-4):
  """Add single state value or a full reference state."""
  tp = Internal.copyRef(t)
  _addState(tp, state, value, adim,
            MInf, alphaZ, alphaY, ReInf, UInf, TInf, PInf, RoInf, LInf,
            Mus, MutSMuInf, TurbLevelInf)
  return tp

def _addState(t, state=None, value=None, adim='adim1',
              MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8,
              UInf=None, TInf=None, PInf=None, RoInf=None, LInf=None,
              Mus=None, MutSMuInf=0.2, TurbLevelInf=1.e-4):
  ntype = Internal.typeOfNode(t)
  if state is not None and value is not None: # single state value
    addState2Node2__(t, ntype, state, value); return None

  import KCore.Adim
  if state is None: # compute state
    if adim == 'adim1':
      if MInf is None: raise ValueError("addState: MInf is missing.")
      state = KCore.Adim.adim1(MInf, alphaZ, alphaY, ReInf,
                               MutSMuInf, TurbLevelInf)
    elif adim == 'adim2' or adim == 'adim2funk':
      state = KCore.Adim.adim2(MInf, alphaZ, alphaY, ReInf,
                               MutSMuInf, TurbLevelInf)
    elif adim == 'adim3':
      state = KCore.Adim.adim3(MInf, alphaZ, alphaY, ReInf, LInf,
                               MutSMuInf, TurbLevelInf)
    elif adim == 'dim1':
      if UInf is None: raise ValueError("addState: UInf is missing.")
      if TInf is None: raise ValueError("addState: TInf is missing.")
      if PInf is None: raise ValueError("addState: PInf is missing.")
      if LInf is None: raise ValueError("addState: LInf is missing.")
      state = KCore.Adim.dim1(UInf, TInf, PInf, LInf, alphaZ, alphaY,
                              MutSMuInf, TurbLevelInf)
    elif adim == 'dim2':
      if UInf is None: raise ValueError("addState: UInf is missing.")
      if TInf is None: raise ValueError("addState: TInf is missing.")
      if RoInf is None: raise ValueError("addState: RoInf is missing.")
      if LInf is None: raise ValueError("addState: LInf is missing.")
      state = KCore.Adim.dim2(UInf, TInf, RoInf, LInf, alphaZ, alphaY,
                              MutSMuInf, TurbLevelInf)
    elif adim == 'dim3':
      if UInf is None: raise ValueError("addState: UInf is missing.")
      if PInf is None: raise ValueError("addState: PInf is missing.")
      if RoInf is None: raise ValueError("addState: RoInf is missing.")
      if LInf is None: raise ValueError("addState: LInf is missing.")
      state = KCore.Adim.dim3(UInf, PInf, RoInf, LInf, alphaZ, alphaY,
                              MutSMuInf, TurbLevelInf)
    elif adim == 'dim4':
      if UInf is None: raise ValueError("addState: UInf is missing.")
      if TInf is None: raise ValueError("addState: TInf is missing.")
      if PInf is None: raise ValueError("addState: PInf is missing.")
      if LInf is None: raise ValueError("addState: LInf is missing.")
      if Mus is None: raise ValueError("addState: Mus is missing.")
      state = KCore.Adim.dim4(UInf, TInf, PInf, LInf, alphaZ, alphaY,
                              Mus, MutSMuInf, TurbLevelInf)

  UInf   = state[1] / state[0]
  VInf   = state[2] / state[0]
  WInf   = state[3] / state[0]
  addState2Node2__(t, ntype, 'VelocityX', UInf)
  addState2Node2__(t, ntype, 'VelocityY', VInf)
  addState2Node2__(t, ntype, 'VelocityZ', WInf)
  addState2Node2__(t, ntype, 'Density', state[0])
  addState2Node2__(t, ntype, 'MomentumX', state[1])
  addState2Node2__(t, ntype, 'MomentumY', state[2])
  addState2Node2__(t, ntype, 'MomentumZ', state[3])
  addState2Node2__(t, ntype, 'EnergyStagnationDensity', state[4])
  addState2Node2__(t, ntype, 'Pressure', state[5])
  addState2Node2__(t, ntype, 'Temperature', state[6])
  addState2Node2__(t, ntype, 'Cv', state[7])
  addState2Node2__(t, ntype, 'Mach', state[8])
  addState2Node2__(t, ntype, 'Reynolds', state[9])
  addState2Node2__(t, ntype, 'Gamma', state[11])
  addState2Node2__(t, ntype, 'Rok', state[12])
  addState2Node2__(t, ntype, 'RoOmega', state[13])
  addState2Node2__(t, ntype, 'TurbulentSANuTildeDensity', state[14])
  addState2Node2__(t, ntype, 'Mus', state[15])
  addState2Node2__(t, ntype, 'Cs', state[16])
  addState2Node2__(t, ntype, 'Ts', state[17])
  addState2Node2__(t, ntype, 'Pr', state[18])
  return None

# Ajoute un noeud state/value suivant le type de t (in place)
def addState2Node2__(t, ntype, state, value):
  if ntype == 1: # add to zone
    addState2Node__(t, state, value)
  elif ntype == 2: # add to all zones
    for i in t: addState2Node__(i, state, value)
  elif ntype == 3: # add to all bases
    bases = Internal.getBases(t)
    for b in bases: addState2Node__(b, state, value)
  elif ntype == 4: # add to base
    addState2Node__(t, state, value)
  elif ntype == 5: # add to all bases
    for b in tp: addState2Node__(b, state, value)
  else: addState2Node__(t, state, value) # direct to node

# Ajoute un noeud state/value au noeud a (in place)
def addState2Node__(a, state, value):
  # Container: FlowEquationSet or ReferenceState
  if state == 'EquationDimension' or state == 'GoverningEquations' or state == 'TurbulenceModel':
    H = []
    for n in a[2]:
      if n[0] == 'FlowEquationSet': H = n; break

    if H == []:
      a[2].append(['FlowEquationSet', None, [], 'FlowEquationSet_t'])
      H = a[2][len(a[2])-1]
  else:
    H = []
    for n in a[2]:
      if n[0] == 'ReferenceState': H = n; break

    if H == []:
      a[2].append(['ReferenceState', None, [], 'ReferenceState_t'])
      H = a[2][len(a[2])-1]

  # EquationDimension
  if state == 'EquationDimension':
    nodes = Internal.getNodesFromName(H, state)
    if nodes == []:
      v = numpy.empty((1), numpy.int32); v[0] = value
      H[2].append([state, v, [], '"int"']) # Better DataArray_t
    else:
      v = numpy.empty((1), numpy.int32); v[0] = value
      nodes[0][1] = v

  # GoverningEquations
  elif state == 'GoverningEquations':
    nodes = Internal.getNodesFromName(H, state)
    v = numpy.fromstring(value, 'c')
    if nodes == []:
      H[2].append([state, v, [], 'GoverningEquations_t'])
    else:
      nodes[0][1] = v

  # TurbulenceModel
  elif state == 'TurbulenceModel':
    nodes = Internal.getNodesFromName(H, state)
    v = numpy.fromstring(value, 'c')
    if nodes == []:
      H[2].append([state, v, [], 'TurbulenceModel_t'])
    else:
      nodes[0][1] = v

  # Reference state
  else:
    nodes = Internal.getNodesFromName(H, state)
    if nodes == []:
      v = numpy.empty((1), numpy.float64); v[0] = value
      H[2].append([state, v, [], 'DataArray_t'])
    else:
      v = numpy.empty((1), numpy.float64); v[0] = value
      nodes[0][1] = v

  return a

# -- getState
# Retourne un vecteur identique a Adim.XXX
# Doit etre l'exact inverse de addState
def getState__(state, name):
  A = Internal.getNodeFromName1(state, name)
  if A is None: raise ValueError("getState: %s is missing in tree ReferenceState."%name)
  return Internal.getValue(A)

def getState(t):
  state = Internal.getNodeFromName(t, 'ReferenceState')
  if state is None: raise ValueError("getState: ReferenceState is missing in tree.")
  RoInf = getState__(state, 'Density')
  RouInf = getState__(state, 'MomentumX')
  RovInf = getState__(state, 'MomentumY')
  RowInf = getState__(state, 'MomentumZ')
  RoeInf = getState__(state, 'EnergyStagnationDensity')
  PInf = getState__(state, 'Pressure')
  TInf = getState__(state, 'Temperature')
  cvInf = getState__(state, 'Cv')
  MInf = getState__(state, 'Mach')
  ReInf = getState__(state, 'Reynolds')
  Gamma = getState__(state, 'Gamma')
  RokInf = getState__(state, 'Rok')
  RoomegaInf = getState__(state, 'RoOmega')
  RonutildeInf = getState__(state, 'TurbulentSANuTildeDensity')
  Mus = getState__(state, 'Mus')
  Cs = getState__(state, 'Cs')
  Ts = getState__(state, 'Ts')
  Pr = getState__(state, 'Pr')
  return [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
          ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
          Mus, Cs, Ts, Pr]

# -- addChimera2Base
# add a chimera user defined node to a base
# this node stores specific Cassiopee Chimera settings
def addChimera2Base(base, setting, value):
  """Add chimera setting as node in base."""
  basep = Internal.copyRef(base)
  _addChimera2Base(basep, setting, value)
  return basep

def _addChimera2Base(base, setting, value):
  node = Internal.getNodeFromName1(base, '.Solver#Chimera')
  if node is not None:
    chimera = node
  else:
    chimera = ['.Solver#Chimera', None, [], 'UserDefinedData_t']
    base[2].append(chimera)

  # Priority
  if setting == 'Priority':
    v = numpy.empty((1,1), numpy.int32); v[0,0] = value
    node = Internal.getNodeFromName1(chimera, 'Priority')
    if node is not None:
      node[1] = v
    else:
      a = ['Priority', v, [], 'DataArray_t']
      chimera[2].append(a)

  # XRayTol
  elif setting == 'XRayTol':
    v = numpy.empty((1,1), numpy.float64); v[0,0] = value
    node = Internal.getNodeFromName1(chimera, 'XRayTol')
    if node is not None:
      node[1] = v
    else:
      a = ['XRayTol', v, [], 'DataArray_t']
      chimera[2].append(a)

  # XRayDelta
  elif setting == 'XRayDelta':
    v = numpy.empty((1,1), numpy.float64); v[0,0] = value
    node = Internal.getNodeFromName1(chimera, 'XRayDelta')
    if node is not None:
      node[1] = v
    else:
      a = ['XRayDelta', v, [], 'DataArray_t']
      chimera[2].append(a)

  # DoubleWallTol
  elif setting == 'DoubleWallTol':
    v = numpy.empty((1,1), numpy.float64); v[0,0] = value
    node = Internal.getNodeFromName1(chimera, 'DoubleWallTol')
    if node is not None:
      node[1] = v
    else:
      a = ['DoubleWallTol', v, [], 'DataArray_t']
      chimera[2].append(a)

  # relations
  elif (setting == '+' or setting == '-' or setting == '0' or setting == 'N'):
    v = numpy.fromstring(value, 'c')
    a = [setting, v, [], 'UserDefinedData_t']
    chimera[2].append(a)

  return None

#==============================================================================
# Retourne un arbre ou certaines bases ont ete eliminees
#==============================================================================
def getPrimaryTree(t):
  import re
  tp = newPyTree()
  excluded = ['CONTOURS', 'BODY#', 'SLICES', 'SURFACES']
  bases = t[1:]
  for b in bases:
    found = False
    if b[3] == 'CGNSBase_t':
      for i in excluded:
        exp = re.compile(i)
        if exp.search(b[0]) is not None:
          found = True
          break
    if not found: tp.append(b)
  return tp

#==============================================================================
# check if name is already in list
# return True or False
#==============================================================================
def checkNameInList(name, list):
  for i in list:
    if name == i: return True
  return False

#==============================================================================
# mergeTrees
# Merge 2 arbres. Les bases de t2 sont ajoutees a t1. Les noms des bases
# sont modifies si elles existent deja dans t1.
#==============================================================================
def mergeTrees(t1, t2):
  t1p = Internal.copyRef(t1)
  # Enregistre les noms des bases de t1
  t1BaseNames = []
  bases = Internal.getBases(t1)
  for b in bases: t1BaseNames.append(b[0])

  bases = Internal.getBases(t2)
  for b in bases:
    ret = checkNameInList(b[0], t1BaseNames)
    if not ret: t1p[2].append(b)
    else:
      c = 0
      while ret == True:
        ret = checkNameInList('%s.%d'%(b[0],c), t1BaseNames)
        c += 1
      b[0] = '%s.%d'%(b[0],c-1)
      t1p[2].append(b)
  return t1p

#==============================================================================
# Retourne les faces des BCs pour les zones non structurees et pour
# les BCs definies par faces (BC physiques) une liste par zone
# [ ['nomBC', numpy(faces), 'nomBC', numpy(faces)...] ...]
# nomBC concatene le nom de la BC et son type (BCName#BCType)
#==============================================================================
def getBCFaces(t):
  BCFaces = []
  zones = Internal.getZones(t)
  for z in zones:
    BCFZ = []
    BCs = Internal.getNodesFromType2(z, 'BC_t')
    for b in BCs:
      name = b[0]
      BCtype = Internal.getValue(b)
      n = Internal.getNodesFromName1(b, 'PointList')
      p = Internal.getNodesFromType1(b, 'GridLocation_t')
      if p != []:
        loc = Internal.getValue(p[0])
        if n != [] and loc == 'FaceCenter':
          if BCtype == 'FamilySpecified':
            familyName = Internal.getNodeFromType1(b, 'FamilyName_t')
            if familyName is not None:
              ft = Internal.getValue(familyName)
              BCFZ += [name+Internal.SEP2+ft, n[0][1]]
            else: BCFZ += [name, n[0][1]]
          else: BCFZ += [name+Internal.SEP2+BCtype, n[0][1]]
    BCFaces.append(BCFZ)
  return BCFaces

#==============================================================================
# Ajoute les BCFaces dans l'arbre
#==============================================================================
def addBCFaces(t, BCFaces):
  tp = Internal.copyRef(t)
  _addBCFaces(t, BCFaces)
  return tp

def _addBCFaces(t, BCFaces):
  if BCFaces == []: return None
  zones = Internal.getZones(t)
  nz = 0; c = 0
  for z in zones:
    myBC = BCFaces[c]
    l = len(myBC)
    for i in xrange(l/2):
      name = myBC[2*i]
      names = name.split(Internal.SEP2)
      if len(names) == 2:
        name = names[0]; bctype = names[1]
        if bctype not in Internal.KNOWNBCS:
         name = myBC[2*i]; bctype = 'FamilySpecified:'+name
      else: name = names[0]; bctype = 'FamilySpecified:'+names[0]
      faces = myBC[2*i+1]
      _addBC2Zone(z, name, bctype, faceList=faces)
    c += 1
  return None

#==============================================================================
# -- Fonctions de preconditionnement (hook) --
#==============================================================================

# -- createHook
def createHook(a, function='None'):
  """Create a hook for a given function.
    Usage: hook = createHook(a, function)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if (function == 'extractMesh' or function == 'adt'):
    return Converter.createHook(fields, function)
  else:
    if len(fields) == 1: return Converter.createHook(fields[0], function)
    else: return Converter.createHook(fields, function)

# -- createGlobalHook
def createGlobalHook(a, function='None',extended=0):
  """Create a global hook for all zones in a.
  Usage: hook = createGlobalHook(a, function)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  return Converter.createGlobalHook(fields, function,extended)

# -- freeHook
def freeHook(hook):
  """Free hook.
  Usage: freeHook(hook)"""
  Converter.freeHook(hook)

#==============================================================================
# -- Fonctions d'identification geometrique --
#==============================================================================

# -- identifyNodes: identifie les noeuds de a dans hook
def identifyNodes(hook, a, tol=1.e-11):
  """Find in a hook nearest points of nodes of a. return identified node indices.
  Usage: identifyNodes(hook, a)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if len(fields) == 1: return Converter.identifyNodes(hook, fields[0], tol)
  else: return Converter.identifyNodes(hook, fields, tol)

# -- identifyFaces: identifie les centres de faces de a dans hook
def identifyFaces(hook, a, tol=1.e-11):
  """Find in a hook nearest points of face centers of a. return identified face indices.
  Usage: identifyFaces(hook, a)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if len(fields) == 1: return Converter.identifyFaces(hook, fields[0], tol)
  else: return Converter.identifyFaces(hook, fields, tol)

# -- identifyElements: identifie le centre des elements de a dans hook
def identifyElements(hook, a, tol=1.e-11):
  """Find in a hook nearest points of element centers of a. return identified element indices.
  Usage: identifyElements(hook, a)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if len(fields) == 1: return Converter.identifyElements(hook, fields[0], tol)
  else: return Converter.identifyElements(hook, fields, tol)

# -- identifySolutions: recopie la solution de tDnr dans tRcv par identification
# des noeuds et des centres
def identifySolutions(tRcv, tDnr, hookN=None, hookC=None, vars=[], tol=1.e-12):
  """Identify points in stored in a global hook to mesh points and set the solution if donor
  and receptor points are distant from tol.
  Usage: identifySolutions(tRcv, tDnr, hookN, hookC, vars, tol)"""
  varsC=[]; varsN=[]
  for v in vars:
    s = v.find('centers:')
    if s != -1: varsC.append(v)
    else: varsN.append(v)

  if varsC == [] and hookC is not None:
    varsC = getVarNames(tDnr, excludeXYZ=True, loc='centers')
    if varsC != []: varsC = varsC[0]
  if varsN == [] and hookN is not None:
    varsN = getVarNames(tDnr, excludeXYZ=True, loc='nodes')
    if varsN != []: varsN = varsN[0]

  if len(varsC) == 0 and len(varsN) == 0: return tRcv
  if varsC != []:
    for nov in xrange(len(varsC)):
      vc = varsC[nov].split(':')[1]
      varsC[nov] = vc

  t2 = Internal.copyRef(tRcv)
  t2p, typen = Internal.node2PyTree(t2)

  zones = Internal.getZones(t2p)
  coordsR = getFields(Internal.__GridCoordinates__, zones)

  if varsN != [] and hookN is not None:
    fnodes = getFields(Internal.__FlowSolutionNodes__, tDnr)
    resn = Converter.identifySolutions(coordsR,fnodes,hookN,vars=varsN,tol=tol)
    setFields(resn, zones, 'nodes')
  if varsC != [] and hookC is not None:
    fcenters = getFields(Internal.__FlowSolutionCenters__,tDnr)
    centersR = Converter.node2Center(coordsR)
    resc=Converter.identifySolutions(centersR,fcenters,hookC,vars=varsC,tol=tol)
    setFields(resc,zones,'centers')

  t2 = Internal.pyTree2Node(t2p,typen)
  return t2

# -- nearestNodes: identifie le noeud de a le plus proche d'un point de hook
def nearestNodes(hook, a):
  """Identify nearest nodes to a in hook. return identified face indices.
  Usage: nearestNodes(hook, a)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if len(fields) == 1: return Converter.nearestNodes(hook, fields[0])
  else: return Converter.nearestNodes(hook, fields)

# -- nearestFaces: identifie la face de a la plus proches d'un point de hook
def nearestFaces(hook, a):
  """Identify nearest face centers to a in hook. return identified face indices.
  Usage: nearestFaces(hook, a)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if len(fields) == 1: return Converter.nearestFaces(hook, fields[0])
  else: return Converter.nearestFaces(hook, fields)

# -- nearestElements: identifie le centre de l'elements de a le plus proche
# d'un point de hook
def nearestElements(hook, a):
  """Identify nearest element centers to a in hook. return identified face indices.
  Usage: nearestFaces(hook, a)"""
  fields = getFields(Internal.__GridCoordinates__, a)
  if len(fields) == 1: return Converter.nearestElements(hook, fields[0])
  else: return Converter.nearestElements(hook, fields)

#==============================================================================
# -- Connectivity management --
#==============================================================================

# -- mergeConnectivity
# IN: z1: zone BE
# IN: z2: zone BE (to be merged in z1) avec un seul type d'elements
# si boundary==1, les noeuds de z2 sont identifies dans z1
# si boundary==0, les noeuds de z2 sont merges dans z1 et reidentifies
# IN: boundary: 0 (not a boundary zone), 1 (a boundary zone, add it as a
# boundary connectivity)
def mergeConnectivity(z1, z2, boundary=0):
  """Gather an additional zone connectivity in z1.
  Usage: mergeConnectivity(z1, z2)"""
  zout = Internal.copyRef(z1)
  _mergeConnectivity(zout, z2, boundary)
  return zout

def _mergeConnectivity(z1, z2, boundary=0):
  # Analyse zone z2
  dims = Internal.getZoneDim(z2)
  neb = dims[2] # nbre d'elts de z2
  eltType, nf = Internal.eltName2EltNo(dims[3]) # type d'elements de z2

  # On cherche l'element max dans les connectivites de z1
  maxElt = 0
  connects = Internal.getNodesFromType(z1, 'Elements_t')
  for cn in connects:
    r = Internal.getNodeFromName1(cn, 'ElementRange')
    m = r[1][1]
    maxElt = max(maxElt, m)

  # Si boundary=0, on fusionne les coordonnees
  if boundary == 0:
    import Transform.PyTree as T
    zn1 = convertArray2Node(z1)
    zn2 = convertArray2Node(z2)
    zn = T.join(zn1, zn2)
    # reset Coordinates in z1
    _deleteFlowSolutions__(z1)
    cont1 = Internal.getNodeFromName(z1, Internal.__GridCoordinates__)
    contn = Internal.getNodeFromName(zn, Internal.__GridCoordinates__)
    for name in ['CoordinateX', 'CoordinateY', 'CoordinateZ']:
      p1 = Internal.getNodeFromName(cont1, name)
      pn = Internal.getNodeFromName(contn, name)
      p1[1] = pn[1]
    np = Internal.getZoneDim(zn)[1]
    z1[1] = numpy.copy(z1[1])
    z1[1][0,0] = np # nouveau nbre de pts

    # Reidentifie la connectivite de z1
    hook = createHook(zn, 'nodes')
    ids = identifyNodes(hook, z1)
    nodes = Internal.getNodesFromType(z1, 'Elements_t')
    for n in nodes:
      node = Internal.getNodeFromName(n, 'ElementConnectivity')
      oldc = node[1]
      newc = numpy.copy(oldc)
      newc[:] = ids[oldc[:]-1]
      node[1] = newc

    # Ajoute les connectivites de z2
    nebb = 0
    z1[1][0,1] += neb # nouveau nbre de cellules

    # on identifie les noeuds de z2 dans zn
    ids = identifyNodes(hook, z2)

    # on cree un nouveau noeud connectivite dans z1 (avec le nom de la zone z2)
    # elts = Internal.getNodesFromType2(z2, 'Elements_t'); c = 0
    # for e in elts:
    #   Internal.createUniqueChild(z1, z2[0]+'.'+str(c), 'Elements_t', value=[eltType,nebb])
    #   node = Internal.getNodeFromName(z1, z2[0])
    #   Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
    #                               value=[maxElt+1,maxElt+neb])
    #   oldc = Internal.getNodeFromName(e, 'ElementConnectivity')[1]
    #   newc = numpy.copy(oldc)
    #   newc[:] = ids[oldc[:]-1]
    #   Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)
    #   c += 1
    Internal.createUniqueChild(z1, z2[0], 'Elements_t', value=[eltType,nebb])
    node = Internal.getNodeFromName(z1, z2[0])
    Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                               value=[maxElt+1,maxElt+neb])
    oldc = Internal.getNodeFromName(z2, 'ElementConnectivity')[1] # first
    newc = numpy.copy(oldc)
    newc[:] = ids[oldc[:]-1]
    Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)

  else: # connectivite boundary (subzone)
    # on identifie les noeuds de z2 dans z1
    hook = createHook(z1, 'nodes')
    ids = identifyNodes(hook, z2)
    z1[1] = numpy.copy(z1[1])
    z1[1][0,2] += neb

    # on cree un nouveau noeud connectivite dans z1 (avec le nom de la zone z2)
    nebb = neb
    Internal.createUniqueChild(z1, z2[0], 'Elements_t', value=[eltType,nebb])
    node = Internal.getNodeFromName(z1, z2[0])
    Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                                value=[maxElt+1,maxElt+neb])
    oldc = Internal.getNodeFromName(z2, 'ElementConnectivity')[1] # first
    newc = numpy.copy(oldc)
    newc[:] = ids[oldc[:]-1]
    Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t',
                               value=newc)
  return None

# -- breakConnectivity
# break a multiple connectivity zone into single elements zones
# don't break boundary connectivity
# IN: t: t to break
def breakConnectivity(t):
    """Break a multi-element zone into single element zones.
    Usage: breakConnectivity(t)"""
    tp, typen = Internal.node2PyTree(t)
    bases = Internal.getBases(tp)
    for b in bases:
        c = 0; l = len(b[2])
        for c in xrange(l):
            z = b[2][c]
            if z[3] == 'Zone_t':
                # compte les connectivites elements (hors boundary)
                connects = Internal.getElementNodes(z)
                N = len(connects)
                if (N <= 1): break # une seule connectivite
                if (N == 2):
                  type1 = connects[0][1][0]; type2 = connects[1][1][0]
                  if ((type1 == 22 and type2 == 23) or (type1 == 23 and type2 == 22)): # pur NGON
                    break;

                iBE = []; iNGon = -1; iNFace = -1; i = 0
                for co in connects:
                  ctype = co[1][0]
                  if (ctype == 22): iNGon = i
                  elif (ctype == 23): iNFace = i
                  else: iBE.append(i)
                  i += 1

                N = len(iBE) # split les BE
                for p in xrange(N):
                  i = iBE[p]
                  zp = Internal.copyRef(z)
                  _deleteFlowSolutions__(zp, 'centers')
                  _deleteGridConnectivity__(zp)
                  _deleteZoneBC__(zp)
                  GEl = Internal.getNodesFromType1(zp, 'Elements_t')
                  GE = Internal.getNodeFromName1(zp, connects[i][0])
                  eltType, nf = Internal.eltNo2EltName(GE[1][0])
                  zp[0] = getZoneName(z[0]+'_'+eltType)
                  # Enleve toutes les connects a part la ieme
                  for GEj in GEl:
                    if (GEj is not GE): Internal._rmNodesByName(zp, GEj[0])
                  # Renumerote la connectivite
                  r = Internal.getNodeFromName(GE, 'ElementRange'); r = r[1]
                  start = r[0]; end = r[1]
                  # modifie le range de la connectivite
                  r[0] = 1; r[1] = end-start+1
                  # modifie le nbre d'elements dans la zone
                  zp[1] = numpy.copy(z[1])
                  zp[1][0,1] = end-start+1
                  if (i == 0): b[2][c] = zp
                  else: b[2] += [zp]
                  #zp = pushBC(z, zp, type='BCC')

                if (iNGon != -1 and iNFace != -1): # NGon additionnel
                  i1 = iNGon; i2 = iNFace
                  zp = Internal.copyRef(z)
                  _deleteFlowSolutions__(zp, 'centers')
                  _deleteGridConnectivity__(zp)
                  _deleteZoneBC__(zp)
                  GEl = Internal.getNodesFromType1(zp, 'Elements_t')
                  GE1 = Internal.getNodeFromName1(zp, connects[i1][0])
                  GE2 = Internal.getNodeFromName1(zp, connects[i2][0])
                  zp[0] = getZoneName(z[0]+'_NGON')
                  # Enleve toutes les connects a part la ieme
                  for GEj in GEl:
                    if (GEj is not GE1 and GEj is not GE2): Internal._rmNodesByName(zp, GEj[0])
                  # Renumerote la connectivite
                  r = Internal.getNodeFromName(GE1, 'ElementRange'); r = r[1]
                  start = r[0]; end = r[1]
                  # modifie le range de la connectivite
                  r[0] = 1; r[1] = end-start+1; fin = end-start+1
                  r = Internal.getNodeFromName(GE2, 'ElementRange'); r = r[1]
                  start = r[0]; end = r[1]
                  # modifie le range de la connectivite
                  r[0] = fin; r[1] = fin+end-start+1
                  # modifie le nbre d'elements dans la zone
                  zp[1] = numpy.copy(z[1])
                  zp[1][0,1] = end-start+1
                  b[2] += [zp]
                  #zp = pushBC(z, zp, type='F')

    if typen == 1: typen = 2 # une zone renvoie une liste de zones
    return Internal.pyTree2Node(tp, typen)

#==============================================================================
# renumerote les ElementRanges pour qu'ils soient dans un ordre contigue
# y compris pour les BCCs
#==============================================================================
def _renumberElementConnectivity(t):
  zones = Internal.getZones(t)
  for z in zones:
    elts = Internal.getNodesFromType1(z, 'Elements_t')
    c = 1
    for e in elts:
      r = Internal.getNodeFromName1(e, 'ElementRange')
      r[1] = numpy.copy(r[1]); r = r[1]
      delta = r[1]-r[0]+1
      r[0] = c; r[1] = c+delta-1; c += delta
  return None

# -- selectOneConnectivity
# Retourne une nouvelle zone avec une seule de ses connectivites
# (passee en connectivite volumique)
# IN: name: nom de la ElementConnectivity a conserver
# or IN: number: no de la ElementConnectivity a conserver (first=0)
# or IN: range=[min,max]
def selectOneConnectivity(z, name=None, number=None, range=None):
  zp = Internal.copyRef(z)
  elts = Internal.getNodesFromType1(zp, 'Elements_t')

  if name is not None:
    for e in elts:
      if e[0] != name: Internal._rmNodesByName(zp, e[0])
      else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
  elif number is not None:
    c = 0
    for e in elts:
      if (c != number): Internal._rmNodesByName(zp, e[0])
      else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
      c += 1
  elif range is not None:
    for e in elts:
      r = Internal.getNodeFromName1(e, 'ElementRange')
      if r is not None:
        r = r[1]
        if (r[0] > range[0] or r[1] < range[1]):
          Internal._rmNodesByName(zp, e[0])
        else:
          if (r[0] == range[0] and r[1] == range[1]): # full
            e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
          else: # slice
            (name, nnodes) = Internal.eltNo2EltName(e[1][0])
            if (name != 'NGON' and name != 'NFACE' and name != 'MIXED'):
              e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
              r = Internal.getNodeFromName1(e, 'ElementRange')
              r[1] = numpy.copy(r[1]); r[1][0] = 1; r[1][1] = range[1]-range[0]+1
              c = Internal.getNodeFromName1(e, 'ElementConnectivity')
              c[1] = c[1][nnodes*(range[0]-1):nnodes*(range[1])+1]
            else: print 'Warning: selectOneConnectivity: slice impossible.'
      else: Internal._rmNodesByName(zp, e[0])
  _renumberElementConnectivity(zp)
  return zp

# -- selectConnectivity
# Retourne une nouvelle zone avec une seule de ses connectivites
# (passee en connectivite volumique)
# IN: name: nom de la ElementConnectivity a conserver
# or IN: number: no de la ElementConnectivity a conserver (first=0)
# or IN: range=[min,max]
def selectConnectivity(z, name=None, number=None, range=None):
  zp = Internal.copyRef(z)
  elts = Internal.getNodesFromType1(zp, 'Elements_t')

  if name is not None:
    for e in elts:
      if e[0] != name: Internal._rmNodesByName(zp, e[0])
      else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
  elif number is not None:
    c = 0
    for e in elts:
      if (c != number): Internal._rmNodesByName(zp, e[0])
      else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
      c += 1
  elif range is not None:
    for e in elts:
      r = Internal.getNodeFromName1(e, 'ElementRange')
      if r is not None:
        r = r[1]
        if (r[0] != range[0] and r[1] != range[1]):
          Internal._rmNodesByName(zp, e[0])
        if (r[0] == range[0] and r[1] == range[1]): # full
          # print 'no slice'
          e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
        if (r[0] == range[0] and r[1] < range[1]): # full
          # print 'slice1',r[0],r[1]
          (name, nnodes) = Internal.eltNo2EltName(e[1][0])
          if (name != 'NGON' and name != 'NFACE' and name != 'MIXED'):
            e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
            r0 = r[0]; r1 = r[1]
            r = Internal.getNodeFromName1(e, 'ElementRange')
            r[1] = numpy.copy(r[1]); r[1][0] = 1; r[1][1] = r1-r0+1
            c = Internal.getNodeFromName1(e, 'ElementConnectivity')
            r = r[1]
            c[1] = c[1][nnodes*(r[0]-1):nnodes*(r[1])+1]
          else: print 'Warning: selectConnectivity: slice impossible.'
        if (r[0] > range[0] and r[1] == range[1]): # full
          # print 'slice2',r[0],r[1]
          (name, nnodes) = Internal.eltNo2EltName(e[1][0])
          if (name != 'NGON' and name != 'NFACE' and name != 'MIXED'):
            e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
            r0 = r[0]; r1 = r[1]
            r = Internal.getNodeFromName1(e, 'ElementRange')
            r[1] = numpy.copy(r[1]); r[1][0] = 1; r[1][1] = r1-r0+1
            c = Internal.getNodeFromName1(e, 'ElementConnectivity')
            r = r[1]
            c[1] = c[1][nnodes*(r[0]-1):nnodes*(r[1])+1]
          else: print 'Warning: selectConnectivity: slice impossible.'
      else: Internal._rmNodesByName(zp, e[0])
  _renumberElementConnectivity(zp)
  return zp


#=============================================================================
# Rm duplicated periodic zones in a according to the node 'TempPeriodicZone'
#=============================================================================
def removeDuplicatedPeriodicZones__(a):
    for z in Internal.getZones(a):
        parent,d = Internal.getParentOfNode(a, z)
        isperiod = Internal.getNodeFromName1(z, 'TempPeriodicZone')
        if isperiod is not None: del parent[2][d]
    return a

#=============================================================================
# Extract periodic zone info and duplicate zone in same base
#=============================================================================
def duplicatePeriodicZones__(a):
    try: import Transform.PyTree as T
    except:
        raise ImportError("duplicatePeriodicZones__: requires Transform module.")
    atype = Internal.typeOfNode(a)
    if atype != 4:  # base
        print 'Warning: duplicatePeriodicZones__: input node must be a CGNS basis.'
        print 'Skipped.'
        return a

    zones = Internal.getNodesFromType1(a, 'Zone_t')
    zonesdup = []
    for z in zones:
      zname = z[0]
      # Periodic match info
      gcmatch = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
      for gc in gcmatch:
        rotationData, translVect = Internal.getPeriodicInfo__(gc)

        if translVect == []: translated = False
        else: translated = True

        if rotationData == []: rotated = False
        else: rotated = True
        if rotated == True and translated == True:
          print 'Warning: duplicatePeriodicZones__: rotation and translation cannot be applied at the same time. %s periodic grid connectivity not taken into account.'%gc[0]
        elif rotated == False and translated == False: pass
        else:
          zdonorname = Internal.getValue(gc)
          zd = Internal.getNodeFromName(a,zdonorname)
          if zd != []:
            if rotated:
              xc = rotationData[0]; yc = rotationData[1]; zc = rotationData[2]
              vx = rotationData[3]; vy = rotationData[4]; vz = rotationData[5]
              angle = rotationData[6]
              if angle != 0.:
                zddup = T.rotate(zd,(xc,yc,zc), (vx,vy,vz),-angle)
                if angle >0.: signangle = -1
                else: signangle = 1
                rotationInfo = Internal.createNode('SignOfRotationAngle','UserDefinedData_t',value=signangle)
                # creation du noeud temporaire le marquant comme periodique
                zddup[0] = getZoneName(zddup[0]+'_dup')
                Internal.createChild(zddup,'TempPeriodicZone','UserDefinedData_t',value=zdonorname,children=[rotationInfo])
                zonesdup.append(zddup)
            elif translated:
              tvx=translVect[0];tvy=translVect[1];tvz=translVect[2]
              zddup = T.translate(zd,(tvx,tvy,tvz))
              # creation du noeud temporaire le marquant comme periodique
              zddup[0] = getZoneName(zddup[0]+'_dup')
              Internal.createChild(zddup,'TempPeriodicZone','UserDefinedData_t',value=zdonorname,children=None)
              zonesdup.append(zddup)

      # Chimere periodique: compatible avec elsA uniquement
      usd = Internal.getNodeFromName2(z,".Solver#Param")
      if usd is not None:
        perdir = Internal.getNodeFromName1(usd,'periodic_dir')
        if perdir is not None:
          perdir = Internal.getValue(perdir)
          ang1 = Internal.getNodeFromName1(usd,'axis_ang_1'); ang1 = Internal.getValue(ang1)
          ang2 = Internal.getNodeFromName1(usd,'axis_ang_2'); ang2 = Internal.getValue(ang2)
          angle = float(ang2)/float(ang1)*360.
          xc = Internal.getNodeFromName1(usd,'axis_pnt_x'); xc = Internal.getValue(xc)
          yc = Internal.getNodeFromName1(usd,'axis_pnt_y'); yc = Internal.getValue(yc)
          zc = Internal.getNodeFromName1(usd,'axis_pnt_z'); zc = Internal.getValue(zc)
          vx = Internal.getNodeFromName1(usd,'axis_vct_x'); vx = Internal.getValue(vx)
          vy = Internal.getNodeFromName1(usd,'axis_vct_y'); vy = Internal.getValue(vy)
          vz = Internal.getNodeFromName1(usd,'axis_vct_z'); vz = Internal.getValue(vz)
          if angle != 0.:
            if perdir == 1 or perdir == 3:
              zdup = T.rotate(z,(xc,yc,zc), (vx,vy,vz),angle)
              signangle = 1
              rotationInfo = Internal.createNode('SignOfRotationAngle','UserDefinedData_t',value=signangle)
              # creation du noeud temporaire le marquant comme periodique
              zdup[0] = getZoneName(zdup[0]+'_dup')
              Internal.createChild(zdup,'TempPeriodicZone','UserDefinedData_t',value=zname,children=[rotationInfo])
              zonesdup.append(zdup)
            if perdir == 2 or perdir == 3:
              zdup = T.rotate(z,(xc,yc,zc), (vx,vy,vz),-angle)
              signangle =-1
              rotationInfo = Internal.createNode('SignOfRotationAngle','UserDefinedData_t',value=signangle)
              # creation du noeud temporaire le marquant comme periodique
              zdup[0] = getZoneName(zdup[0]+'_dup')
              Internal.createChild(zdup,'TempPeriodicZone','UserDefinedData_t',value=zname,children=[rotationInfo])
              zonesdup.append(zdup)
    a[2] += zonesdup
    return a
        
#==============================================================================
def convertPyTree2FFD(zone, RefStat, FlowEq, nd):
  print "nd=", nd
  print "RefStat=", RefStat
  print "FlowEq =", FlowEq
  Converter.converter.convertPyTree2FFD(zone,RefStat,FlowEq,nd,
                                        Internal.__GridCoordinates__,
                                        Internal.__FlowSolutionNodes__,
                                        Internal.__FlowSolutionCenters__ )
  return None

#==============================================================================
# - client/server -
#==============================================================================
def send(data, host, rank=0, port=15555):
    """Send data to server."""
    Converter.send(data, host, rank, port)

def createSockets(nprocs=1, port=15555):
    """Create sockets."""
    return Converter.createSockets(nprocs, port)

def listen(s):
    """Listen to one socket."""
    return Converter.listen(s)

