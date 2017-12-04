"""Find grid connectivity.
"""
#
# Python Interface to define grid connectivity
#
import Connector
import connector
import numpy
__version__ = Connector.__version__
try:
    import KCore
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter.cgnslib as CGL
    import Converter.cgnskeywords as CGK
except: raise ImportError("Connector.PyTree: requires Converter module.")

# Variables IBM pour le post traitement
__PRESSURE__='Pressure'
__UTAU__    ='utau'
__YPLUS__   ='yplus'
__DENSITY__ ='Density'

#==============================================================================
# connectMatch between NGON zones
#==============================================================================
def _connectMatchNGON__(a, tol, dim, glob, allExtFaces=None, allExtIndices=None, periodic=0,
                        rotationCenter=None,rotationAngle=None,Translation=None,signT=0,signR=0):
    try: import Post.PyTree as P; import Transform.PyTree as T
    except: raise ImportError("connectMatchNGON__ requires Post and Transform modules.")

    # -------------------------
    # Exterior faces + indices
    # -------------------------
    zones = Internal.getZones(a)
    zonesp = []; indirZones = []
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[3] == 'NGON':
            zonesp.append(z)
            indirZones.append(noz)
        noz += 1
    if len(zonesp)==0: return glob

    # Get exterior faces
    if allExtFaces is None or allExtIndices is None:
        allExtIndices = []
        allExtFaces = P.exteriorFaces(zonesp, indices=allExtIndices)
        C._initVars(allExtFaces,'{centers:tag1}=-1.')# defines the opposite window
        C._initVars(allExtFaces,'{centers:tag2}=-1.') # defines the opp index in opp window

    # identify matching exterior faces
    if allExtFaces != []:
        nzones = len(zonesp)
        tagsF = C.node2Center(allExtFaces)
        tagsF = C.getAllFields(tagsF, 'nodes')
        tagsF = Connector.identifyMatching(tagsF, tol)
        infos = Connector.gatherMatchingNGon__(tagsF, allExtIndices)
        rcvZones = infos[0]
        dnrZones = infos[1]
        allListRcvFaces = infos[2]
        allListDnrFaces = infos[3]
        nmatch = rcvZones.shape[0]
        for nm in xrange(nmatch):
            noz1p = rcvZones[nm]; noz2p = dnrZones[nm]
            isok = 1
            if periodic == 1:
                if noz1p >= nzones and noz2p < nzones: noz1p = noz1p-nzones
                elif noz2p >= nzones and noz1p < nzones: noz2p = noz2p-nzones
                else: isok = 0

            if isok == 1:
                noz1 = indirZones[noz1p]; noz2 = indirZones[noz2p]
                z1OppName = zones[noz1][0]
                z2OppName = zones[noz2][0]
                name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
                name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
                faceListR = allListRcvFaces[nm]
                faceListD = allListDnrFaces[nm]
                if periodic == 1:
                    tsign = numpy.array([signT,signT,signT])
                    rsign = numpy.array([signR,signR,signR])
                    C._addBC2Zone(zones[noz1],name1,'BCMatch',faceList=faceListR,\
                                      zoneDonor=z2OppName, faceListDonor=faceListD,tol=tol,\
                                      rotationCenter=rotationCenter, rotationAngle=rsign*rotationAngle,\
                                      translation=tsign*Translation)
                    tsign = numpy.array([-signT,-signT,-signT])
                    rsign = numpy.array([-signR,-signR,-signR])
                    C._addBC2Zone(zones[noz2],name2,'BCMatch',faceList=faceListD,\
                                      zoneDonor=z1OppName, faceListDonor=faceListR,tol=tol,\
                                      rotationCenter=rotationCenter, rotationAngle=rsign*rotationAngle,\
                                      translation=tsign*Translation)
                else:
                    C._addBC2Zone(zones[noz1],name1,'BCMatch',faceList=faceListR,\
                                      zoneDonor=z2OppName, faceListDonor=faceListD,tol=tol)
                    C._addBC2Zone(zones[noz2],name2,'BCMatch',faceList=faceListD,\
                                      zoneDonor=z1OppName, faceListDonor=faceListR,tol=tol)

        C.setFields(tagsF, allExtFaces, 'centers')
    # Sortie
    return glob
#==============================================================================
# connectMatch between NGON zones
#==============================================================================
def _connectMatchHybrid__(a, tol, dim, glob):
    try: import Post.PyTree as P; import Transform.PyTree as T
    except: raise ImportError("connectMatchHybrid__ requires Post and Transform modules.")

    # -------------------------
    # Exterior faces + indices
    # -------------------------
    zones = []; indirZones = []
    noz = 0

    # identifyMatching runs over structured then unstructure. Same order must be kept
    # that s why 2 lists are first built and merged
    nzonesS = 0; nzonesU = 0
    for z in Internal.getZones(a):
        dimZ = Internal.getZoneDim(z)
        if dimZ[0]=='Structured':
            zones.append(z)
            indirZones.append(noz)
            noz += 1
            nzonesS+=1
    for z in Internal.getZones(a):
        dimZ = Internal.getZoneDim(z)
        if dimZ[3] == 'NGON':
            zones.append(z)
            indirZones.append(noz)
            noz+=1
            nzonesU+=1
    nzones = len(zones)
    if len(zones)==0: return glob
    if nzonesS == 0 or nzonesU==0: return glob
    # Get exterior faces
    allExtIndices=[]; allExtFaces=[]
    for z in zones:        
        extIndices=[]
        extFaces = P.exteriorFaces(z, indices=extIndices)        
        Internal._rmNodesByType(extFaces,'ZoneBC_t')
        Internal._rmNodesByType(extFaces,'ZoneGridConnectivity_t')
        if Internal.getZoneType(z)==1:
            extFaces = C.convertArray2NGon(extFaces)
            C._initVars(extFaces,'{centers:tag1}=-2.')
        else: C._initVars(extFaces,'{centers:tag1}=-1.')
        C._initVars(extFaces,'{centers:tag2}=-1.') # defines the opp index in opp window
        allExtFaces.append(extFaces); allExtIndices+=extIndices

    # identify matching exterior faces
    if allExtFaces != []:
        tagsF = C.node2Center(allExtFaces)
        tagsF = C.getAllFields(tagsF, 'nodes')
        tagsF = Connector.identifyMatching(tagsF, tol)
        infos = Connector.gatherMatchingNGon__(tagsF, allExtIndices)
        rcvZones = infos[0]
        dnrZones = infos[1]
        allListRcvFaces = infos[2]
        allListDnrFaces = infos[3]
        nmatch = rcvZones.shape[0]
        for nm in xrange(nmatch):
            noz1p = rcvZones[nm]; noz2p = dnrZones[nm]
            noz1 = indirZones[noz1p]; noz2 = indirZones[noz2p]
            z1OppName = zones[noz1][0]
            z2OppName = zones[noz2][0]
            name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
            name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
            faceListR = allListRcvFaces[nm]
            faceListD = allListDnrFaces[nm]
            C._addBC2Zone(zones[noz1],name1,'BCMatch',faceList=faceListR,\
                          zoneDonor=z2OppName, faceListDonor=faceListD,tol=tol)
            C._addBC2Zone(zones[noz2],name2,'BCMatch',faceList=faceListD,\
                          zoneDonor=z1OppName, faceListDonor=faceListR,tol=tol)

        C.setFields(tagsF, allExtFaces, 'centers')
    # Sortie
    return glob
#==============================================================================
# Computes the 1 to 1 grid connectivity between structured zones
#==============================================================================
def _connectMatchStruct__(a, tol, dim, glob):
    zones = []
    for z in Internal.getZones(a):
        if Internal.getZoneType(z)==1: zones.append(z)
    # extract empty windows
    structTags,structWins,structIndirBlkOfWins,typeOfWins,dimsI,dimsJ,dimsK = \
        getEmptyWindowsInfoStruct__(zones,dim)

    model ='Euler'
    bases = Internal.getBases(a)
    if bases != []:
        c = Internal.getNodeFromName2(bases[0], 'GoverningEquations')
        if c is not None: model = Internal.getValue(c)

    # Identify matching cells for structured zones
    if structTags != []:
        structTags = Connector.identifyMatching(structTags,tol)
        structTags = Converter.extractVars(structTags,['tag1','tag2'])
        # Gather into structured patches [[[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac]]
        infos = Connector.gatherMatching(structWins,structTags,typeOfWins,structIndirBlkOfWins,
                                         dimsI, dimsJ, dimsK, dim, tol)
        for info in infos:
            noz1 = info[0][0]; noz2 = info[0][1]
            range1 = info[1]; range2 = info[2]
            topp0 = info[3]
            dimZ = Internal.getZoneDim(zones[noz1])
            dimzone = dimZ[4]
            if dimzone == 3: topp = [1,2,3]
            else:
                topp = [1,2]
                topp0 = [topp0[0], topp0[1]]

            if topp0[0] > 0: topp[topp0[0]-1] = 1
            else: topp[-topp0[0]-1] = -1
            if topp0[1] > 0: topp[topp0[1]-1] = 2
            else: topp[-topp0[1]-1] = -2
            if dimzone == 3:
                if topp0[2] > 0: topp[topp0[2]-1] = 3
                else: topp[-topp0[2]-1] = -3
            #------------------------------------------
            # addBC2Zone...
            name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
            name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
            C._addBC2Zone(zones[noz1],name1,'BCMatch',range1,zones[noz2],range2,topp0)
            C._addBC2Zone(zones[noz2],name2,'BCMatch',range2,zones[noz1],range1,topp)

            # couplage RANS/laminar ou euler
            model_z1 = model; model_z2 = model
            eq= Internal.getNodeFromName2(zones[noz1], 'GoverningEquations')
            if eq is not None: model_z1 = Internal.getValue(eq)
            eq= Internal.getNodeFromName2(zones[noz2], 'GoverningEquations')
            if eq is not None: model_z2 = Internal.getValue(eq)

            if (model_z1=='NSTurbulent' and model_z1 != model_z2):
                # creation flag pour tranfert rans/LES
                datap1 = numpy.ones(1, numpy.int32)
                datap2 = numpy.ones(1, numpy.int32)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz1], name1), 'RANSLES', 'DataArray_t', datap1)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz2], name2), 'RANSLES', 'DataArray_t', datap2)
                name_extrap = 'RANS_LES%d_%d'%(noz1,noz2)
                C._addBC2Zone(zones[noz1],name_extrap,'BCExtrapolateRANS',range1)

            if (model_z2=='NSTurbulent'  and  model_z1 != model_z2):
                datap1 = numpy.ones(1, numpy.int32)
                datap2 = numpy.ones(1, numpy.int32)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz2], name2), 'RANSLES', 'DataArray_t', datap2)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz1], name1), 'RANSLES', 'DataArray_t', datap1)
                name_extrap = 'RANS_LES%d_%d'%(noz2,noz1)
                C._addBC2Zone(zones[noz2], name_extrap, 'BCExtrapolateRANS', range2)

    return glob

#==============================================================================
# Returns the list of undefined windows as 2D zones
#         the tag defining matching blks and indices
#         the indirection tab to get the blk from which a window comes
#         the type of windows (i=1: 1, i=imax: 2, j=1: 3 etc)
#         the dimensions of the original blocks
#==============================================================================
def getEmptyWindowsInfoNGON__(t, dim=3):
    try: import Transform.PyTree as T
    except: raise ImportError("Connector.PyTree.getEmptyWindowsInfo__ requires Transform.PyTree module.")
    zones = Internal.getZones(t)
    nzones = len(zones)
    indirBlkOfWins=[]; allTags = []
    for noz in xrange(nzones):
        z = zones[noz]
        zp = C.extractVars(z,['CoordinateX','CoordinateY','CoordinateZ'])# pour eviter des subzones de champs inutiles
        dimZ = Internal.getZoneDim(zp)
        if dimZ[3] == 'NGON':
            faceList = C.getEmptyBC(zp, dim=dim)
            if faceList != []:
                for fl in faceList:
                    winp = T.subzone(zp, fl, type='faces')
                    C._initVars(winp,'centers:tag1',-1.) # defines the opposite window
                    C._initVars(winp,'centers:tag2',-2.) # defines the opp index in opp window
                    allTags += [winp]; indirBlkOfWins += [noz]

    return allTags, indirBlkOfWins

#==============================================================================
def getEmptyWindowsInfoStruct__(t, dim=3):
    try:import Transform.PyTree as T
    except:raise ImportError("getEmptyWindowsInfo__ requires Transform.PyTree modules.")
    zones = Internal.getZones(t)
    nzones = len(zones)
    allWins=[]; typeOfWins=[]; indirBlkOfWins=[]; allTags=[]
    dimsI=[]; dimsJ=[]; dimsK=[]
    for noz in xrange(nzones):
        z = zones[noz]
        dimsZ = Internal.getZoneDim(z)
        if dimsZ[0] == 'Structured':
            ni = dimsZ[1]; nj = dimsZ[2]; nk = dimsZ[3]; dimZone = dimsZ[4]
            nic = max(1,ni-1); njc = max(1,nj-1); nkc = max(1,nk-1)
            nicnjc = nic*njc
            dimsI.append(ni); dimsJ.append(nj); dimsK.append(nk)
            ranges = C.getEmptyBC(z, dim=dim)
            if ranges != []:
                locWins = []; locTypes = []; locIndir=[]
                winp = T.subzone(z,(1,1,1),(1,nj,nk))
                win = C.getFields(Internal.__GridCoordinates__,winp)[0]
                locWins.append(win); locTypes.append(1); locIndir.append(noz)

                winp = T.subzone(z,(ni,1,1),(ni,nj,nk))
                win = C.getFields(Internal.__GridCoordinates__,winp)[0]
                locWins.append(win); locTypes.append(2); locIndir.append(noz)

                winp = T.subzone(z,(1,1,1),(ni,1,nk))
                win = C.getFields(Internal.__GridCoordinates__,winp)[0]
                locWins.append(win); locTypes.append(3); locIndir.append(noz)

                winp = T.subzone(z,(1,nj,1),(ni,nj,nk))
                win = C.getFields(Internal.__GridCoordinates__,winp)[0]
                locWins.append(win); locTypes.append(4); locIndir.append(noz)

                if dim == 3:
                    winp = T.subzone(z,(1,1,1),(ni,nj,1))
                    win = C.getFields(Internal.__GridCoordinates__,winp)[0]
                    locWins.append(win); locTypes.append(5); locIndir.append(noz)

                    winp = T.subzone(z,(1,1,nk),(ni,nj,nk))
                    win = C.getFields(Internal.__GridCoordinates__,winp)[0]
                    locWins.append(win); locTypes.append(6); locIndir.append(noz)

                locTags = Converter.node2Center(locWins)
                locTags = Converter.initVars(locTags,'tag1',-2.) # defines the opposite window
                locTags = Converter.initVars(locTags,'tag2',-2.) # defines the opposite index in opposite window
                for r in ranges:
                    imin = r[0]; jmin = r[2]; kmin = r[4]
                    imax = r[1]; jmax = r[3]; kmax = r[5]
                    now = 0
                    if imin == imax:
                        if imin == 1: now = 1
                        else: now = 2
                    elif jmin == jmax:
                        if jmin == 1: now = 3
                        else: now = 4
                    elif kmin == kmax:
                        if kmin == 1: now = 5
                        else: now = 6

                    tag = locTags[now-1]
                    postag = KCore.isNamePresent(tag,'tag1')
                    taga = tag[1][postag,:]
                    imax = min(imax-1,nic)
                    jmax = min(jmax-1,njc)
                    kmax = min(kmax-1,nkc)
                    if dim == 2: kmax = max(1,kmax)

                    if now == 1 or now == 2:
                        for k in xrange(kmin-1,kmax):
                            taga[jmin-1+k*njc:jmax+k*njc] = -1.
                    elif now == 3 or now == 4:
                        for k in xrange(kmin-1,kmax):
                            taga[imin-1+k*nic:imax+k*nic] = -1.
                    elif now == 5 or now == 6:
                        for j in range(jmin-1,jmax):
                            taga[imin-1+j*nic:imax+j*nic] = -1.
                allWins += locWins
                allTags += locTags
                indirBlkOfWins += locIndir
                typeOfWins += locTypes
    return allTags, allWins, indirBlkOfWins, typeOfWins, dimsI, dimsJ, dimsK
#=============================================================================
# Computes the 1-to-1 connectivity between:
# 1. structured zones
# 2. NGON zones
# 3. structured/NGON 
#=============================================================================
def connectMatch(t, tol=1.e-6, dim=3, type='all'):
    a,typen = Internal.node2PyTree(t)
    glob = 0
    if type == 'structured':        
        glob = _connectMatchStruct__(a,tol,dim,glob)
    elif type == 'unstructured':
        glob = _connectMatchNGON__(a,tol,dim,glob)
    elif type == 'hybrid':
        glob = _connectMatchHybrid__(a,tol,dim,glob)
    else: # all
        glob = _connectMatchStruct__(a,tol,dim,glob)
        glob = _connectMatchNGON__(a,tol,dim,glob)
        glob = _connectMatchHybrid__(a,tol,dim,glob)
    return Internal.pyTree2Node(a,typen)

#===============================================================================
# Detection auto des raccords coincidents dans le cas periodique en
# rotation ou translation
# Standard: l'angle de rotation ou le vecteur de translation est oriente de
# l'interface courante a l'interface connectee.
# Dans le standard, si 3 angles de rotations non nuls, la composition se fait
# selon x, puis y, puis z
#===============================================================================
def connectMatchPeriodic(t, rotationCenter=[0.,0.,0.], 
                         rotationAngle=[0.,0.,0.],
                         translation=[0.,0.,0.], tol=1.e-6, dim=3,
                         unitAngle=None):
    a = Internal.copyRef(t)
    a = C.tagDefinedBC(a)
    a = connectMatchPeriodicStruct__(a, rotationCenter, rotationAngle, translation, tol, dim, unitAngle)
    a = connectMatchPeriodicNGON__(a, rotationCenter, rotationAngle, translation, tol, dim, unitAngle)
    C._rmVars(a, 'definedBC')
    return a

def connectMatchPeriodicNGON__(a, rotationCenter, rotationAngle, translation, tol, dim, unitAngle):
    try: import Post.PyTree as P; import Transform.PyTree as T
    except: raise ImportError("connectMatchPeriodicNGON__ requires Transform and Post modules.")

    if unitAngle == CGK.Radian_s: rotationAngle=[v*180./numpy.pi for v in rotationAngle]
    elif unitAngle not in CGK.AngleUnits_l+[None]: raise ValueError('Value for keyword argument "unitAngle" should be %s'%str(CGK.AngleUnits_l))

    typeOfInputNode = Internal.typeOfNode(a)
    if typeOfInputNode == -1: raise ValueError("connectMatchPeriodicNGON__: invalid input node.")
    zones = Internal.getZones(a)
    indirZones = []; zonesU = []
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[3]=='NGON': zonesU.append(z); indirZones.append(noz)
        noz += 1
    if len(zonesU) == 0: return a

    # get exterior faces
    allExtIndices = []; allExtFaces0 = []
    for z in zonesU:
        indicesF = []
        f = P.exteriorFaces(z, indices=indicesF)
        indicesF = indicesF[0]
        C._initVars(f, 'centers:__tag__', 1.)

        # BC classique
        bnds = Internal.getNodesFromType2(z, 'BC_t')
        # BC Match
        bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        # BC Overlap/NearMatch/NoMatch
        bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')

        zp = Internal.copyRef(z)
        C._deleteZoneBC__(zp); C._deleteFlowSolutions__(zp)

        defined = [] # BC deja definies
        for bc in bnds:
            flist = Internal.getNodeFromName1(bc, Internal.__FACELIST__)
            if flist is not None: defined.append(T.subzone(zp, flist[1], type='faces'))
            erange = Internal.getNodeFromName1(bc, Internal.__ELEMENTRANGE__)
            if erange is not None:
                r = erange[1]
                defined.append(C.selectOneConnectivity(zp,range=[r[0,0],r[0,1]]))

        hook = C.createHook(f, 'elementCenters')
        if defined != []:
            tag = Internal.getNodeFromName2(f, '__tag__')[1]
            defined = C.convertArray2NGon(defined)
            defined = T.join(defined)
            id0 = C.identifyElements(hook, defined)
            tag[id0[:]-1] = 0
        sel = P.selectCells2(f, 'centers:__tag__')
        id1 = C.identifyElements(hook,sel)
        id2 = numpy.empty(id1.size, numpy.int32)
        id2[:] = indicesF[id1[:]-1]
        C.freeHook(hook)
        if id2 != []:
            allExtIndices.append(id2)
            C._initVars(sel,'{centers:tag1}=-1.')# defines the opposite window
            C._initVars(sel,'{centers:tag2}=-1.') # defines the opp index in opp window
            allExtFaces0.append(sel)
        else:
            allExtIndices.append(indicesF)
            C._rmVars(f, ['centers:__tag__'])
            C._initVars(f,'{centers:tag1}=-2.')# defines the opposite window
            C._initVars(f,'{centers:tag2}=-2.') # defines the opp index in opp window
            allExtFaces0.append(f)

    # duplicate exterior faces
    infoPer = duplicatePeriodicZones__(allExtFaces0,rotationCenter,rotationAngle,translation,tol,dim)
    nzonesU = len(zonesU)
    typePeriodic = infoPer[0]
    if typePeriodic==1: signT = [-1,1]; signR=[0,0]
    elif typePeriodic==2: signT=[0,0]; signR=[-1,1]
    elif typePeriodic==3: signT = [-1,-1,1,1]; signR=[-1,1,-1,1]
    dupname = 'DUPPER_' # prefix for duplicated zones
    for i in xrange(1, len(infoPer)):
        # renommage des zones dupliquees
        for noz in xrange(nzonesU):
            zname = infoPer[i][noz][0]
            infoPer[i][noz][0] = dupname+zname
        glob = 0
        glob = _connectMatchNGON__(zonesU,tol,dim,glob,allExtFaces=allExtFaces0+infoPer[i],\
                                       allExtIndices=allExtIndices+allExtIndices, periodic=1,\
                                       rotationCenter=rotationCenter,rotationAngle=rotationAngle,\
                                       Translation=translation, signT=signT[i-1],signR=signR[i-1])
        # update centers:tag1
        if glob > 0:
            tag1p = C.getField('centers:tag1',infoPer[i])
            tag1o = C.getField('centers:tag1',allExtFaces0)
            for noz in xrange(len(tag1p)): tag1p[noz][0]='tag1p'
            tag1o = Converter.addVars([tag1o,tag1p])
            tag1o = Converter.initVars(tag1o,'tag1=maximum({tag1},{tag1p})')
            C.setFields(tag1o,allExtFaces0,loc='centers')
            C._rmVars(allExtFaces0,'centers:tag1p')

        infoPer[i] = []
    # update the original tree
    if typeOfInputNode == 1: # une zone
        a = zonesU[0]
    elif typeOfInputNode == 2: # une liste de zones
        noi = 0
        for indirZ in indirZones:
            zones[indirZ] = zonesU[noi]
            noi += 1
        return zones
    else: # base, liste de bases, arbre
        noi = 0
        for indirZ in indirZones:
            parent,d2 = Internal.getParentOfNode(a,zones[indirZ])
            parent[2][d2] = zonesU[noi]
            noi += 1
    return a

def connectMatchPeriodicStruct__(a,rotationCenter,rotationAngle,translation,tol,dim,unitAngle):
    typeOfInputNode = Internal.typeOfNode(a)
    if typeOfInputNode == -1: raise ValueError("connectMatchPeriodicStruct__: invalid input node.")
    zones = Internal.getZones(a)
    indirZones = []; zonesS = []
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[0] == 'Structured': zonesS.append(z); indirZones.append(noz)
        noz += 1
    if len(zonesS) == 0: return a

    if unitAngle in [CGK.Null_s,CGK.Degree_s,None]: rotationAngleD=rotationAngle
    elif unitAngle == CGK.Radian_s: rotationAngleD=[v*180./numpy.pi for v in rotationAngle]
    else: raise ValueError('Value for keyword argument "unitAngle" should be %s'%str(CGK.AngleUnits_l))

    infoPer = duplicatePeriodicZones__(zonesS,rotationCenter,rotationAngleD,translation,tol,dim)
    nzonesS = len(zonesS)
    typePeriodic = infoPer[0]
    if typePeriodic==1: signT = [-1,1]; signR=[0,0]
    elif typePeriodic==2: signT=[0,0]; signR=[-1,1]
    elif typePeriodic==3: signT = [-1,-1,1,1]; signR=[-1,1,-1,1]
    dupname = 'DUPPER_' # prefix for duplicated zones
    for i in xrange(1, len(infoPer)):
        # renommage des zones dupliquees
        for noz in xrange(nzonesS):
            zname = infoPer[i][noz][0]
            infoPer[i][noz][0] = dupname+zname
        glob = 0
        glob = _connectMatchStruct__(zonesS+infoPer[i],tol,dim,glob)
        for noz in xrange(nzonesS):
            _addPeriodicInfoStruct__(zonesS[noz],rotationCenter,rotationAngle,translation,\
                                         signT[i-1],signR[i-1],dupname,unitAngle=unitAngle)
        infoPer[i] = []
    # update the original tree
    if typeOfInputNode == 1: # une zone
        a = zonesS[0]
    elif typeOfInputNode == 2: # une liste de zones
        noi = 0
        for indirZ in indirZones:
            zones[indirZ] = zonesS[noi]
            noi += 1
        return zones
    else: # base, liste de bases, arbre
        noi = 0
        for indirZ in indirZones:
            parent,d2 = Internal.getParentOfNode(a,zones[indirZ])
            parent[2][d2] = zonesS[noi]
            noi += 1
    return a

#===============================================================================
def _addPeriodicInfoStruct__(z,rotationCenter,rotationAngle,translation,signT,signR, dupname='DUPPER_',unitAngle=None):
    nodes = Internal.getNodesFromType2(z,'GridConnectivity1to1_t')
    for info in nodes:
        if len(info[1]) > 8:
            donorNamePref = info[1][0:7]; donorName = info[1][7:]
            if isinstance(donorNamePref, numpy.ndarray): donorNamePref = donorNamePref.tostring()
            if donorNamePref == dupname: # cas periodique
                info[1] = donorName
                info[2].append(['GridConnectivityProperty', None, [], 'GridConnectivityProperty_t'])
                info = info[2][len(info[2])-1]
                info[2].append(['Periodic', None, [], 'Periodic_t'])
                info = info[2][0]
                v = numpy.zeros( (3), numpy.float64)
                v[0] = rotationCenter[0]; v[1] = rotationCenter[1]; v[2] = rotationCenter[2]
                info[2].append(['RotationCenter', v, [], 'DataArray_t'])
                v = numpy.zeros((3), numpy.float64)
                v[0] = signR*rotationAngle[0]; v[1] = signR*rotationAngle[1]; v[2] = signR*rotationAngle[2]
                info[2].append(['RotationAngle', v, [], 'DataArray_t'])
                if unitAngle is not None: CGL.newDimensionalUnits(info[2][-1],[CGK.Null_s]*4+[unitAngle])
                v = numpy.zeros((3), numpy.float64)
                v[0] = signT*translation[0]; v[1] = signT*translation[1]; v[2] = signT*translation[2]
                info[2].append(['Translation', v, [], 'DataArray_t'])

    return None

#===============================================================================
def duplicatePeriodicZones__(t, rotationCenter=[0.,0.,0.], 
                             rotationAngle=[0.,0.,0.],
                             translation=[0.,0.,0.], tol=1.e-6, dim=3):
    try: import Transform.PyTree as T
    except: raise ImportError("connectMatchPeriodic: requires Transform module.")
    a = Internal.copyRef(t)

    typePeriodic = 0 # 1: translation, 2: rotation, 3: les deux
    for i in xrange(3):
        if float(translation[i]) != 0.: typePeriodic = 1; break
    for i in xrange(3):
        if float(rotationAngle[i]) != 0.: typePeriodic += 2; break

    # Periodicite par rotation ou translation separement
    zones = Internal.getZones(a)
    dupZones = None
    if typePeriodic == 0: return None
    elif typePeriodic == 1:
        zonesdupP = T.translate(zones,(translation[0],translation[1],translation[2]))# periodicite en +transVect
        zonesdupM = T.translate(zones,(-translation[0],-translation[1],-translation[2]))# periodicite en -transVect
        dupZones = [typePeriodic,zonesdupP,zonesdupM]
    elif typePeriodic == 2:
        zonesdupP = T.rotate(zones,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),\
                                 (rotationAngle[0],rotationAngle[1],rotationAngle[2]))
        zonesdupM = T.rotate(zones,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),\
                                 (-rotationAngle[0],-rotationAngle[1],-rotationAngle[2]))
        dupZones = [typePeriodic,zonesdupP,zonesdupM]

    elif typePeriodic == 3:
        zonesdup0 = T.translate(zones,(translation[0],translation[1],translation[2]))# periodicite en +transVect

        # 1. Premier sens de translation
        # 1.1 dans le premier sens de rotation
        zonesdupP = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(rotationAngle[0],rotationAngle[1],rotationAngle[2]))
        # 1.2. dans le sens oppose
        zonesdupM = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(-rotationAngle[0],-rotationAngle[1],-rotationAngle[2]))
        dupZones = [typePeriodic,zonesdupP,zonesdupM]

        zonesdup0 = T.translate(zones,(-translation[0],-translation[1],-translation[2]))# periodicite en -transVect
        zonesdupP = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(rotationAngle[0],rotationAngle[1],rotationAngle[2]))
        # 2.2. dans le sens oppose
        zonesdupM = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(-rotationAngle[0],-rotationAngle[1],-rotationAngle[2]))
        dupZones += [zonesdupP,zonesdupM]
    return dupZones

#==============================================================================
# Set degenerated BC
#==============================================================================
def setDegeneratedBC(t,dim=3, tol=1.e-10):
    try: import Generator
    except: raise ImportError("setDegeneratedBC: requires Generator module.")

    a = Internal.copyRef(t)
    # Is t a topTree, a list of bases etc ?
    # listzones = 0: basis or toptree
    # listzones = 1: list of zones
    listzones = 1 # Is a  a list of zones or a basis or a toptree
    toptree = Internal.isTopTree(a)
    if toptree: listzones = 0
    else:
        # a base ou non
        base = Internal.getBases(a)
        if base != []: listzones = 0
        else: # liste zones ou zone ?
            stdNode = Internal.isStdNode(a)
            if stdNode != 0: return a

    allTags,allWins,indirBlkOfWins,typeOfWins,dimsI,dimsJ,dimsK=getEmptyWindowsInfoStruct__(a,dim)
    if allTags != []:
        volumes = Generator.getVolumeMap(allWins)
        allTags = Converter.addVars([allTags,volumes])
        allTags = Connector.identifyDegenerated(allTags, tol)
        allTags = Converter.extractVars(allTags,['tag1','tag2'])
        # Gather degenerated cells into structured patches [[[noz1],[imin1,imax1,...]],... ]
        infos = Connector.gatherDegenerated(allTags,typeOfWins,indirBlkOfWins,dimsI, dimsJ, dimsK, dim)
    else: infos = []
    zones = Internal.getZones(a)
    glob = 0
    for info in infos:
        noz1 = info[0]; range1 = info[1]
        # Get Parent Nodes
        if listzones == 0:
            r1 = Internal.getParentOfNode(a, zones[noz1])
            parent1 = r1[0]; d1 = r1[1]
        # addBC
        name1 = 'degen%d_%d'%(noz1+1,glob); glob+=1
        if dim == 2: zones[noz1] = C.addBC2Zone(zones[noz1],name1,'BCDegeneratePoint',range1)
        else: zones[noz1] = C.addBC2Zone(zones[noz1],name1,'BCDegenerateLine',range1)

        if listzones == 0: parent1[2][d1] = zones[noz1]
    if listzones == 1: a = zones
    return a

#=============================================================================
# Computes the n to p grid connectivity between zones defined in t
#=============================================================================
def connectNearMatch(t, ratio=2,tol=1.e-6, dim=3):
    try: import Generator.PyTree as G; import Transform.PyTree as T
    except: raise ImportError("connectNearMatch: requires Generator and Transform modules.")
    a,typen = Internal.node2PyTree(t)

    allRatios = []
    if not isinstance(ratio, list):
        if dim == 3:
            allRatios.append([ratio,1,1])#ij
            allRatios.append([1,ratio,1])#ik
            allRatios.append([1,1,ratio])#jk
            allRatios.append([ratio,ratio,1])#ij
            allRatios.append([ratio,1,ratio])#ik
            allRatios.append([1,ratio,ratio])#jk
            allRatios.append([ratio,ratio,ratio])#ijk
        else:
            allRatios.append([ratio,1,1])#ij
            allRatios.append([1,ratio,1])#ik
            allRatios.append([ratio,ratio,1])#ij
    else: allRatios = [ratio]

    model ='Euler'
    bases = Internal.getNodesFromType2(t, 'CGNSBase_t')
    if bases != []:
        eq = Internal.getNodeFromName2(bases[0], 'GoverningEquations')
        if eq is not None: model = Internal.getValue(eq)

    glob = 0
    zones = Internal.getZones(a)
    nzones = len(zones)
    for ratios in allRatios:
        zones2 = []
        for noz1 in xrange(nzones):
            z1 = zones[noz1]
            # modification des GC pour garder les frontieres definies par une BC lors du passage en oneovern
            bcmatch = Internal.getNodesFromType2(z1, 'GridConnectivity1to1_t')
            for i in bcmatch:
                r = Internal.getNodeFromName1(i, 'PointRange')
                if r is not None:
                    ranger = r[1]
                    w = Internal.range2Window(ranger)
                    z1 = C.addBC2Zone(z1, 'dummy', 'BCWall', w)

            bcnearmatch = Internal.getNodesFromType2(z1, 'GridConnectivity_t')
            for i in bcnearmatch:
                type = Internal.getNodeFromName1(i, 'GridConnectivityType')
                if type is not None:
                    val = Internal.getValue(type)
                    if (val == 'Abutting'):
                        r = Internal.getNodeFromName1(i, 'PointRange')
                        ranger = r[1]
                        w = Internal.range2Window(ranger)
                        z1 = C.addBC2Zone(z1, 'dummy', 'BCWall', w)
            # oneovern
            zones2.append(T.oneovern(z1,(ratios[0],ratios[1],ratios[2])))

        # Abutting info for oneovern zones
        allTags2, allWins2, indirBlkOfWins2, typeOfWins2,dimsI2,dimsJ2,dimsK2=getEmptyWindowsInfoStruct__(zones2,dim)
        #
        allTags1, allWins1, indirBlkOfWins1, typeOfWins1,dimsI1,dimsJ1,dimsK1=getEmptyWindowsInfoStruct__(zones,dim)
        #
        ntags1 = len(allTags1)
        for now1 in xrange(ntags1): indirBlkOfWins1[now1] +=nzones
        #
        if len(allTags2) < len(allTags1):
            for iblk in xrange(len(indirBlkOfWins2)):
                if indirBlkOfWins2[iblk] != indirBlkOfWins1[iblk]-len(zones):
                    allTags1.pop(iblk)
                    allWins1.pop(iblk)
                    indirBlkOfWins1.pop(iblk)
                    typeOfWins1.pop(iblk)

        # Identify matching windows of zones with windows of zones2
        if allTags1 != [] and len(allTags1) == len(allTags2):
            allTags = Connector.identifyMatchingNM(allTags2, allTags1,tol)
        else: allTags = []
        allWins=allWins2+allWins1
        typeOfWins=typeOfWins2+typeOfWins1
        indirBlkOfWins=indirBlkOfWins2+indirBlkOfWins1
        dimsI = dimsI2+dimsI1
        dimsJ = dimsJ2+dimsJ1
        dimsK = dimsK2+dimsK1
        #
        # Gather Matching cells into structured patches [ [[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac] ]
        if allTags != []:
            infos = Connector.gatherMatchingNM(allWins, allTags, typeOfWins, indirBlkOfWins, dimsI, dimsJ, dimsK, dim, tol)
        else: infos = []
        for info in infos:
            noz1 = info[0][0]; noz2 = info[0][1]
            dimZ = Internal.getZoneDim(zones[noz1])[4]
            noz2 = noz2-nzones
            dims1 = Internal.getZoneDim(zones[noz1])
            ni1=dims1[1]; nj1=dims1[2]; nk1=dims1[3]
            dimsp1 = Internal.getZoneDim(zones2[noz1])
            nip1=dimsp1[1]; njp1=dimsp1[2]; nkp1=dimsp1[3]
            now1 = info[4][0]; now2 = info[4][1]
            if now1 != now2 and noz1 != noz2:
                range1 = info[1]; range2 = info[2]
                topp0 = info[3];
                now1 = info[4][0]; now2 = info[4][1]
                if (dimZ == 3): topp = [1,2,3]
                else:
                    topp = [1,2]
                    topp0 = [topp0[0], topp0[1]]

                if (topp0[0] > 0): topp[topp0[0]-1] = 1
                else: topp[-topp0[0]-1] = -1
                if (topp0[1] > 0): topp[topp0[1]-1] = 2
                else: topp[-topp0[1]-1] = -2
                if (dimZ == 3):
                    if (topp0[2] > 0): topp[topp0[2]-1] = 3
                    else: topp[-topp0[2]-1] = -3
                #
                # addBC2Zone...
                name1 = 'nmatch%d_%d'%(noz1+1,glob); glob+=1
                name2 = 'nmatch%d_%d'%(noz2+1,glob); glob+=1
                # Get Parent Nodes
                r1 = Internal.getParentOfNode(a, zones[noz1]); parent1 = r1[0]; d1 = r1[1]
                r2 = Internal.getParentOfNode(a, zones[noz2]); parent2 = r2[0]; d2 = r2[1]
                # addBCNearMatch
                rangenm1 = Internal.window2Range(range1)# copy
                if ratios[0] != 1: rangenm1 = G.refineBCRanges__(rangenm1, nip1, njp1, nkp1, ni1, nj1, nk1, 1, ratios[0])
                if ratios[1] != 1: rangenm1 = G.refineBCRanges__(rangenm1, nip1, njp1, nkp1, ni1, nj1, nk1, 2, ratios[1])
                if ratios[2] != 1: rangenm1 = G.refineBCRanges__(rangenm1, nip1, njp1, nkp1, ni1, nj1, nk1, 3, ratios[2])
                rangenm1 = Internal.range2Window(rangenm1)
                zones[noz1] = C.addBC2Zone(zones[noz1],name1,'BCNearMatch',rangenm1,zones[noz2],range2  , topp0)
                zones[noz2] = C.addBC2Zone(zones[noz2],name2,'BCNearMatch',range2  ,zones[noz1],rangenm1, topp )
 

                #couplage RANS/laminar ou euler
                model_z1 = model; model_z2 = model
                eq= Internal.getNodeFromName2(zones[noz1], 'GoverningEquations')
                if eq is not None: model_z1 = Internal.getValue( eq )
                eq= Internal.getNodeFromName2(zones[noz2], 'GoverningEquations')
                if eq is not None: model_z2 = Internal.getValue( eq )

                if (model_z1=='NSTurbulent'  and  model_z1 != model_z2):
                   #creation flag pour tranfert rans/LES
                   datap1 = numpy.ones(1, numpy.int32)
                   datap2 = numpy.ones(1, numpy.int32)
                   Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz1], name1) , 'RANSLES', 'DataArray_t', datap1)
                   Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz2], name2) , 'RANSLES', 'DataArray_t', datap2)
                   name1 = 'RANS_LES%d_%d'%(noz1,noz2)
                   C._addBC2Zone(zones[noz1],name1,'BCExtrapolateRANS',rangenm1)

                if (model_z2=='NSTurbulent'  and  model_z1 != model_z2):
                   datap1 = numpy.ones(1, numpy.int32)
                   datap2 = numpy.ones(1, numpy.int32)
                   Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz2], name2) , 'RANSLES', 'DataArray_t', datap2)
                   Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz1], name1) , 'RANSLES', 'DataArray_t', datap1)
                   name2 = 'RANS_LES%d_%d'%(noz2,noz1)
                   C._addBC2Zone(zones[noz2],name2,'BCExtrapolateRANS',range2  )
                  
                parent1[2][d1] = zones[noz1]; parent2[2][d2] = zones[noz2]

    return Internal.pyTree2Node(a,typen)

#=============================================================================
# Retourne la liste des domaines d'interpolation et la liste des celln associes
# en centres pour une zone noz0 donnee d une base nob0 donnee a partir des
# autres bases. Si sameBase=1: recherche aussi sur sa base
# Si periodicite, il faut avoir duplique les zones periodiques avant
# car la recherche des zones periodiques se fait par 'TempPeriodicZone'
# Fonction appelee par setInterpolations
#=============================================================================
def getInterpolationDomains__(bases, nob0, noz0, zn0, surfacest=[],
                              sameBase=0, intersectionDict=None):
    surf = 0; surfaces = []; periodiczones = []
    if surfacest != []: surf = 1
    nbases = len(bases)
    names = []; zones = []; cellns = []
    for nob in xrange(nbases):
        if nob != nob0:
            nodes = Internal.getNodesFromType1(bases[nob], 'Zone_t')
            for noz in xrange(len(nodes)):
                zname = nodes[noz][0]; intersect = 1
                zn = C.getFields(Internal.__GridCoordinates__,nodes[noz])[0]
                if intersectionDict is not None:
                    if zname not in intersectionDict: intersect = 0

                if intersect == 1:
                    cellN = C.getField('centers:cellN', nodes[noz])[0]
                    donorName = nodes[noz][0]
                    dupzoneinfo = Internal.getNodeFromName(nodes[noz],'TempPeriodicZone')
                    isper = 0
                    if dupzoneinfo is not None:
                        # periodicite par rotation ? On regarde le signe de la rotation
                        rotinfo = Internal.getNodeFromName(dupzoneinfo,'SignOfRotationAngle')
                        if rotinfo is not None:
                            sign = Internal.getValue(rotinfo)
                            if sign == 1: isper=2
                            else: isper=3
                        donorName = Internal.getValue(dupzoneinfo)
                    periodiczones.append(isper)
                    names.append(donorName); zones.append(zn); cellns.append(cellN)
                    if (surf == 1): surfaces.append(surfacest[nob][noz])
        else:
            if sameBase == 1:
                nodes = Internal.getNodesFromType1(bases[nob], 'Zone_t')
                for noz in xrange(len(nodes)):
                    if (noz != noz0):
                        zname = nodes[noz][0]; intersect = 1
                        zn = C.getFields(Internal.__GridCoordinates__,nodes[noz])[0]
                        if intersectionDict is not None:
                            if zname not in intersectionDict: intersect = 0

                        if intersect == 1:
                            cellN = C.getField('centers:cellN', nodes[noz])[0]
                            donorName = nodes[noz][0]
                            dupzoneinfo = Internal.getNodeFromName(nodes[noz],'TempPeriodicZone')
                            if dupzoneinfo is not None:
                                donorName = Internal.getValue(dupzoneinfo)
                                periodiczones.append(1)
                            else: periodiczones.append(0)
                            names.append(donorName); zones.append(zn); cellns.append(cellN)
                            if surf == 1: surfaces.append(surfacest[nob][noz])
    return [names, zones, cellns, surfaces, periodiczones]

#=============================================================================
# blank intersecting and negative volume cells in a mesh. Useful for
# strand grids.
# Set the cellN to 0 for invalid cells, 2 for neighbouring cells
#=============================================================================
def blankIntersectingCells(t, tol=1.e-10, depth=2):
    """Set the cellN at 0 for intersecting cells in the normal direction
    to the wall.
    Usage: blankIntersectingCells(t, tol, depth)"""
    a = Internal.copyRef(t)
    a = addCellN__(a, loc='centers')
    coords = C.getFields(Internal.__GridCoordinates__, a)
    cellN = C.getField('centers:cellN', a)
    res = Connector.blankIntersectingCells(coords, cellN, tol)
    C.setFields(res, a, 'centers')
    return a

#==============================================================================
# Masquage
#==============================================================================
def blankCells(t, bodies, blankingMatrix=[], depth=2,
               blankingType='cell_intersect', delta=1.e-10, dim=3,
               tol=1.e-8, XRaydim1=1000, XRaydim2=1000):
    try: import Transform as T
    except: raise ImportError("blankCells: requires Transform module.")
    if (depth != 1 and depth != 2):
        print 'Warning: blankCells: depth must be equal to 1 or 2. Set to default value (2).'
        depth = 2

    if (blankingType != 'cell_intersect' and \
        blankingType != 'cell_intersect_opt' and \
        blankingType != 'center_in' and \
        blankingType != 'node_in'):
        print 'Warning: blankCells: blankingType must be cell_intersect, cell_intersect_opt, center_in or node_in.'
        print 'Set to default (cell_intersect).'
        blankingType = 'cell_intersect'

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    elif blankingType == 'cell_intersect_opt':
        if depth == 2: blankType = -2
        else: blankType = -1
    else:
        print 'Warning: blankCells: blankingType must be cell_intersect, cell_intersect_opt, center_in or node_in.'
        print 'Set to default (cell_intersect).'
        blankType = 1
                
    nb = 0
    a = Internal.copyRef(t)    
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    a = addCellN__(a, loc=loc)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCells: no basis found in input tree.")

    if blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(bodies)), numpy.int32)
    for b in bases:
        coords = C.getFields(Internal.__GridCoordinates__, b)
        if coords != []:
            if loc == 'centers': cellN = C.getField('centers:cellN', b)
            else: cellN = C.getField('cellN', b)
            for nb2 in xrange(len(bodies)):
                blanking = blankingMatrix[nb, nb2]
                if (bodies[nb2] != [] and \
                    (blanking == 1 or blanking == -1)):
                    bc = []
                    for z in bodies[nb2]:
                        c = C.getFields(Internal.__GridCoordinates__, z)
                        if c != []:
                            c = c[0]
                            if len(c) == 5: # structure
                                # pour le 2D
                                if c[2] == 2: c = T.reorder(c, (-3,1,2))
                                elif c[3] == 2: c = T.reorder(c, (1,-3,2))
                            bc.append(c)
                    masknot = 0
                    if blanking == -1: masknot = 1
                    cellN = Connector.blankCells(
                        coords, cellN, bc, blankingType=blankType, \
                        delta=delta, dim=dim, masknot=masknot, tol=tol, \
                        XRaydim1=XRaydim1, XRaydim2=XRaydim2)
            C.setFields(cellN, b, loc, False)
        nb += 1
    return a

#==============================================================================
# Masquage par Tetra
#==============================================================================
def blankCellsTetra(t, mT4, blankingMatrix=[], blankingType='node_in',
                    tol=1.e-12, cellnval=0, overwrite=0):
    try: import Transform as T
    except: raise ImportError("blankCells: requires Transform module.")

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    else:
        print 'Warning: blankCellsTetra: blankingType must be cell_intersect, center_in or node_in.'
        print 'Set to default (cell_intersect).'
        blankType = 1

    nb = -1
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    a = addCellN__(a, loc=loc)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCellsTetra: no basis found in input tree.")

    if blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(mT4)), numpy.int32)
    for b in bases:
        nb += 1
        coords = C.getFields(Internal.__GridCoordinates__, b)
        if coords == []: continue

        if len(coords[0]) == 5: coords = Converter.convertArray2Hexa(coords) # STRUCT -> HEXA
        else:
          if len(coords[0]) == 4:
            if coords[0][3] == 'NGON': coords = Converter.convertArray2Hexa(coords)
        #print 'TYPE IS %s' %(coords[0][3])

        if loc == 'centers': cellN = C.getField('centers:cellN', b)
        else: cellN = C.getField('cellN', b)
        bc = []
        for nb2 in xrange(len(mT4)):
            blanking = blankingMatrix[nb, nb2]
            #if (mT4[nb2] == []): print 'empty'
            if (mT4[nb2] == [] or \
                (blanking != 1 and blanking != -1)): continue
            i = 0
            for z in mT4[nb2]:
                c = C.getFields(Internal.__GridCoordinates__, z)
                if c != []:
                    c = c[0]
                    bc.append(c)
                i += 1
        if bc == []:
            #print 'Warning : nothing to mask for base %d'%(nb)
            continue
        bc = T.join(bc)
        cellN = Connector.blankCellsTetra(
                coords, cellN, bc, blankingType=blankType, tol=tol, cellnval=cellnval, overwrite=overwrite)
        bc = None
        coords = None
        C.setFields(cellN, b, loc, False)
    return a

#==============================================================================
# Masquage par Tri (surface Tri)
#==============================================================================
def blankCellsTri(t, mT3, blankingMatrix=[], blankingType='node_in',
                    tol=1.e-12, cellnval=0, overwrite=0):
    try: import Transform as T
    except: raise ImportError("blankCells: requires Transform module.")

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    else:
        print 'Warning: blankCellsTri: blankingType must be cell_intersect, center_in or node_in.'
        print 'Set to cell_intersect.'
        blankType = 1

    nb = -1
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    a = addCellN__(a, loc=loc)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCellsTri: no basis found in input tree.")

    if blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(mT3)), numpy.int32)
    for b in bases:
        nb += 1
        coords = C.getFields(Internal.__GridCoordinates__, b)
        if coords == []: continue

        if len(coords[0]) == 5: coords = Converter.convertArray2Hexa(coords) # STRUCT -> HEXA
        else:
          if len(coords[0]) == 4:
            if coords[0][3] == 'NGON': coords = Converter.convertArray2Hexa(coords)
        #print 'TYPE IS %s' %(coords[0][3])

        if loc == 'centers': cellN = C.getField('centers:cellN', b)
        else: cellN = C.getField('cellN', b)
        bc = []
        for nb2 in xrange(len(mT3)):
            blanking = blankingMatrix[nb, nb2]
            #if (mT3[nb2] == []): print 'empty'
            if (mT3[nb2] == [] or \
                (blanking != 1 and blanking != -1)): continue
            i = 0
            for z in mT3[nb2]:
                c = C.getFields(Internal.__GridCoordinates__, z)
                if c != []:
                    c = c[0]
                    bc.append(c)
                i += 1
        if bc == []:
            #print 'Warning : nothing to mask for base %d'%(nb)
            continue
        bc = Converter.convertArray2Tetra(bc); bc = T.join(bc);
        cellN = Connector.blankCellsTri(
                coords, cellN, bc, blankingType=blankType, tol=tol, cellnval=cellnval, overwrite=overwrite)
        bc = None
        coords = None
        C.setFields(cellN, b, loc, False)
    return a


#==============================================================================
# optimisation du recouvrement
#==============================================================================
def optimizeOverlap(t, double_wall=0, priorities=[], intersectionsDict=None):
    try: import Generator.PyTree as G
    except: raise ImportError, 'optimizeOverlap requires Generator module.'
    try: import Post.PyTree as P
    except: raise ImportError, 'optimizeOverlap requires Post module.'
    if double_wall == 1: import DoubleWall
    tol = 1.e-10
    a = Internal.copyRef(t)

    #=====================================================
    # 1-ajout du celln si n'existe pas pour une zone
    #=====================================================
    a = addCellN__(a, loc='centers')
    a = G.getVolumeMap(a)
    bases = Internal.getBases(a)
    nbases = len(bases)

    #=====================================================
    # 2-traitement priorites
    #=====================================================
    nprios = len(priorities)/2
    prios = []
    size = 0
    if nprios == 0:
        for nob in xrange(nbases): prios.append(0)
    else:
        max_prio = 0
        for nob in xrange(nbases):
            baseName = bases[nob][0]
            prio = -1
            for nop in xrange(nprios):
                if priorities[2*nop] == baseName:
                    prio = priorities[2*nop+1]
                    max_prio = max(max_prio, prio)
                    break
            prios.append(prio)

        max_prio += 1
        nprios = len(prios)
        for nop in xrange(nprios):
            if prios[nop] == -1: prios[nop] = max_prio

    #=======================================================================
    # 2 - Recherche des periodicites:
    #     duplication des blocs periodiques dans les bases associees
    #     creation d'un noeud fils au niveau de la zone dupliquee de nom 'TemporaryPeriodicZone'
    #=======================================================================
    for nob in xrange(nbases):
        parentb,db = Internal.getParentOfNode(a,bases[nob])
        bases[nob] = C.duplicatePeriodicZones__(bases[nob])
        parentb[2][db] = bases[nob]

    # Updates the intersection Dictionnary with new duplicated temporary zones:
    if not intersectionsDict:
        intersectionsDict = getIntersectingDomains(a, method='AABB')
    else:
        NewAndOldZones = Internal.getZones(a)
        #NAOZqty = len(NewAndOldZones)
        NewAndOldZonesNames = []
        #[NewAndOldZonesNames.append(NewAndOldZones[i][0]) for i in xrange(NAOZqty)]
        [NewAndOldZonesNames.append(z[0]) for z in NewAndOldZones]
        OldNames = intersectionsDict.keys()
        NewNames = []
        #for i in xrange(NAOZqty):
        #    if NewAndOldZonesNames[i] not in OldNames: NewNames.append(NewAndOldZonesNames[i])
        for zn in NewAndOldZonesNames:
            if zn not in OldNames: NewNames.append(zn)
        if not not NewNames:
            tDuplZones = C.newPyTree(['DuplicatedZones'])
            #for i in xrange(len(NewNames)): tDuplZones[2][1][2].append(Internal.getNodeFromName2(a,NewNames[i]))
            for newname in NewNames: tDuplZones[2][1][2].append(Internal.getNodeFromName2(a,newname))
            periodicZonesIntersectionDict = getIntersectingDomains(tDuplZones,a,method='hybrid')
            intersectionsDict.update(periodicZonesIntersectionDict)

    # Creation of extended centers meshes
    allExtCenters = []; allCenters = []; zones = [];
    for nob1 in xrange(nbases):
        zones.append(Internal.getNodesFromType1(bases[nob1], 'Zone_t'))
        nodesPerBase = C.getFields(Internal.__GridCoordinates__,zones[nob1])
        if nodesPerBase == []:
            allExtCenters.append([[]]*len(zones[nob1]))
            allCenters.append([[]]*len(zones[nob1]))
        else:
            allExtCenters.append(Converter.node2ExtCenter(nodesPerBase))
            allCenters.append(Converter.node2Center(nodesPerBase))
        nodesPerBase = []

    #=======================================================
    # 4-Donor cell search: bbox intersection + adt creation
    #=======================================================
    # on cree par zone de chq base la liste des noms des domaines intersectants
    nobOfIntersectBasesAndZones = []; allHooks = []
    for nob1 in xrange(nbases):
        nobOfIntersectBasesAndZonesForBase1 = []
        allHooksForBase1 = []
        for noz1 in xrange(len(zones[nob1])):
            isIntersected = False
            nobOfIntersectBasesAndZonesForZone1=[]
            for nob2 in xrange(nbases):
                if nob2 != nob1 :
                    nobOfIntersectZonesOfBase2 = [] # numero des zones de base2 intersectant z1
                    for noz2 in xrange(len(zones[nob2])):
                        if zones[nob2][noz2][0] in intersectionsDict[zones[nob1][noz1][0]]:
                            isIntersected = True
                            nobOfIntersectZonesOfBase2.append(noz2)
                    # par base2 opposee, dit si intersecte et les numeros dans base2 des zones intersectantes
                    nobOfIntersectBasesAndZonesForZone1 += [nob2, nobOfIntersectZonesOfBase2]
            nobOfIntersectBasesAndZonesForBase1.append(nobOfIntersectBasesAndZonesForZone1)
            if isIntersected:
                ae1 = allExtCenters[nob1][noz1]
                hook = Converter.createHook([ae1],'extractMesh')
            else: hook = None
            allHooksForBase1.append(hook)
        allHooks.append(allHooksForBase1)
        nobOfIntersectBasesAndZones.append(nobOfIntersectBasesAndZonesForBase1)

    #=====================================================
    # 5-optimisation du recouvrement
    #=====================================================
    isDW = 0
    if double_wall == 0:
        for nob1 in xrange(nbases-1):
            base1 = bases[nob1]
            zones1 = Internal.getNodesFromType1(base1, 'Zone_t')
            prio1 = prios[nob1]; noz1 = 0
            for noz1 in xrange(len(zones1)):
                z1 = zones1[noz1]
                isTempPeriodicZone1 = 0
                if Internal.getNodeFromName(z1,'TempPeriodicZone') == []:
                    r1 = Internal.getParentOfNode(a, z1); parent1 = r1[0]; d1 = r1[1]
                else: isTempPeriodicZone1 = 1
                ae1 = allExtCenters[nob1][noz1]
                ac1 = allCenters[nob1][noz1]
                sol1 = C.getField('centers:cellN', z1)[0]
                vol1 = C.getField('centers:vol', z1)[0]
                ac1 = Converter.addVars([ac1,sol1,vol1])
                adt1 = allHooks[nob1][noz1]
                nobOfIntersectBasesAndZonesForZone1 = nobOfIntersectBasesAndZones[nob1][noz1]
                for nobi in xrange(len(nobOfIntersectBasesAndZonesForZone1)/2):
                    nob2 = nobOfIntersectBasesAndZonesForZone1[nobi*2]
                    if nob2 > nob1:
                        prio2 = prios[nob2]
                        base2 = bases[nob2]
                        zones2 = Internal.getNodesFromType1(base2,'Zone_t')
                        nobOfIntersectZones2 = nobOfIntersectBasesAndZonesForZone1[nobi*2+1]
                        for noz2 in nobOfIntersectZones2:
                            z2 = zones2[noz2]
                            adt2 = allHooks[nob2][noz2]
                            isTempPeriodicZone2 = 0
                            if Internal.getNodeFromName(z2,'TempPeriodicZone') == []:
                                r2 = Internal.getParentOfNode(a, z2); parent2 = r2[0]; d2 = r2[1]
                            else: isTempPeriodicZone2 = 1
                            ae2 = allExtCenters[nob2][noz2]
                            ac2 = allCenters[nob2][noz2]
                            sol2 = C.getField('centers:cellN',z2)[0]
                            vol2 = C.getField('centers:vol',z2)[0]
                            ac2 = Converter.addVars([ac2,sol2,vol2])
                            res = Connector.optimizeOverlap__(ae1, ac1, ae2, ac2, prio1=prio1, prio2=prio2, \
                                                                    isDW=isDW, hook1=adt1, hook2=adt2)
                            cellN1 = Converter.extractVars(res[0],['cellN'])
                            cellN2 = Converter.extractVars(res[1],['cellN'])
                            C.setFields([cellN1], z1, 'centers', False)
                            C.setFields([cellN2], z2, 'centers', False)
                            if isTempPeriodicZone1==0: parent1[2][d1] = z1
                            if isTempPeriodicZone2==0: parent2[2][d2] = z2

    else: #double wall
        dwInfo = DoubleWall.extractDoubleWallInfo__(a)
        firstWallCenters = dwInfo[0]; surfacesExtC = dwInfo[1]
        # liste des surfaces en centres etendus de toutes les zones de toutes les bases
        for nob1 in xrange(nbases-1):
            base1 = bases[nob1]
            zones1 = Internal.getNodesFromType1(base1, 'Zone_t')
            prio1 = prios[nob1]
            for noz1 in xrange(len(zones1)):
                z1 = zones1[noz1]
                isTempPeriodicZone1 = 0
                if Internal.getNodeFromName(z1,'TempPeriodicZone') == []:
                    r1 = Internal.getParentOfNode(a, z1); parent1 = r1[0]; d1 = r1[1]
                else: isTempPeriodicZone1 = 1
                ae1 = allExtCenters[nob1][noz1]
                ac1 = allCenters[nob1][noz1]
                sol1 = C.getField('centers:cellN',z1)[0]
                vol1 = C.getField('centers:vol',z1)[0]
                ac1 = Converter.addVars([ac1,sol1,vol1])
                adt1 = allHooks[nob1][noz1]

                firstWallCenters1 = firstWallCenters[nob1][noz1]
                surfacesExtC1 = surfacesExtC[nob1][noz1]

                # parcours des bases intersectantes
                nobOfIntersectBasesAndZonesForZone1 = nobOfIntersectBasesAndZones[nob1][noz1]
                for nobi in xrange(len(nobOfIntersectBasesAndZonesForZone1)/2):
                    nob2 = nobOfIntersectBasesAndZonesForZone1[nobi*2]
                    if nob2 > nob1:
                        prio2 = prios[nob2]
                        base2 = bases[nob2]
                        zones2 = Internal.getNodesFromType1(base2,'Zone_t')
                        nobOfIntersectZones2 = nobOfIntersectBasesAndZonesForZone1[nobi*2+1]
                        for noz2 in nobOfIntersectZones2:
                            z2 = zones2[noz2]
                            isTempPeriodicZone2 = 0
                            if Internal.getNodeFromName(z2,'TempPeriodicZone') == []:
                                r2 = Internal.getParentOfNode(a, z2); parent2 = r2[0]; d2 = r2[1]
                            else: isTempPeriodicZone2 = 1
                            ae2 = allExtCenters[nob2][noz2]
                            ac2 = allCenters[nob2][noz2]
                            sol2 = C.getField('centers:cellN',z2)[0]
                            vol2 = C.getField('centers:vol',z2)[0]
                            ac2 = Converter.addVars([ac2,sol2,vol2])
                            adt2 = allHooks[nob2][noz2]
                            #
                            isDW=0
                            firstWallCenters2 = firstWallCenters[nob2][noz2]
                            surfacesExtC2 = surfacesExtC[nob2][noz2]
                            if firstWallCenters1 != [] and firstWallCenters2 != []:
                                isDW = 1
                                acp1 = Converter.initVars(ac1,'{cellN}=minimum(2.,2*{cellN})')
                                acp2 = Converter.initVars(ac2,'{cellN}=minimum(2.,2*{cellN})')
                                acn1 = Connector.changeWall__(acp1, firstWallCenters1, surfacesExtC2)
                                acn2 = Connector.changeWall__(acp2, firstWallCenters2, surfacesExtC1)
                                cellN1 = Converter.extractVars(ac1,['cellN','vol'])
                                cellN2 = Converter.extractVars(ac2,['cellN','vol'])
                                acn1 = Converter.addVars([acn1,cellN1])
                                acn2 = Converter.addVars([acn2,cellN2])
                                res = Connector.optimizeOverlap__(ae1, acn1, ae2, acn2, prio1,prio2,isDW,\
                                                                        hook1=adt1, hook2=adt2)
                                cellN1 = Converter.extractVars(res[0],['cellN'])
                                cellN2 = Converter.extractVars(res[1],['cellN'])
                                ac1 = Converter.rmVars(ac1,['cellN']); ac1 = Converter.addVars([ac1,cellN1])
                                ac2 = Converter.rmVars(ac2,['cellN']); ac2 = Converter.addVars([ac2,cellN2])
                                C.setFields([cellN1], z1, 'centers', False)
                                C.setFields([cellN2], z2, 'centers', False)
                                if isTempPeriodicZone1==0: parent1[2][d1] = z1
                                if isTempPeriodicZone2==0: parent2[2][d2] = z2
                            else:
                                res = Connector.optimizeOverlap__(ae1, ac1, ae2, ac2, prio1,prio2,isDW,\
                                                                        hook1=adt1, hook2=adt2)
                                cellN1 = Converter.extractVars(res[0],['cellN'])
                                cellN2 = Converter.extractVars(res[1],['cellN'])
                                ac1 = Converter.rmVars(ac1,['cellN']); ac1 = Converter.addVars([ac1,cellN1])
                                ac2 = Converter.rmVars(ac2,['cellN']); ac2 = Converter.addVars([ac2,cellN2])
                                C.setFields([cellN1], z1, 'centers', False)
                                C.setFields([cellN2], z2, 'centers', False)

                                if isTempPeriodicZone1==0: parent1[2][d1] = z1
                                if isTempPeriodicZone2==0: parent2[2][d2] = z2
    # Delete duplicated periodic zones
    a = C.removeDuplicatedPeriodicZones__(a)
    #
    C._rmVars(a, 'centers:vol')
    # remise a jour du celln : de 3 passe a 2
    C._initVars(a,'{centers:cellN}=minimum(2.,{centers:cellN})')
    #
    for nob1 in xrange(len(allHooks)):
        allHooksForBase1 = allHooks[nob1]
        for hookZ1 in allHooksForBase1:
            if hookZ1 is not None: Converter.freeHook(hookZ1)
    return a

#==============================================================================
def maximizeBlankedCells(t, depth=2, dir=1, cellNName='centers:cellN'):
    var = cellNName.split(':')
    if len(var)==2 :
        if var[0] == 'centers': loc = 'centers'
        varC = var[1]
    else: 
        loc = 'nodes'; varC = var[0]

    ghost = Internal.getNodeFromName3(t, 'ZoneRind')
    if ghost is None:
        a = Internal.addGhostCells(t, t, depth, modified=[cellNName])
    else: a = Internal.copyRef(t)
    cellN = C.getField(cellNName, a)
    cellN = Connector.maximizeBlankedCells(cellN, depth, dir, cellNName=varC)
    C.setFields(cellN, a, loc, False)
    if ghost is None:
        a = Internal.rmGhostCells(a, a, depth, modified=[cellNName])
    return a

#==============================================================================
# Apply overlap BCs on cell nature field inside zone
# compatible avec une rangee de cellules d'interpolation
# Seulement pour les grilles structurees (no check)
#==============================================================================
def _applyBCOverlapsStructured(z, depth, loc):
    varc = 'cellN'
    if loc == 'centers': varc = 'centers:cellN'; shift = 0
    else: shift = 1
    cellN = C.getField(varc, z)[0]
    ni = cellN[2]; nj = cellN[3]; nk = cellN[4]

    overlaps = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for o in overlaps:
        n = Internal.getNodeFromType(o, 'GridConnectivityType_t')
        if n is not None:
            val = Internal.getValue(n)
            if val == 'Overset':
                isDD = 0
                userDef = Internal.getNodesFromName(o, 'UserDefinedData')
                if userDef != []:
                    if len(userDef[0]) == 4:
                        info = userDef[0][2][0]
                        if info[0] == 'doubly_defined': isDD = 1
                if isDD == 0:
                    r = Internal.getNodesFromType(o, 'IndexRange_t')
                    l = Internal.getNodesFromType(o, 'IndexArray_t')
                    if (r == [] and l == []):
                        print "Warning: applyBCOverlaps: BCOverlap is ill-defined."
                    elif r != []:
                        rangew = r[0][1]
                        w = Internal.range2Window(rangew)
                        imin = w[0]; jmin = w[2]; kmin = w[4]
                        imax = w[1]; jmax = w[3]; kmax = w[5]
                        cellN = Connector.applyBCOverlapsStruct__(cellN,(imin,jmin,kmin),(imax,jmax,kmax),depth,loc)
                        C.setFields([cellN], z, loc, False)
    return None

def _applyBCOverlapsUnstructured(z, depth, loc):
    varc = 'cellN'
    if loc == 'centers': varc = 'centers:cellN'
    cellN = C.getField(varc, z)[0]
    zoneBC = Internal.getNodesFromType2(z, 'BC_t')
    for bc in zoneBC:
        val = Internal.getValue(bc)
        if val == 'BCOverlap':
            faceListN = Internal.getNodesFromName(bc, Internal.__FACELIST__)
            if faceListN != []:
                faceList = Internal.getValue(faceListN[0])
                cellN = Connector.applyBCOverlapsNG__(cellN,faceList,depth,loc)
                C.setFields([cellN],z,loc,False)
    return None

def applyBCOverlaps(t, depth=2, loc='centers'):
  a = Internal.copyRef(t)
  # ajout du celln si n'existe pas pour une zone
  a = addCellN__(a, loc=loc)
  zones = Internal.getZones(a)
  for z in zones:
      dimZ = Internal.getZoneDim(z)
      if dimZ[0] == 'Structured': _applyBCOverlapsStructured(z,depth,loc)
      else:
          if dimZ[3] == 'NGON':
              _applyBCOverlapsUnstructured(z,depth,loc)
          else:
              print 'Warning: applyBCOverlaps: only for NGON unstructured zones.'
  return a

#==============================================================================
# IN: a: contains the cellN located at nodes or centers
# IN: depth can be positive or negative
# Return depth layers of interpolated points at the fringe of blanked points
#==============================================================================
def setHoleInterpolatedPoints(a, depth=2, dir=0, loc='centers', 
                              cellNName='cellN'):
    count = 0
    a = setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName)
    count +=1

    ghost = Internal.getNodeFromName(a, 'ZoneRind')
    if loc == 'centers': varcelln = 'centers:'+cellNName
    else: varcelln = cellNName
    if ghost is None:
        a = Internal.addGhostCells(a, a, 2, modified=[varcelln])
        a = setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName)
        a = Internal.rmGhostCells(a,a,2, modified=[varcelln])
    return a

def setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName='cellN'):
    if depth == 0: return a
    if loc == 'centers': varcelln = 'centers:'+cellNName
    else: varcelln = cellNName
    for z in Internal.getZones(a):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Unstructured' and count == 1: pass
        else: # passage ghost cells
            cellN = C.getField(varcelln, z)[0]
            if cellN != []:# cellN existe
                if depth < 0:
                    cellN = Converter.initVars(cellN,'%s = 1-{cellN}+({cellN}>1.5)*3'%cellNName)
                    if loc == 'centers':
                        cellN = Connector.getOversetHolesInterpCellCenters__(cellN, -depth, dir, cellNName)
                    else:
                        cellN = Connector.getOversetHolesInterpNodes__(cellN,-depth, dir, cellNName)
                    cellN = Converter.initVars(cellN, '%s = 1-{cellN}+({cellN}>1.5)*3'%cellNName)

                else: # depth > 0 
                    if loc =='centers': cellN = Connector.getOversetHolesInterpCellCenters__(cellN, depth, dir, cellNName)
                    else: cellN = Connector.getOversetHolesInterpNodes__(cellN, depth, dir, cellNName)
                C.setFields([cellN], z, loc, False)
    return a

#==============================================================================
# IN: ni, nj, nk
# IN: pointList: liste d'index i,j,k
# OUT: liste d'index ind = i + j*ni + k*ni*nj
#==============================================================================
def globalIndex__(ni, nj, nk, pointList):
    s = pointList.shape[1]
    a = numpy.empty((s), dtype=numpy.int32)
    nij = ni*nj
    a[:] = (pointList[0,:]-1) + (pointList[1,:]-1)*ni + (pointList[2,:]-1)*nij
    return a

#==============================================================================
# Cette fonction remplace le noeud oversetHoles
#==============================================================================
def cellN2OversetHoles(t, append=False):
    if append:
        try: import Compressor
        except: raise ImportError("cellN2OversetHoles: requires Compressor module.")
    a = Internal.copyRef(t)
    zones = Internal.getZones(a)
    for z in zones:
        # First find cellN
        cellN = C.getField('cellN', z)[0] # cellN en noeud
        loc = 'nodes'
        if cellN == []:
            loc = 'centers'
            cellN = C.getField('centers:cellN', z)[0] # cellN en centres
        if cellN != []:
            cellN = cellN[1]
            cellN = cellN.reshape(cellN.size, order='Fortran' )
        # dim
        dim = Internal.getZoneDim(z); cellDim = dim[4]
        # Structured zone
        if dim[0] == 'Structured': # output im,jm,km
            im = dim[1]; jm = dim[2]; km = dim[3]
            if loc == 'centers': im = im-1; jm = jm-1; km = km-1
            nb = 0
            for i in xrange(cellN.size):
                if cellN[i] == 0: nb += 1
            pointList = numpy.empty((cellDim,nb), dtype=numpy.int32, order='Fortran' )
            connector.cellN2OversetHolesStruct(pointList, cellN, im, jm, cellDim, cellN.size)
        # Unstructured zone
        else:
            nb = 0
            for i in xrange(cellN.size):
                if cellN[i] == 0: nb += 1
            pointList = numpy.empty((1,nb), dtype=numpy.int32)
            connector.cellN2OversetHolesUnStruct(pointList, cellN, cellN.size)
        if nb != 0:
            # Push OversetHoles Node
            GC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
            if GC == []:
                z[2].append(['ZoneGridConnectivity', None, [],
                             'ZoneGridConnectivity_t'])
                info = z[2][len(z[2])-1]
            else: info = GC[0]

            OH = Internal.getNodesFromType1(info, 'OversetHoles_t')
            if OH == []:
                info[2].append(['OversetHoles', None, [], 'OversetHoles_t'])
                info = info[2][len(info[2])-1]
                if loc == 'centers':
                    v = numpy.fromstring('CellCenter', 'c')
                    info[2].append(['GridLocation', v, [], 'GridLocation_t'])
                info[2].append(['PointList', pointList, [], 'IndexArray_t'])
            else:
                PL = Internal.getNodesFromName(OH[0], 'PointList') # prev
                if PL == []:
                    OH[0][2].append(['PointList', pointList, [],
                                     'IndexArray_t'])
                else:
                    if not append:
                        PL[0][1] = pointList
                        DL = Internal.getNodesFromName(OH[0], 'DeltaList')
                        if DL != []:
                            (p, c) = Internal.getParentOfNode(OH[0], DL[0])
                            del p[2][c]
                    else: # append
                        ref = PL[0][1]
                        PL[0][1] = pointList
                        if dim[0] == 'Structured': # conversion en indices globaux
                            im = dim[1]; jm = dim[2]; km = dim[3]
                            if loc == 'centers': im = im-1; jm = jm-1; km = km-1
                            ref = globalIndex__(im, jm, km, ref)
                            pointList = globalIndex__(im, jm, km, pointList)

                        delta = Compressor.deltaIndex(pointList, ref)
                        DL = Internal.getNodesFromName(OH[0], 'DeltaList')
                        if DL == []:
                            delta = numpy.concatenate([ref, delta])
                            OH[0][2].append(['DeltaList', delta, [],
                                             'IndexArray_t'])
                        else: # concatenate
                            prev = DL[0][1]
                            delta = numpy.concatenate([prev, delta])
                            DL[0][1] = delta
    return a

#==============================================================================
def setOversetHoles(z, ind):
    if z[3] != 'Zone_t': raise TypeError("setOversetHoles: only for a zone.")
    zp = Internal.copyRef(z)
    type = 0
    if isinstance(ind, numpy.ndarray): type = 1
    elif isinstance(ind, list):
        if ind == []: return zp
        if isinstance(ind[0], int): type = 2
        elif isinstance(ind[0], numpy.ndarray): type = 3
        elif isinstance(ind[0], list): type = 4
    else: raise TypeError("setOversetHoles: argument must be a list or a numpy of indices.")

    # cree la deltaList
    if type == 1: pointList = ind
    elif type == 2: pointList = numpy.array(ind, dtype=numpy.int32)
    elif type == 3: pointList = numpy.concatenate(ind)
    elif type == 4:
        b = []
        for i in ind: b.append(numpy.array(ind, dtype=numpy.int32))
        pointList = numpy.concatenate(b)

    # Push OversetHoles Node
    GC = Internal.getNodesFromType1(zp, 'ZoneGridConnectivity_t')
    if GC == []:
        zp[2].append(['ZoneGridConnectivity', None, [],
                      'ZoneGridConnectivity_t'])
        info = zp[2][len(z[2])-1]
    else: info = GC[0]
    OH = Internal.getNodesFromType(info, 'OversetHoles_t')
    if OH == []:
        info[2].append(['OversetHoles', None, [], 'OversetHoles_t'])
        info = info[2][len(info[2])-1]
        v = numpy.fromstring('CellCenter', 'c')
        info[2].append(['GridLocation', v, [], 'GridLocation_t'])
        if (type == 1 or type == 2):
            info[2].append(['PointList', pointList, [], 'IndexArray_t'])
        else: info[2].append(['DeltaList', pointList, [], 'IndexArray_t'])
    else:
        info = OH[0][2][len(info[2])-1]
        v = numpy.fromstring('CellCenter', 'c')
        info[2].append(['GridLocation', v, [], 'GridLocation_t'])
        if (type == 1 or type == 2):
            info[2].append(['PointList', pointList, [], 'IndexArray_t'])
        else: info[2].append(['DeltaList', pointList, [], 'IndexArray_t'])
    return zp

#------------------------------------------------------------------------------
def getBBIntersectingDomains(A, B, tol=1.e-6):
    """Returns the list zones of B that intersect zones of A.
    Usage: getBBIntersectingDomains(A, B)"""
    try: import Generator.PyTree as G
    except: raise ImportError("getBBIntersectingDomains: requires Generator module.")

    BB_A=G.BB(A)
    BB_B=G.BB(B)

    doms = []
    zones2 = Internal.getZones(BB_B)
    for z2 in zones2:
        zones1 = Internal.getZones(BB_A)
        for z1 in zones1:
            if (z1 != z2 and G.bboxIntersection(z1, z2, tol=tol, isBB=True) == 1):
                doms.append(z1)
    return doms

#------------------------------------------------------------------------------
def getIntersectingDomains(t, t2=None, method='AABB', taabb=None, tobb=None,
                           taabb2=None, tobb2=None):
    try: import Generator.PyTree as G
    except: raise ImportError("getIntersectingDomains: requires Generator module.")
    TotInter = 0
    zones = Internal.getZones(t)
    # Gets all zones names and will iterate over them
    zonesNames = []
    for z in zones: zonesNames.append(z[0])
    N = len(zonesNames)
    if not t2:
        zonesNames2 = zonesNames
        zones2 = zones
        N2 = N
    else:
        zones2 = Internal.getZones(t2)
        # Gets all zones names and will iterate over them
        zonesNames2 = []
        for z in zones2: zonesNames2.append(z[0])
        N2 = len(zonesNames2)
    # Preallocates Dictionnary
    IntDict = dict.fromkeys(zonesNames)

    if method == 'AABB':
        if not taabb: taabb = G.BB(t)
        if not t2: taabb2 = taabb
        elif not taabb2: taabb2 = G.BB(t2)

        for z1 in zonesNames:
            aabb1 = Internal.getZones(taabb)            
            aabb1 = Internal.getNodeFromName(aabb1,z1)
            IntDict[z1] = []
            for z2 in zonesNames2:
                if z1 != z2:
                    aabb2 = Internal.getZones(taabb2)
                    aabb2 = Internal.getNodeFromName(aabb2,z2)
                    if (G.bboxIntersection(aabb1, aabb2, tol=1.e-10, isBB=True, method='AABB') == 1):
                        IntDict[z1].append(z2) # saves the intersected zones names
                        TotInter += 1

    elif method == 'OBB':
        if not tobb: tobb = G.BB(t,method='OBB')
        if not t2: tobb2 = tobb
        elif not tobb2: tobb2 = G.BB(t2,method='OBB')
        for z1 in zonesNames:
            obb1 = Internal.getZones(tobb)
            obb1 = Internal.getNodeFromName(obb1,z1)
            IntDict[z1] = []
            for z2 in zonesNames2:
                if z1 != z2:
                    obb2 = Internal.getZones(tobb2)
                    obb2 = Internal.getNodeFromName(obb2,z2)
                    if (G.bboxIntersection(obb1, obb2, isBB=True, method='OBB') == 1):
                        IntDict[z1].append(z2) # saves the intersected zones names
                        TotInter += 1

    elif method == 'hybrid':
        if not taabb: taabb = G.BB(t)
        if not tobb: tobb = G.BB(t,method='OBB')
        if not t2:
            taabb2 = taabb
            tobb2 = tobb
        if not not t2 and not taabb2: taabb2 = G.BB(t2)
        if not not t2 and not tobb2: tobb2 = G.BB(t2,method='OBB')
        for z1 in zonesNames:
            aabb1 = Internal.getZones(taabb)
            aabb1 = Internal.getNodeFromName(aabb1,z1)
            obb1 = Internal.getZones(tobb)
            obb1 = Internal.getNodeFromName(obb1,z1)
            IntDict[z1] = []
            for z2 in zonesNames2:
                if z1 != z2:
                    obb2 = Internal.getZones(tobb2)
                    obb2 = Internal.getNodeFromName(obb2,z2)
                    aabb2 = Internal.getZones(taabb2)
                    aabb2 = Internal.getNodeFromName(aabb2,z2)
                    if (G.bboxIntersection(obb1, obb2, isBB=True,method='OBB') == 0):
                        continue
                    elif (G.bboxIntersection(aabb1, aabb2, tol=1.e-10, isBB=True,method='AABB') == 0):
                        continue
                    elif (G.bboxIntersection(aabb1, obb2, isBB=True,method='AABBOBB') == 0):
                        continue
                    elif (G.bboxIntersection(aabb2, obb1, isBB=True,method='AABBOBB') == 0):
                        continue
                    else:
                        IntDict[z1].append(z2) # saves the intersected zones names
                        TotInter += 1
    else:
        print 'Warning getIntersectingDomains: method',method,'not implemented. Switched to AABB.'
        return getIntersectingDomains(t, method='AABB', taabb=taabb, tobb=tobb)

    print 'Total zone/zone intersections:', TotInter
    return IntDict
#------------------------------------------------------------------------------
def getCEBBIntersectingDomains(basis0, bases0, sameBase=0):
    """Return the list of interpolation domains defined in bases
    for any zone in basis. If sameBase=1 : interpolation domains
    can be zones from the same basis, but different from the current zone.
    Usage: getCEBBIntersectingDomains(basis, bases, sameBase)"""
    try: import Generator.PyTree as G
    except: raise ImportError("getCEBBIntersectingDomains: requires Generator module.")

    basis = Internal.getBases(basis0)
    bases = Internal.getBases(bases0)
    if len(basis) != 1: raise TypeError("getCEBBIntersectingDomains: not a CGNS base.")
    if len(bases) < 1: raise TypeError("getCEBBIntersectingDomains: not a list of CGNS bases.")
    basis = basis[0]
    # precond : recherche des bases intersectantes
    tol = 1e-6
    [xmin1,ymin1,zmin1,xmax1,ymax1,zmax1] = G.bbox(basis)
    basename1 = basis[0]
    bases2 = []

    for b in bases:
        basename2 = b[0]
        if ((sameBase == 1) or (sameBase == 0 and basename1 != basename2)):
            [xmin2,ymin2,zmin2,xmax2,ymax2,zmax2] = G.bbox(b)
            if (xmax1  > xmin2-tol and xmin1 < xmax2+tol and
                ymax1  > ymin2-tol and ymin1 < ymax2+tol and
                zmax1  > zmin2-tol and zmin1 < zmax2+tol):
                bases2.append(b)

    # recherche des zones intersectantes par base
    doms = []
    zones1 = Internal.getNodesFromType1(basis, 'Zone_t')
    for z1 in zones1:
        name1 = z1[0]
        doms1 = []
        for b2 in bases2:
            basename2 = b2[0]
            zones2 = Internal.getNodesFromType1(b2, 'Zone_t')
            if (sameBase == 1 and basename1 == basename2):
                for z2 in zones2:
                    name2 = z2[0]
                    if name1 != name2:
                        if G.CEBBIntersection(z1, z2) == 1:
                            doms1.append(z2[0])
            else:
                for z2 in zones2:
                    if G.CEBBIntersection(z1, z2) == 1:
                        doms1.append(z2[0])
        doms.append(doms1)
    return doms

#==============================================================================
def getCEBBTimeIntersectingDomains(base0, func, bases0, funcs, \
                                   inititer=0, niter=1, dt=1, \
                                   sameBase=0):
    """Return the list domains defined in bases, that intersect in time
    with base. Func defines a python motion of base with respect to time.
    Usage: getCEBBTimeIntersectingDomains(base, func, bases, funcs,
                                          inititer,niter, dt, sameBase)"""
    try: import RigidMotion.PyTree as R
    except: raise ImportError("getCEBBTimeIntersectingDomains: requires RigidMotion module.")
    base = Internal.getBases(base0)
    bases = Internal.getBases(bases0)
    if len(base) != 1: raise TypeError("getCEBBIntersectingDomains: not a CGNS base.")
    if len(bases) < 1: raise TypeError("getCEBBIntersectingDomains: not a list of CGNS bases.")
    basem = base[0]
    zones = Internal.getNodesFromType1(basem,'Zone_t'); nzones = len(zones)
    zones = []
    t = 0.
    doms = []
    for i in xrange(inititer, niter):
        t = i*dt
        if func != []: basep = R.evalPosition(basem, t, func)
        else: basep = Internal.copyRef(basem)
        # mouvement relatif des bases
        nob = 0
        basesp = []
        for b in bases:
            if funcs[nob] != []: bp = R.evalPosition(b, t, funcs[nob])
            else: bp = Internal.copyRef(b)
            basesp.append(bp)
            nob += 1
        domsp = getCEBBIntersectingDomains(basep, basesp, sameBase)
        if domsp != [[]]*nzones:
            # premieres intersections
            if doms == []: doms = domsp
            else:
                nod = 0
                for dom1 in doms:
                    dom2 = domsp[nod]
                    doms1 = []
                    for name2 in dom2:
                        found = 0
                        for name1 in dom1:
                            if name1 == name2:
                                found = 1
                                break

                        if found == 0: doms1.append(name2)
                    dom1 += doms1
                    doms[nod] = dom1
                    nod += 1

    return doms


#=============================================================================
# application de la CL doublement definie:
# remet a 1 le celln si cellule non interpolable
#=============================================================================
def getDoublyDefinedDonorZones__(t,oversetgcnode):
    listOfDnrZones=[]; listOfDnrCellN=[]
    donorNamesA = oversetgcnode[1]
    donorNames=[]; donorName = ''
    for i in donorNamesA:
        if i == ',':
            donorNames.append(donorName)
            donorName = ''
        else: donorName+=i
    # last one
    donorNames.append(donorName)
    # Duplicated periodic zones
    donorNamesPer = []
    for z in Internal.getZones(t):
        dupzoneinfo = Internal.getNodeFromName(z,'TempPeriodicZone')
        if dupzoneinfo is not None:
            donorName = Internal.getValue(dupzoneinfo)
            if donorName in donorNames:
                donorNamesPer.append(z[0])
    donorNames+= donorNamesPer
    for donorName in donorNames:
        dnrZone = Internal.getNodeFromName2(t,donorName)
        if dnrZone is not None:
            coords = C.getFields(Internal.__GridCoordinates__,dnrZone)[0]
            listOfDnrZones.append(coords)
            cellN2 = C.getField('centers:cellN',dnrZone)[0]
            if cellN2 == []:
                print 'Warning: setDoublyDefined: cellN init to 1 for zone %s.'%dnrZone[0]
                C._initVars(dnrZone,'centers:cellN',1.)
                cellN2 = C.getField('centers:cellN',dnrZone)[0]
            listOfDnrCellN.append(cellN2)

    return listOfDnrZones,listOfDnrCellN

#=============================================================================
def setDoublyDefinedBC(t, depth=2):
    a = Internal.copyRef(t)
    a = addCellN__(a, loc='centers')
    C._initVars(a,'centers:cellN_dd',1.)
    #=======================================================================
    # 2 - Recherche des periodicites :
    #     duplication des blocs periodiques dans les bases associees
    #     creation d un noeud fils au niveau de la zone dupliquee de nom 'TemporaryPeriodicZone'
    #=======================================================================
    bases = Internal.getBases(a)
    nbases = len(bases)
    for nob in xrange(nbases):
        parentb,db = Internal.getParentOfNode(a,bases[nob])
        bases[nob] = C.duplicatePeriodicZones__(bases[nob])
        parentb[2][db] = bases[nob]

    zones = Internal.getZones(a)
    for z in zones:
        if Internal.getNodeFromName(z,'TempPeriodicZone') is not None: pass
        else:
            cellNDD = C.getField('centers:cellN',z)[0]
            (parent, d2) = Internal.getParentOfNode(a, z)
            overlaps = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for o in overlaps:
                n = Internal.getNodesFromType(o, 'GridConnectivityType_t')
                if n != []:
                    val = Internal.getValue(n[0])
                    if val == 'Overset':
                        userDef = Internal.getNodesFromName(o, 'UserDefinedData')
                        if userDef != []:
                            if (len(userDef[0]) == 4):
                                info = userDef[0][2][0]
                                if (info[0] == 'doubly_defined'):
                                    r = Internal.getNodesFromType(o,'IndexRange_t')
                                    win = Internal.range2Window(r[0][1])
                                    zone =  C.getFields(Internal.__GridCoordinates__,z)[0]
                                    listOfInterpZones,cellns=getDoublyDefinedDonorZones__(a,o)
                                    # detection des pts non interpolables
                                    imin = win[0]; jmin = win[2]; kmin = win[4]
                                    imax = win[1]; jmax = win[3]; kmax = win[5]
                                    rangew = [int(imin),int(imax),int(jmin),int(jmax),int(kmin),int(kmax)]
                                    cellNDD = Connector.setDoublyDefinedBC(zone, cellNDD, listOfInterpZones, cellns, \
                                                                               rangew, depth)
            cellNDD[0] = 'cellN_dd'
            C.setFields([cellNDD], z, 'centers', False)
            parent[2][d2] = z

    # Delete duplicated periodic zones
    a = C.removeDuplicatedPeriodicZones__(a)

    C._initVars(a,'centers:cellN=minimum({centers:cellN}*{centers:cellN_dd},2.)')
    C._rmVars(a,['centers:cellN_dd'])
    return a

#=============================================================================
# masque XRay: pierce pts
#=============================================================================
def maskXRay__(body, delta=0., dim=3, isNot=0, tol=1.e-8):
    """Create the pierce points of a X-Ray mask defined by body."""
    body = C.convertArray2Tetra(body)
    surf = C.getFields(Internal.__GridCoordinates__, body)
    pts = Connector.maskXRay__(surf, delta, dim, isNot, tol)
    return C.convertArrays2ZoneNode('XRayPts', [pts])

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

#=============================================================================
# Calcul des donnees pour les interpolations IBC
# IN: aR: arbre, zone, bases des receveurs
#     cellN=2 definit les points a corriger notes PC
#     vect(distance) definit le vecteur normal a la paroi * la distance du PC a la paroi
# IN: aD: arbre, zone, bases des donneurs
# IN: order: ordre des interpolations (2, 3, 5)
# IN: method: 'lagrangian', 'leastsquares'
# IN: penalty=1: penalise une cellule donneuse en terme de volume si elle est au bord
# IN: nature=0: aucun sommet de la cellule d'interpolation ne doit etre avec un cellN=0
#     nature=1: tous les sommets de la cellule d interpolation doivent etre de cellN=1
# IN: hook: hook sur l'arbre de recherche (pas reconstruit dans setInterpIBC), l'ordre doit suivre celui de zonesD
# IN: storage: type de stockage (direct: sur le bloc interpole, inverse: sur le bloc d'interpolation)
# IN: loc='nodes','cells','faces': interpolation appliquee pour les receveurs (localises en noeuds/centres/faces)
# IN: hi: decalage des pts IBC interieurs au corps pour avoir le pt interpole - hi peut etre constant ou variable
# IN: he: decalage des pts IBC exterieurs au corps pour avoir le pt interpole - he peut etre constant ou variable
#IN: image=0: algo basique (point image determine par hi,he ou symetrique)
#IN: image=1: algo par projection sur l isosurface de la distance max du front
# OUT: stockage direct: retourne tR avec les donnees IBC
# OUT: stockage indirect: retourne tD avec les donnees IBC
#=============================================================================
def setIBCData(tR, tD, order=2, penalty=0, nature=0,
               method='lagrangian', loc='nodes', storage='direct',
               hook=None, sameName=0, he=0., hi=0., dim=3, image=0):
    try: import ToolboxIBM as IBM
    except: raise ImportError, "setIBCData requires ToolboxIBM module."

    locR = loc
    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)

    # Si pas de cellN receveur, on retourne
    if loc == 'nodes': cellNPresent = C.isNamePresent(aR, 'cellN')
    else: cellNPresent = C.isNamePresent(aR, 'centers:cellN')
    if cellNPresent == -1:
        if storage == 'direct': return aR
        else: return aD

    aD = addCellN__(aD,loc='nodes')
    zonesRcv = Internal.getZones(aR); nzonesRcv = len(zonesRcv)
    zonesDnr = Internal.getZones(aD); nzonesDnr = len(zonesDnr)

    #---------------------------------------------------------------------------
    # Extraction des pts IBC a corriger + pts parois + interpoles
    #---------------------------------------------------------------------------
    if hook is not None:
        allHooks = hook[:]
        allHooksL = allHooks[:]
    else: allHooks = None

    zonesDnrL = zonesDnr[:]
    nzonesDnrL = nzonesDnr
    nozr = -1
    #-------------------------------------------
    # 1. Get the list of IBC pts
    #-------------------------------------------
    res = IBM.getAllIBMPoints(zonesRcv,loc=locR, hi=hi, he=he)
    correctedPts=res[0]; wallPts=res[1];  interpPts=res[2]

    #-------------------------------------------
    # 2. Interpolation of IBC interp pts
    #-------------------------------------------
    for z in zonesRcv:
        nozr += 1
        if loc == 'nodes': cellNPresent = C.isNamePresent(z,'cellN')
        else: cellNPresent = C.isNamePresent(z, 'centers:cellN')
        if cellNPresent != -1:
            if sameName == 1:
                if hook is not None: allHooks = allHooksL[:]
                zonesDnr = zonesDnrL[:]
                nzonesDnr = nzonesDnrL
                cL = 0; found = 0
                for zo in zonesDnr:
                    if zo[0] == z[0]: found = 1; break
                    cL += 1
                if found == 1:
                    if hook is not None: allHooks.pop(cL)
                    zonesDnr.pop(cL); nzonesDnr = nzonesDnr-1


            _setIBCDataForZone__(z, zonesDnr, correctedPts[nozr], wallPts[nozr], interpPts[nozr], \
                                    loc=locR, order=order, \
                                    penalty=penalty,nature=nature,method=method,storage=storage,hook=hook,dim=dim)

    # fin parcours des zones receveuses
    if storage == 'direct': return aR
    else:
        ztD = Internal.getZones(tD)
        zaD = Internal.getZones(aD)
        for i in xrange(len(ztD)): # enleve le cellN si interpData l'a ajoute
            ret = C.isNamePresent(ztD[i], 'cellN')
            if ret == -1: C._rmVars(zaD[i],['cellN'])
        return aD

def _setIBCDataForZone__(z, zonesDnr, correctedPts, wallPts, interpPts, loc='nodes', \
                             order=2, penalty=0, nature=0, method='lagrangian', storage='direct', hook=None, dim=3, bcType =0):

    arraysD = C.getFields(Internal.__GridCoordinates__, zonesDnr)
    cellND  = C.getField('cellN', zonesDnr); arraysD = Converter.addVars([arraysD,cellND])

    model = "Euler"
    a = Internal.getNodeFromName2(zonesDnr[0], 'model')
    if a is not None: model = Internal.getValue(a)
    if model != "Euler": bcType = 3

    nzonesDnr = len(arraysD)
    #-------------------------------------------
    # 3. Interpolation of IBC interp pts
    #-------------------------------------------
    # resInterp = [rcvInd1D,donorInd1D,donorType,coefs,extrap,orphan, EXdirs]
    resInterp = Connector.setInterpData__(interpPts, arraysD, order=order, penalty=penalty, nature=nature, method=method, hook=hook, dim=dim)

    if resInterp is not None:
        # Bilan
        nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
        nbinterpolated0 = interpPts[1].shape[1]
        if len(resInterp[4]) > 0:
            for r in resInterp[4]: nbextrapolated += r.shape[0]
        nborphan = resInterp[5].size
        nbinterpolated = nbinterpolated0-nbextrapolated-nborphan

        print 'Zone %s: interpolated=%d; extrapolated=%d; orphan=%d'%(z[0],nbinterpolated,nbextrapolated,nborphan)
        if  nborphan>0: print 'Warning: zone %s has %d orphan points !'%(z[0], nborphan)

        #-------------------------------------------------------------------------
        # 4. Update the indices of interpolated pts -> indices of corrected points
        #-------------------------------------------------------------------------
        posindcell = KCore.isNamePresent(interpPts, 'indcell')
        indcells = interpPts[1][posindcell,:]

        #-------------------------------------------------------------------------
        # 5. Sort IBC coordinates wrt donor zones
        #-------------------------------------------------------------------------
        allCorrectedPts=[[]]*nzonesDnr; allWallPts=[[]]*nzonesDnr; allMirrorPts=[[]]*nzonesDnr
        xPC0 = correctedPts[1][0,:]; yPC0 = correctedPts[1][1,:]; zPC0 = correctedPts[1][2,:]
        xPW0 = wallPts[1][0,:]; yPW0 = wallPts[1][1,:]; zPW0 = wallPts[1][2,:]
        xPI0 = interpPts[1][0,:]; yPI0 = interpPts[1][1,:]; zPI0 = interpPts[1][2,:]

        #-------------------------------------------------------------------------
        # 6. Set the original indices of corrected pts in receptor zones
        #-------------------------------------------------------------------------
        # Orphan
        if nborphan>0:
            for noi in xrange(nborphan):
                noind = resInterp[5][noi]
                resInterp[5][noi] = indcells[noind]

        # Interpolated/Extrapolated
        for noz in xrange(nzonesDnr):
            ninterploc = resInterp[0][noz].size
            if ninterploc>0: # domaine d'interpolation
                correctedPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPC = correctedPtsL[1][0,:]; yPC = correctedPtsL[1][1,:]; zPC = correctedPtsL[1][2,:]
                wallPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPW = wallPtsL[1][0,:]; yPW = wallPtsL[1][1,:]; zPW = wallPtsL[1][2,:]
                mirrorPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPI = mirrorPtsL[1][0,:]; yPI = mirrorPtsL[1][1,:]; zPI = mirrorPtsL[1][2,:]

                for noi in xrange(ninterploc):
                    index = resInterp[0][noz][noi]
                    # Indices of receptor pts
                    resInterp[0][noz][noi] = indcells[index]
                    # coordinates of receptor/wall/mirror pts
                    xPC[noi] = xPC0[index]; yPC[noi] = yPC0[index]; zPC[noi] = zPC0[index]
                    xPW[noi] = xPW0[index]; yPW[noi] = yPW0[index]; zPW[noi] = zPW0[index]
                    xPI[noi] = xPI0[index]; yPI[noi] = yPI0[index]; zPI[noi] = zPI0[index]

                allCorrectedPts[noz] = correctedPtsL
                allWallPts[noz] = wallPtsL
                allMirrorPts[noz] = mirrorPtsL

        if len(resInterp[4])>0: # Sort extrapolated points wrt donor zones
            for noz in xrange(nzonesDnr):
                nextraploc = resInterp[4][noz].size
                if nextraploc>0: # Extrapoles
                    for noi in xrange(nextraploc):
                        index = resInterp[4][noz][noi]
                        resInterp[4][noz][noi] = indcells[index]

        #----------------------------------
        # 7. Stockage dans l' arbre
        # direct: on stocke dans aR
        # inverse: on stocke dans aD
        #----------------------------------   
        PL = numpy.array([],numpy.int32)
        PLD = numpy.array([],numpy.int32)
        EXTRAP = numpy.array([],numpy.int32)
        INTERPTYPE = numpy.array([],numpy.int32)
        CFS = numpy.array([],numpy.float64)
        VOL = numpy.array([],numpy.float64)
        ORPHAN = resInterp[5]
        for noz in xrange(nzonesDnr):
            # ajout des donnees d'interpolation
            ninterploc = resInterp[0][noz].size
            if ninterploc>0:# est un domaine donneur
                if resInterp[4][noz].size > 0:  EXTRAP = resInterp[4][noz]

                if storage == 'direct':
                    createInterpRegion__(z, zonesDnr[noz][0],resInterp[0][noz],resInterp[1][noz],resInterp[3][noz],\
                                             resInterp[2][noz], VOL, EXTRAP, ORPHAN, \
                                             tag='Receiver', loc=loc,itype='ibc')
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(z, zonesDnr[noz][0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType)
                else: # inverse
                    createInterpRegion__(zonesDnr[noz],z[0],resInterp[1][noz],resInterp[0][noz],resInterp[3][noz],\
                                         resInterp[2][noz], VOL, EXTRAP, ORPHAN, \
                                         tag='Donor', loc=loc,itype='ibc')
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(zonesDnr[noz], z[0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType)
                    
            elif nborphan==nbinterpolated0: # Only orphan pts: orphan pt list is stored in any candidate donor zone
                if storage == 'direct':
                    createInterpRegion__(z, zonesDnr[noz][0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                         tag='Receiver', loc=loc,itype='ibc')
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(z, zonesDnr[noz][0], correctedPts, wallPts, interpPts, bcType)
                else: # inverse
                    createInterpRegion__(zonesDnr[noz],z[0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                         tag='Donor', loc=loc,itype='ibc')
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(zonesDnr[noz], z[0], correctedPts, wallPts, interpPts, bcType)    
      
    return None

def _addIBCCoords__(z, zname,correctedPts, wallPts, interpolatedPts, bcType):
    nameSubRegion = 'IBCD_'+zname
    zsr = Internal.getNodesFromName1(z, nameSubRegion)
    coordsPC = Converter.extractVars(correctedPts,['CoordinateX','CoordinateY','CoordinateZ'])
    coordsPW = Converter.extractVars(wallPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    coordsPI = Converter.extractVars(interpolatedPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    zsr[0][2].append(['CoordinateX_PC',coordsPC[1][0,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateY_PC',coordsPC[1][1,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateZ_PC',coordsPC[1][2,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateX_PI',coordsPI[1][0,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateY_PI',coordsPI[1][1,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateZ_PI',coordsPI[1][2,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateX_PW',coordsPW[1][0,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateY_PW',coordsPW[1][1,:], [], 'DataArray_t'])
    zsr[0][2].append(['CoordinateZ_PW',coordsPW[1][2,:], [], 'DataArray_t'])

    # Creation des numpy d extraction
    nIBC    = coordsPC[1].shape[1]
    #print   nIBC
    pressNP = numpy.zeros((nIBC),numpy.float64)
    densNP  = numpy.zeros((nIBC),numpy.float64)

    zsr[0][2].append(['Pressure', pressNP, [], 'DataArray_t'])
    zsr[0][2].append(['Density' , densNP , [], 'DataArray_t'])
    utauNP  = numpy.zeros((nIBC),numpy.float64)
    yplusNP = numpy.zeros((nIBC),numpy.float64)
    zsr[0][2].append(['utau' , utauNP , [], 'DataArray_t'])
    zsr[0][2].append(['yplus', yplusNP, [], 'DataArray_t'])
    return None

#==============================================================================
# Uses a Dictionnary of intersections and calls setInterpData when appropriate
#==============================================================================
'''     ################### APPROVAL REQUIRED ####################
def setInterpDataDict(tR, tD, double_wall=0, order=2, penalty=1, nature=0,
                  method='lagrangian', loc='nodes', storage='direct',
                  hook=None,
                  topTreeRcv=None, topTreeDnr=None, sameName=1, dim=3,
                  intersectionsDict=None):
    if not intersectionsDict:
        intersectionsDict = getIntersectingDomains(tR, tD, mode='hybrid')
    for z in Internal.getZones(tR):
        nobr,nozr = C.getNobNozOfZone(z, tR)
        dnrZones = []; nobOfDnrBases = []; nobOfDnrZones=[];
        for zd in Internal.getZones(tD):
            if zd[0] in intersectionsDict[z[0]]:
                dnrZones.append(zd)
                nobd,nozd = C.getNobNozOfZone(zd, tD)
                nobOfDnrBases.append(nobd)
                nobOfDnrZones.append(nozd)

        if storage=='inverse':
            dnrZones = setInterpData(z,dnrZones, double_wall, order, penalty, nature,
                          method, loc, storage, hook,
                          topTreeRcv, topTreeDnr, sameName, dim)
            for i in xrange(0,len(dnrZones)):
                tD[2][nobOfDnrBases[i]][2][nobOfDnrZones[i]] = dnrZones[i]
        else:
            z = setInterpData(z,dnrZones, double_wall, order, penalty, nature,
                          method, loc, storage, hook,
                          topTreeRcv, topTreeDnr, sameName, dim)
            tR[2][nobr][2][nozr] = z
    if storage=='inverse':
        return tD
    else:
        return tR
'''

# Adapt aD pour RANS/LES
def _adaptForRANSLES__(tR, aD):
    zrdict = {}
    zones = Internal.getNodesFromType2(tR, 'Zone_t')
    for z in zones: zrdict[z[0]] = z

    zonesD = Internal.getNodesFromType2(aD, 'Zone_t')
    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:
            zrcvname = Internal.getValue(s)
            try:
                zr = zrdict[zrcvname]
                model_z1 = 'Euler'; model_z2 = 'Euler'
                a = Internal.getNodeFromName2(zr, 'GoverningEquations')
                if a is not None: model_z1 = Internal.getValue(a)
                a = Internal.getNodeFromName2(zd, 'GoverningEquations')
                if a is not None: model_z2 = Internal.getValue(a)

                if ((model_z2=='NSTurbulent' or model_z1=='NSTurbulent') and  model_z1 != model_z2):
                   datap = numpy.ones(1, numpy.int32)
                   Internal.createUniqueChild(s, 'RANSLES', 'DataArray_t', datap)
            except: pass
    return None

#==============================================================================
# Calcul des donnees d'interpolation et stockage dans l'arbre CGNS/Python
# ----------------------------------------------------------------------------
# -------------------
# Si la variable cellN n'existe pas, toute la zone receveuse est interpolee
# Si la variable cellN existe, seuls les pts de cellN 2 sont interpoles (noeuds ou centres)
# Si loc='faces' et cellN existe pour des centres, interpolations aux pts EX effectues
# -------------------
# IN: aR: arbre, zone, bases des receveurs
# IN: aD: arbre, zone, bases des donneurs
# IN: order: ordre des interpolations (2, 3, 5)
# IN: method='lagrange', 'leastsquares','conservative'
# IN: penalty=1: penalise une cellule donneuse en terme de volume si elle est au bord
# IN: nature=0: aucun sommet de la cellule d'interpolation ne doit etre avec un cellN=0
#     nature=1: toutes les sommets de la cellule d'interpolation doivent etre de cellN=1
# IN: hook: hook sur l'adt (pas reconstruit dans setInterpData), l'ordre doit suivre celui de zonesD
# IN: double_wall=1: activation de la technique double wall
# IN: storage: type de stockage (direct: sur le bloc interpole, inverse: sur le bloc d'interpolation)
# IN: loc='nodes','cells','faces': interpolation appliquee pour les receveurs (localises en noeuds/centres/faces)
# IN: topTreeR: top tree des receveurs, sert pour extraire les FamilyBC du double wall
# IN: topTreeD: top tree des donneurs, sert pour extraire les FamilyBC du double wall
# IN: sameName: si 1, n'interpole pas a partir des zones donneuses portant le meme nom que les zones receveuses
# IN: itype='both','chimera','abutting': calcule toutes les donnees de transferts, seules les chimere, seules les match/nearmatch
# OUT: stockage direct: retourne tR avec les donnees d'interpolation
# OUT: stockage indirect: retourne tD avec les donnees d'interpolation
# RMQ : method='conservative' -> tout le domaine receveur est pour l instant considere a interpoler (maquette)
#==============================================================================
def setInterpData(tR, tD, double_wall=0, order=2, penalty=1, nature=0,
                  method='lagrangian', loc='nodes', storage='direct',
                  hook=None,
                  topTreeRcv=None, topTreeDnr=None, sameName=1, dim=3, itype='both'):

    locR = loc
    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)

    # Recherche pour les pts coincidents (base sur les GridConnectivity)
    if itype != 'chimera':
        if storage == 'direct': aR = setInterpDataForGhostCells__(aR,aD,storage,loc)
        else: 
            aD = setInterpDataForGhostCells__(aR,aD,storage,loc)

            # Determination du model pour RANS/LES : Stef2Ivan: pourquoi pas aussi en stockage direct
            _adaptForRANSLES__(tR, aD)

    # Si pas de cellN receveur, on retourne
    if loc == 'nodes': cellNPresent = C.isNamePresent(aR, 'cellN')
    elif loc=='centers': cellNPresent = C.isNamePresent(aR, 'centers:cellN')
    else: 
        raise ValueError("setInterpData: invalid loc provided.")
    if cellNPresent == -1 or itype == 'abutting':
        if storage == 'direct': return aR
        else: return aD

    locCellND = 'nodes'
    aD = addCellN__(aD, loc=locCellND)
    
    if method == 'conservative' and itype != 'abutting':  
        if loc != 'centers':
            raise ValueError("setInterpData: conservative type is available only for loc='centers'.")
        else: return setInterpDataConservative__(aR, aD, storage=storage)

    zonesRcv = Internal.getZones(aR); nzonesRcv = len(zonesRcv)
    zonesDnr = Internal.getZones(aD); nzonesDnr = len(zonesDnr)

    #---------------------------------------------------------------------------
    # CAS DOUBLE WALL
    # Extraction des parois de projection issues des zones donneuses
    # Extraction des premiers centres ou premiers noeuds (selon locR) des zones receveuses
    #---------------------------------------------------------------------------
    donorSurfs = []; interpWallPts = []
    noWallsInDnr = 1 # parametre qui vaut 1 si pas de parois trouvees dans les zones donneuses, 0 si au moins une
    noWallsInRcv = 1 # idem mais pour les receveurs
    if double_wall == 1:
        import DoubleWall
        try: import Geom.PyTree as D
        except: raise ImportError("setInterpData+double wall requires Geom.PyTree module.")

        ttreeR = []
        if topTreeRcv is not None: ttreeR = topTreeRcv
        else:
            if Internal.isTopTree(aR): ttreeR = aR
            else:
                if Internal.getBases(aR) != []: ttreeR = aR# c est une base on peut recuperer les familles
                else: print 'Warning: setInterpData+double wall: receptor zones may require a top tree.'
        ttreeD = []
        if topTreeDnr is not None: ttreeD = topTreeDnr
        else:
            if Internal.isTopTree(aD): ttreeD = aD
            else:
                if Internal.getBases(aD) != []: ttreeD = aD # c'est une base on peut recuperer les familles
                else: print 'Warning: setInterpData+double wall: donors zones may require a top tree.'

        # Zones donneuses : on recupere les surfaces de projection double wall
        for zd in zonesDnr:
            walls = C.extractBCOfType(zd,'BCWall',topTree=ttreeD)
            if walls != []:
                noWallsInDnr = 0
                walls = D.getCurvatureHeight(walls)
                walls = C.convertArray2Tetra(walls, split="withBarycenters")
                walls = C.getAllFields(walls,loc='nodes')
            donorSurfs.append(walls)

        # Zones receveuses : on determine les 1ers points paroi (centres ou noeuds)
        # recup des familles de type paroi dans les zones receveuses
        famwallsR = C.getFamilyBCNamesOfType(ttreeR, 'BCWall')
        famwallsR += C.getFamilyBCNamesOfType(ttreeR, 'BCWallViscous')
        famwallsR += C.getFamilyBCNamesOfType(ttreeR, 'BCWallInviscid')

        # interpWallPts : selon la loc, on recupere les premiers centres ou noeuds
        for zr in zonesRcv:
            wallRanges = DoubleWall.getBCWallRanges__(zr,famwallsR)
            if wallRanges != []:
                noWallsInRcv = 0
                if locR == 'nodes': interpWallPts.append(DoubleWall.getFirstPointsInfo__(zr, wallRanges,loc='nodes'))
                else: interpWallPts.append(DoubleWall.getFirstPointsInfo__(zr, wallRanges,loc='centers'))
            else: interpWallPts.append([])

    if noWallsInDnr == 1 or noWallsInRcv == 1: double_wall = 0 # on desactive le double wall
    arraysD = C.getFields(Internal.__GridCoordinates__, zonesDnr)
    cellND = C.getField('cellN', zonesDnr)
    arraysD = Converter.addVars([arraysD,cellND])
    cellND = []
    #---------------------------------------------------------------------------
    # 1. Extraction des points a interpoler
    #    interpPts : un array si pas de double wall
    #              : une liste d arrays avec les coordonnees des pts a interpoler modifiees
    # 2. Calcul des coefs et donneurs
    #---------------------------------------------------------------------------
    if hook is not None:
        allHooks = hook[:]
        allHooksL = allHooks[:]
    else: allHooks = None

    arraysDL = arraysD[:]
    donorSurfsL = donorSurfs[:]
    zonesDnrL = zonesDnr[:]
    nzonesDnrL = nzonesDnr
    nozr = -1

    for z in zonesRcv:
        nozr += 1
        if loc == 'nodes': cellNPresent = C.isNamePresent(z, 'cellN')
        else: cellNPresent = C.isNamePresent(z, 'centers:cellN')
        if cellNPresent != -1:
            if sameName == 1:
                arraysD = arraysDL[:]
                if hook is not None: allHooks = allHooksL[:]
                donorSurfs = donorSurfsL[:]
                zonesDnr = zonesDnrL[:]
                nzonesDnr = nzonesDnrL
                cL = 0; found = 0
                for zo in zonesDnr:
                    if zo[0] == z[0]: found = 1; break;
                    cL += 1
                if found == 1:
                    arraysD.pop(cL)
                    if hook is not None: allHooks.pop(cL)
                    zonesDnr.pop(cL); nzonesDnr = nzonesDnr-1
                    if double_wall == 1: donorSurfs.pop(cL)

            #-------------------------------------------
            # Etape 1: on recupere les pts a interpoler
            #-------------------------------------------
            interpPts = [];
            isdw = 0 # si double wall effectivement active, interpPts est une liste d'arrays
            if locR == 'nodes':
                an = C.getFields(Internal.__GridCoordinates__, z)[0]
                cellN = C.getField('cellN',z)[0]
                an = Converter.addVars([an,cellN])
                if double_wall == 1 and interpWallPts[nozr] != []: # dw: liste d arrays
                    isdw = 1
                    for nozd in xrange(nzonesDnr):
                        an2 = Connector.changeWall__(an, interpWallPts[nozr], donorSurfs[nozd])
                        interpPts.append(Connector.getInterpolatedPoints__(an2))
                else: # pas de dw: un seul array
                    interpPts = Connector.getInterpolatedPoints__(an)

            elif locR == 'centers':
                an = C.getFields(Internal.__GridCoordinates__, z)[0]
                ac = Converter.node2Center(an)
                cellN = C.getField('centers:cellN',z)[0]
                ac = Converter.addVars([ac,cellN])
                if double_wall == 1 and interpWallPts[nozr] != []:# dw : liste d arrays
                    isdw = 1
                    for nozd in xrange(nzonesDnr):
                        ac2 = Connector.changeWall__(ac, interpWallPts[nozr], donorSurfs[nozd])
                        interpPts.append(Connector.getInterpolatedPoints__(ac2))
                else:  # pas de dw : un seul array
                    interpPts = Connector.getInterpolatedPoints__(ac)

            #---------------------------------------------
            # Etape 2 : calcul des donnees d'interpolation
            #---------------------------------------------
            # resInterp = [rcvInd1D,donorInd1D,donorType,coefs,extrap,orphan, EXdirs]
            resInterp = Connector.setInterpData__(interpPts, arraysD, order=order, penalty=penalty, nature=nature, method=method, hook=allHooks, dim=dim)
            if resInterp is not None:
                # Bilan
                nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
                if double_wall==0: nbinterpolated = interpPts[1].shape[1]
                else: nbinterpolated = interpPts[0][1].shape[1]
                if len(resInterp[4])>0:
                    for r in resInterp[4]: nbextrapolated += r.shape[0]
                nborphan = resInterp[5].size
                nbinterpolated = nbinterpolated-nbextrapolated-nborphan

                print 'Zone %s: interpolated=%d ; extrapolated=%d ; orphan=%d'%(z[0], nbinterpolated, nbextrapolated, nborphan)
                if  nborphan>0: print 'Warning: zone %s has %d orphan points !'%(z[0], nborphan)
                # on remet a une seule zone, attention si x,y,z sont necessaires ensuite
                # les coordonnees peuvent etre fausses a cause du double walls
                indcells=[]
                if loc == 'faces': vari = 'indcell1'
                else: vari = 'indcell'
                if isdw == 1:
                    posindcell = KCore.isNamePresent(interpPts[0],vari)
                    if posindcell != -1: indcells = interpPts[0][1][posindcell,:]
                else:
                    posindcell = KCore.isNamePresent(interpPts,vari)
                    if posindcell != -1: indcells = interpPts[1][posindcell,:]

                # on recupere les bons indices de pts interpoles (indcell ds interpPts), de EXdir
                # pour ceux obtenus en sortie de setInterpData
                # Orphelins
                if nborphan > 0:
                    for noi in xrange(nborphan):
                        noind = resInterp[5][noi]
                        resInterp[5][noi] = indcells[noind]
                # Interpoles/Extrapoles
                for noz in xrange(nzonesDnr):
                    ninterploc = resInterp[0][noz].size
                    if ninterploc > 0:# domaine d'interpolation
                        # Indices des receveurs
                        for noi in xrange(ninterploc):
                            index = resInterp[0][noz][noi]
                            resInterp[0][noz][noi] = indcells[index]
                if len(resInterp[4])>0: # pts extrapoles par zone donneuse
                    for noz in xrange(nzonesDnr):
                        nextraploc = resInterp[4][noz].size
                        if nextraploc > 0:# Extrapoles
                            for noi in xrange(nextraploc):
                                index = resInterp[4][noz][noi]
                                resInterp[4][noz][noi] = indcells[index]

                #----------------------------------
                # Etape 3: Stockage dans l'arbre
                # direct: on stocke dans aR
                # inverse: on stocke dans aD
                #----------------------------------
                for noz in xrange(nzonesDnr):
                    # ajout des donnees d interpolation
                    ninterploc = resInterp[0][noz].size
                    if ninterploc>0: # domaine d'interpolation
                        if storage == 'direct':
                            if resInterp[4][noz].size == 0: extrapPts = numpy.array([],numpy.int32)
                            else: extrapPts = resInterp[4][noz]
                            if resInterp[6] == []: EXdir = numpy.array([],numpy.int32)
                            else: EXdir = resInterp[6][noz]
                            createInterpRegion__(z, zonesDnr[noz][0], resInterp[0][noz], resInterp[1][noz], resInterp[3][noz], \
                                                     resInterp[2][noz], numpy.array([],numpy.float64), \
                                                     extrapPts, resInterp[5], tag='Receiver', loc=locR, EXDir=EXdir)

                        else: # inverse
                            if resInterp[4][noz].size == 0: extrapPts = numpy.array([],numpy.int32)
                            else: extrapPts = resInterp[4][noz]
                            if resInterp[6] == []: EXdir = numpy.array([],numpy.int32)
                            else: EXdir = resInterp[6][noz]
                            createInterpRegion__(zonesDnr[noz], z[0], resInterp[1][noz], resInterp[0][noz], resInterp[3][noz], \
                                                     resInterp[2][noz], numpy.array([],numpy.float64), \
                                                     extrapPts, resInterp[5], tag='Donor', loc=locR, EXDir=EXdir)
    if storage != 'direct':
        _adaptForRANSLES__(tR, aD)

    # fin parcours des zones receveuses
    if storage == 'direct': return aR
    else:
        ztD = Internal.getZones(tD)
        zaD = Internal.getZones(aD)
        for i in xrange(len(ztD)): # enleve le cellN is interpData l'a ajoute
            ret = C.isNamePresent(ztD[i], 'cellN')
            if ret == -1: C._rmVars(zaD[i],['cellN'])
        return aD

#-------------------------------------------------------------------------
# setInterpDataConservative__
#-------------------------------------------------------------------------
def setInterpDataConservative__(tR, tD, storage='direct'):
    try: import Post.PyTree as P; import Generator.PyTree as G
    except: raise ImportError("setInterpDataConservative__: requires Post module.")
    locR = 'centers'
    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)

    tRBB = G.BB(aR,method='AABB')
    tDBB = G.BB(aD,method='AABB')
    intersectionDict = getIntersectingDomains(aR, aD, method='AABB',taabb=tRBB, taabb2=tDBB)
    # hook sur les zones donneuses
    allDnrCells={}; allDnrZones = {}; allIndicesDnrOrig={}
    for zd in Internal.getZones(aD):
        hookD = C.createHook(zd,'elementCenters')
        zdnr = P.selectCells(zd,'({centers:cellN}>0.)*({centers:cellN}<2.)>0.')
        zdnr = C.convertArray2NGon(zdnr); zdnr = G.close(zdnr)
        allDnrCells[zd[0]]=zdnr; allDnrZones[zd[0]]=zd
        indicesDnrOrig = C.identifyElements(hookD,zdnr)
        allIndicesDnrOrig[zd[0]]=indicesDnrOrig 
        C.freeHook(hookD)

    for zr in Internal.getZones(aR):
        #1. recuperation des points a interpoler
        interpPts = P.selectCells(zr,'{centers:cellN}>1.')
        if interpPts != []:
            interpPts=C.convertArray2NGon(interpPts); interpPts = G.close(interpPts)
            #indices correspondant
            hookR = C.createHook(zr,'elementCenters')
            indicesRcvOrig = C.identifyElements(hookR,interpPts)
            C.freeHook(hookR)
        
            #2. calcul des donnees d interpolation
            # SP : pour optimiser on pourra prendre le dictionnaire d intersection des pts interpoles uniquement
            arraysD=[];  donorZoneNames=[]; listOfIndicesDnrOrig=[]
            for zdnrname in intersectionDict[zr[0]]:
                zdnr = allDnrCells[zdnrname]
                if zdnr != []:
                    ad = C.getFields(Internal.__GridCoordinates__,zdnr)[0]
                    arraysD.append(ad)
                    listOfIndicesDnrOrig.append(allIndicesDnrOrig[zdnrname])
                    donorZoneNames.append(zdnrname)
            nzonesDnr = len(arraysD)
            interpPtsA = C.getFields(Internal.__GridCoordinates__,interpPts)[0]

            resInterp = Connector.connector.setInterpDataCons(interpPtsA, arraysD, indicesRcvOrig, listOfIndicesDnrOrig)
            if resInterp is not None:
                # Bilan
                nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
                for indicesR in resInterp[0]: nbinterpolated += indicesR.shape[0]
                nborphan = resInterp[4].size
                print 'Zone %s: interpolated=%d ; orphan=%d'%(zr[0], nbinterpolated, nborphan)
                if  nborphan>0: print 'Warning: zone %s has %d orphan points !'%(zr[0], nborphan)

                # Orphelins
                if nborphan > 0:
                    for noi in xrange(nborphan):
                        index = resInterp[4][noi]
                        resInterp[4][noi]=indicesRcvOrig[index]
                    orphanPts = resInterp[4]
                else: orphanPts = numpy.array([],numpy.int32)

                extrapPts = numpy.array([],numpy.int32)
                EXdir = numpy.array([],numpy.int32)                   

                # 3. stockage dans l arbre                            
                # on recupere les bons indices de pts interpoles (indcell ds interpPts)
                for noz in xrange(nzonesDnr):
                    dnrname=donorZoneNames[noz]
                    # ajout des donnees d interpolation
                    ninterploc = resInterp[0][noz].size
                    if ninterploc>0: # domaine d'interpolation
                        if storage == 'direct':
                            createInterpRegion__(zr, dnrname, resInterp[0][noz], resInterp[1][noz], resInterp[3][noz], \
                                                 resInterp[2][noz], numpy.array([],numpy.float64), \
                                                 extrapPts, orphanPts, tag='Receiver', loc=locR, EXDir=EXdir)

                        else: # inverse
                            extrapPts = numpy.array([],numpy.int32)
                            EXdir = numpy.array([],numpy.int32)
                            zd = allDnrZones[dnrname]
                            createInterpRegion__(zd, zr[0], resInterp[1][noz], resInterp[0][noz], resInterp[3][noz], \
                                                     resInterp[2][noz], numpy.array([],numpy.float64), \
                                                     extrapPts, orphanPts, tag='Donor', loc=locR, EXDir=EXdir)


    if storage != 'direct': _adaptForRANSLES__(tR, aD)

    # fin parcours des zones receveuses
    if storage == 'direct': return aR
    else:
        ztD = Internal.getZones(tD)
        zaD = Internal.getZones(aD)
        for i in xrange(len(ztD)): # enleve le cellN is interpData l'a ajoute
            ret = C.isNamePresent(ztD[i], 'cellN')
            if ret == -1: C._rmVars(zaD[i],['cellN'])
        return aD

#===============================================================================
# Uses the topological information to set the interpolation data in the tree
# for 1-to-1 grid connectivity. Receptor points are ghost points (centers or nodes)
# IN: tR: pyTree with ghost cells already created
#          BCMatch ranges correspond to the original pytree without ghost cells
#          must contain a Rind_t node for each zone
# IN: tD: pyTree of donors (consistent with setInterpTransfers)
# IN: d : number of ghost cells that have been created
# IN: loc='nodes' or 'centers'
# IN: 'storage'='direct' or 'inverse'
# OUT: t: with interpolation data stored in ZoneSubRegion_t nodes of name 'ID_*'
#===============================================================================
def setInterpDataForGhostCells__(tR, tD, storage='direct', loc='nodes'):
    try: import Converter.GhostCells as GhostCells
    except: raise ImportError("setInterpDataForGhostCells__ requires Converter.GhostCells module.")
    # empty numpy arrays for zonesubregion nodes
    indicesExtrap = numpy.array([],numpy.int32)
    indicesOrphan = numpy.array([],numpy.int32)
    vols =  numpy.array([],numpy.float64)
    EXdir = numpy.array([],numpy.int32)

    if loc == 'nodes': locR = 0; locS = 'Vertex'
    else: locR = 1; locS = 'CellCenter'

    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)
    for zp in Internal.getZones(aR):
        zname = zp[0]
        zoneDimR = Internal.getZoneDim(zp)
        if zoneDimR[0] == 'Unstructured':
            print 'Warning: setInterpDataForGC not yet implemented for unstructured zones.'
        else: # Structured
            dimPb = zoneDimR[4]
            rindnode = Internal.getNodeFromType1(zp, 'Rind_t')
            if rindnode is not None: # rind indices exist : ghost cell data to be computed
                rindnode = rindnode[1]
                rindimin = rindnode[0]; rindimax = rindnode[1]
                rindjmin = 0; rindjmax = 0; rindkmin = 0; rindkmax = 0
                if dimPb > 1:
                    rindjmin = rindnode[2]; rindjmax = rindnode[3]
                    if dimPb == 3: rindkmin = rindnode[4]; rindkmax = rindnode[5]

                rindrcv = [rindimin,rindimax,rindjmin,rindjmax,rindkmin,rindkmax]
                imr = zoneDimR[1]; jmr = zoneDimR[2]; kmr = zoneDimR[3]
                if locR == 1: imr=imr-1; jmr=jmr-1; kmr=kmr-1
                listofjoins = Internal.getNodesFromType2(zp, 'GridConnectivity1to1_t')
                for join in listofjoins:
                    # receptor window ranges
                    prange = Internal.getNodeFromName1(join,'PointRange')[1]
                    # donor window ranges
                    prangedonor = Internal.getNodeFromName1(join,'PointRangeDonor')[1]

                    # trirac
                    if dimPb == 3: trirac = [1,2,3]
                    elif dimPb == 2: trirac = [1,2]
                    else: trirac = [1]
                    transfo = Internal.getNodeFromName1(join, 'Transform')
                    if transfo is not None:
                        trirac[0] = transfo[1][0]
                        if dimPb != 1: trirac[1] = transfo[1][1]
                        if dimPb == 3: trirac[2] = transfo[1][2]
                    Periodic = Internal.getNodeFromType2(join,'Periodic_t')
                    RotationAngle=None; RotationCenter=None
                    if Periodic is not None:
                        RotationAngle = Internal.getNodeFromName1(Periodic,'RotationAngle')
                        RotationCenter = Internal.getNodeFromName1(Periodic,'RotationCenter')
                        if RotationAngle is not None:
                           RotationAngle[1][0]=-RotationAngle[1][0]
                           RotationAngle[1][1]=-RotationAngle[1][1]
                           RotationAngle[1][2]=-RotationAngle[1][2]

                    # donor zone name
                    zdonorname = Internal.getValue(join)
                    zdonor = Internal.getZones(aR)
                    zdonor = Internal.getNodesFromName(zdonor, zdonorname)
                    zdonor = zdonor[0] # Warning: if a zone name is not unique, the first zone is kept
                    # check if the donor zone is defined in aD
                    zdonorp = Internal.getNodesFromName2(aD, zdonorname)
                    if zdonorp == []:
                        raise ValueError("setInterpDataForGhostCells: donor zone not found in donor pyTree.")
                    zdonorp = Internal.getNodesFromType2(zdonorp, 'Zone_t')
                    if zdonorp == []:
                        raise ValueError("setInterpDataForGhostCells: donor zone not found in donor pyTree.")
                    if len(zdonorp)  > 1 :
                        print 'Warning: setInterpDataForGhostCells: zone name %s defined several times in donor pyTree.'%zdonorname
                    zdonorp = zdonorp[0]
                    # donor zone dimensions
                    zoneDimD = Internal.getZoneDim(zdonor)
                    imd = zoneDimD[1]; jmd = zoneDimD[2]; kmd = zoneDimD[3]
                    if locR == 1: imd=imd-1; jmd=jmd-1; kmd=kmd-1

                    # rind for donor
                    rindnode = Internal.getNodeFromType1(zdonor, 'Rind_t')
                    if rindnode is not None: # rind indices exist: ghost cell data to be computed
                        rindnode = rindnode[1]
                        rindimin = rindnode[0]; rindimax = rindnode[1]
                        rindjmin = 0; rindjmax = 0; rindkmin = 0; rindkmax = 0
                        if dimPb > 1:
                            rindjmin = rindnode[2]; rindjmax = rindnode[3]
                            if dimPb == 3: rindkmin = rindnode[4]; rindkmax = rindnode[5]
                        #
                        rinddnr = [rindimin,rindimax,rindjmin,rindjmax,rindkmin,rindkmax]
                        # get directions of receptor and donor borders
                        dirR = GhostCells.getDirBorderStruct__(prange,dimPb)
                        dirD = GhostCells.getDirBorderStruct__(prangedonor,dimPb)
                        # get list of border nodes and direction of border
                        if dimPb != 1:
                            # arrayOfIndicesR : indices globaux des 1ers points pres de la frontiere
                            # definie par le prange dimensionnes dans le maillage ghost cells
                            shift = 0 # les indices sont ceux des pts en frontiere max GC
                            arrayOfIndicesR = GhostCells.getBorderIndicesStruct__(prange,zoneDimR,dirR,0,locS,dimPb,shift)

                            # listOfIndicesD : indices globaux des 1ers pts donneurs associes a ceux definis par
                            # arrayOfIndicesR
                            shift = 0
                            if dirD  ==-1: shift += rinddnr[0]
                            elif dirD== 1: shift += rinddnr[1]
                            elif dirD==-2: shift += rinddnr[2]
                            elif dirD== 2: shift += rinddnr[3]
                            elif dirD==-3: shift += rinddnr[4]
                            else: shift+= rinddnr[5]
                            if dirR  ==-1: shift += rindimin
                            elif dirR== 1: shift += rindimax
                            elif dirR==-2: shift += rindjmin
                            elif dirR== 2: shift += rindjmax
                            elif dirR==-3: shift += rindkmin
                            else: shift+= rindkmax
                            if locR == 1: shift -= 1
                            listOfIndicesD = GhostCells.getJoinDonorIndicesStruct__(prange,prangedonor,zoneDimD,dirD,
                                                                                    trirac,0,locS,dimPb, shift)
                    # Increments
                    absdirR=abs(dirR); absdirD=abs(dirD)
                    # list of directions of receptor zone
                    dir_rcv = [0]*dimPb; dir_rcv[absdirR-1] = 1
                    # list of direction of donor zone
                    dir_dnr = [0]*dimPb; dir_dnr[absdirD-1] = 1
                    if dimPb == 3:
                        incrR = GhostCells.increment__(dir_rcv[0],dir_rcv[1],dir_rcv[2], imr, jmr, kmr, 0)
                        incrD = GhostCells.increment__(dir_dnr[0],dir_dnr[1],dir_dnr[2], imd, jmd, kmd, 0)
                    elif dimPb == 2:
                        incrR = GhostCells.increment__(dir_rcv[0],dir_rcv[1], imr, jmr, 0)
                        incrD = GhostCells.increment__(dir_dnr[0],dir_dnr[1], imd, jmd, 0)
                    if dirR > 0: incrR = -incrR
                    if dirD < 0: incrD = -incrD

                    if dirR  ==-1: depth = rinddnr[0]
                    elif dirR== 1: depth = rinddnr[1]
                    elif dirR==-2: depth = rinddnr[2]
                    elif dirR== 2:  depth = rinddnr[3]
                    elif dirR==-3: depth = rinddnr[4]
                    else: depth = rinddnr[5]

                    res = Connector.connector.setInterpDataForGC(arrayOfIndicesR, listOfIndicesD, \
                                                                     dimPb, locR, depth, incrR, incrD)
                    # Stockage
                    vol = numpy.array([],numpy.float64)
                    prefix = 'ID'
                    if RotationAngle is not None:
                        val = Internal.getValue(RotationAngle)
                        if val[0]>0. or val[1]>0. or val[2]>0.: prefix='IDPERP'
                        else: prefix = 'IDPERM' 
                    if storage == 'direct':
                        createInterpRegion__(zp, zdonorname,res[0],res[1],res[3],res[2],vols,indicesExtrap,\
                                                 indicesOrphan,tag = 'Receiver',loc=loc,EXDir=EXdir,itype='abutting',\
                                                 prefix=prefix, RotationAngle=RotationAngle, RotationCenter=RotationCenter)

                    else:
                        createInterpRegion__(zdonorp, zname, res[1], res[0], res[3], res[2], vols, indicesExtrap,\
                                                 indicesOrphan, tag = 'Donor',loc=loc,EXDir=EXdir,itype='abutting',\
                                                 prefix=prefix,RotationAngle=RotationAngle, RotationCenter=RotationCenter)
    if storage == 'direct': return aR
    else: return aD

#===============================================================================
# General transfers: Match + Chimera + IBC
# Interpolation is applied to aR
# Beware: variables must be defined in topTreeD at nodes in order to
# be consistent with the computation of connectivity by setInterpData and
# setIBCData
# loc='nodes','centers' defines the location in aR of transfered values
# IN: variablesID=['var1','var2',...]: variables to be used in Chimera transfers
#                =[] : the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers
# IN: bcType (IBC only)  0: glissement
#                        1: adherence
#                        2: loi de paroi log
#                        3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE(,ronultideSA))
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t(,nutildeSA))
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p(,nutildeSA))
# IN: storage=-1/0/1: unknown/direct/inverse
# IN: loc = 'nodes' or 'centers': location of receiver zone field
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
#===============================================================================
def setInterpTransfers(aR, topTreeD,
                       variables=[],
                       variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                       bcType=0, varType=1, storage=-1, cellNVariable='',
                       Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                       Cs=0.3831337844872463, Ts=1.0):
    if variablesIBC is not None:
        nvarIBC = len(variablesIBC)
        if nvarIBC != 5 and nvarIBC != 6:
            raise ValueError("setInterpTransfers: length of variablesIBC must be equal to 5.")

    tR = Internal.copyRef(aR)
    # # checks
    # variablesD = C.getVarNames(topTreeD, excludeXYZ=True, loc='nodes')[0] # donor fields
    # variablesR = C.getVarNames(tR, excludeXYZ=True, loc=loc)[0]# receiver fields
    # if variablesD == []:
    #     raise ValueError("setInterpTransfers: no field to transfer.")
    # if variablesR == []:
    #     if loc == 'nodes': tR = C.addVars(tR,variablesD)
    #     else:
    #         variablesR = []
    #         for v in variablesD: variablesR+=['centers:'+v]
    #         tR = C.addVars(tR,variablesR)
    #

    compact = 0
    _setInterpTransfers(tR, topTreeD, variables=variables, variablesIBC=variablesIBC, 
                        bcType=bcType, varType=varType, storage=storage, compact=compact, 
                        cellNVariable=cellNVariable, Gamma=Gamma, Cv=Cv, MuS=MuS, Cs=Cs, Ts=Ts)

    return tR

#===============================================================================
# General transfers: Chimera + IBC - inplace version
# Interpolation is applied to aR
# Beware: variables must be defined in topTreeD at nodes in order to be
# consistent with the computation of connectivity by setInterpData and
# setIBCData
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variablesI =['var1','var2',...]: variables to be used in Chimera transfers
#                = None: the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers
# IN: bcType (IBC only) 0: glissement
#                       1: adherence
#                       2: loi de paroi log
#                       3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE)+ronultideSA
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t)+nutildeSA
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p)+nutildeSA
# IN: storage=-1/0/1: unknown/direct/inverse
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
#===============================================================================
def _setInterpTransfers(aR, topTreeD,
                        variables=[],
                        variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                        bcType=0, varType=1, storage=-1, compact=0, cellNVariable='',
                        Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                        Cs=0.3831337844872463, Ts=1.0):

    # Recup des donnees a partir des zones receveuses
    if storage < 1:
        # Dictionnaire pour optimisation
        znd = {}
        zones = Internal.getZones(topTreeD)
        for z in zones: znd[z[0]] = z

        zonesR = Internal.getZones(aR)
        for zr in zonesR:
            subRegions = Internal.getNodesFromType1(zr, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:2]                
                # test pour eviter parcours arbre inutile 
                if ((sname == 'ID' or sname == 'AD') and variables is not None) or (sname == 'IB' and variablesIBC is not None):
                   idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                   if idn is not None: # la subRegion decrit des interpolations
                       zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                       zoneRole = Internal.getValue(zoneRole)
                       if zoneRole == 'Receiver':
                           location = Internal.getNodeFromName1(s, 'GridLocation')
                           if location is not None: location = Internal.getValue(location)
                           Coefs = idn[1]
                           DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1] 
                           ListRcv   = Internal.getNodeFromName1(s,'PointList')[1]
                           ListDonor = Internal.getNodeFromName1(s,'PointListDonor')[1]
                           # Recup des champs du receveur
                           zdnrname = Internal.getValue(s)
                           zd = znd[zdnrname]

                           if location == 'CellCenter': loc = 1
                           else: loc = 0
                           # Transferts
                           if (sname == 'ID' or sname == 'AD'):
                                RotAngleNode = Internal.getNodeFromName1(s,'RotationAngle')
                                RotAngleX = 0.; RotAngleY = 0.; RotAngleZ = 0. 
                                if RotAngleNode is not None:
                                    RotAngleNode=RotAngleNode[1]
                                    RotAngleX = RotAngle[0]
                                    RotAngleY = RotAngle[1]
                                    RotAngleZ = RotAngle[2]

                                connector._setInterpTransfers(zr,zd,variables, ListRcv, ListDonor, DonorType, Coefs, loc, varType, compact, cellNVariable,
                                                              Internal.__GridCoordinates__, 
                                                              Internal.__FlowSolutionNodes__, 
                                                              Internal.__FlowSolutionCenters__, 
                                                              RotAngleX, RotAngleY, RotAngleZ)

                           elif sname == 'IB':
                               xPC = Internal.getNodeFromName1(s,'CoordinateX_PC')[1]                              
                               yPC = Internal.getNodeFromName1(s,'CoordinateY_PC')[1]
                               zPC = Internal.getNodeFromName1(s,'CoordinateZ_PC')[1]
                               xPW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                               yPW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                               zPW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                               xPI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                               yPI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                               zPI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                               density = Internal.getNodeFromName1(s,'Density')[1]
                               pressure = Internal.getNodeFromName1(s,'Pressure')[1]
                               utau = Internal.getNodeFromName1(s, 'utau')[1]
                               yplus = Internal.getNodeFromName1(s,'yplus')[1]
                               # Transferts
                               #print 'transfert IBC : zr ', zr[0], ' et donor : ', zd[0]
                               connector._setIBCTransfers(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs, 
                                                          xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                                                          bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts,
                                                          Internal.__GridCoordinates__, 
                                                          Internal.__FlowSolutionNodes__, 
                                                          Internal.__FlowSolutionCenters__)
                            
    # Recup des donnees a partir des zones donneuses
    if storage != 0:
        # Dictionnaire pour optimisation
        znr = {}
        zones = Internal.getZones(aR)
        for z in zones: znr[z[0]] = z
        zonesD = Internal.getZones(topTreeD)
        for zd in zonesD:
            subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:2]
                # test pour eviter parcours arbre inutile 
                if ((sname=='ID' or sname=='AD') and variables is not None) or (sname == 'IB' and variablesIBC is not None):
                   idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                   if idn is not None:
                       zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                       zoneRole = Internal.getValue(zoneRole)
                       if zoneRole == 'Donor':
                           location = Internal.getNodeFromName1(s, 'GridLocation') #localisation des donnees des rcvr
                           if location is not None: location = Internal.getValue(location)
                           Coefs = idn[1]
                           DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                           ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                           ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                           # Recup des champs du receveur
                           zrcvname = Internal.getValue(s)
                           ##zr = Internal.getNodesFromName2(aR, zrcvname)
                           zr = znr.get(zrcvname, None)
                           if zr is not None:
                                if location == 'CellCenter': loc = 1
                                else: loc = 0
                                # Transferts
                                if sname == 'ID' or sname == 'AD':
                                    RotAngleNode = Internal.getNodeFromName1(s,'RotationAngle')
                                    RotAngleX = 0.; RotAngleY = 0.; RotAngleZ = 0. 
                                    if RotAngleNode is not None:
                                        RotAngleNode=RotAngleNode[1]
                                        RotAngleX = RotAngleNode[0]
                                        RotAngleY = RotAngleNode[1]
                                        RotAngleZ = RotAngleNode[2]

                                    connector._setInterpTransfers(zr, zd, variables, ListRcv, ListDonor, DonorType, Coefs,loc, varType, compact, cellNVariable,
                                                                  Internal.__GridCoordinates__, 
                                                                  Internal.__FlowSolutionNodes__, 
                                                                  Internal.__FlowSolutionCenters__,
                                                                  RotAngleX, RotAngleY, RotAngleZ)
                                elif sname == 'IB':
                                    xPC = Internal.getNodeFromName1(s,'CoordinateX_PC')[1]
                                    yPC = Internal.getNodeFromName1(s,'CoordinateY_PC')[1]
                                    zPC = Internal.getNodeFromName1(s,'CoordinateZ_PC')[1]
                                    xPW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                                    yPW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                                    zPW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                                    xPI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                                    yPI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                                    zPI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                                    Density = Internal.getNodeFromName1(s,'Density')[1]
                                    Pressure= Internal.getNodeFromName1(s,'Pressure')[1]
                                    utau    = Internal.getNodeFromName1(s, 'utau')[1]
                                    yplus   = Internal.getNodeFromName1(s,'yplus')[1]
                                    #  print 'transfert IBC : zr ', zr[0], ' et donor : ', zd[0] 
                                    connector._setIBCTransfers(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs, 
                                                               xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, Density, Pressure, utau, yplus,
                                                               bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts,
                                                               Internal.__GridCoordinates__, 
                                                               Internal.__FlowSolutionNodes__, 
                                                               Internal.__FlowSolutionCenters__)
    return None

#===============================================================================
# General transfers: Chimera + IBC - inplace version optimiser par arbre tc compacte par zone donneuse
# Interpolation is applied to aR 
# Beware: variables must be defined in topTreeD at nodes in order to be 
# consistent with the computation of connectivity by setInterpData and 
# setIBCData 
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variablesI =['var1','var2',...]: variables to be used in Chimera transfers
#                = None: the whole FlowSolutionNodes variables in topTreeD are transferred 
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers 
# IN: bcType (IBC only) 0: glissement
#                       1: adherence
#                       2: loi de paroi log
#                       3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE)+ronultideSA
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t)+nultideSA
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p)+nultideSA
# IN: storage=-1/0/1: unknown/direct/inverse
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
#===============================================================================
def __setInterpTransfers(aR, topTreeD, 
                         variables=[], 
                         variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'], 
                         bcType=0, varType=1, storage=-1, compact=0,
                         Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08, 
                         Cs=0.3831337844872463, Ts=1.0):

    # Recup des donnees a partir des zones receveuses    
    if (storage != 1 or compact==0):
        print 'erreur __setInterpTransfers: Mode receveur a coder. Mode compact obligatoire: compact=1'
                            
    #test pour savoir si appel transfert ou ibc
    flagibc = 1
    if variables is not None: flagibc = 0

    # Recup des donnees a partir des zones donneuses
    zones  = Internal.getZones(aR)
    zonesD = Internal.getZones(topTreeD)
    for zd in zonesD:
            param_int = Internal.getNodeFromName1(zd, 'Parameter_int')[1]
            param_real= Internal.getNodeFromName1(zd, 'Parameter_real')[1]
            if   (param_int[0] != 0 and flagibc == 0):
                connector.__setInterpTransfers(zones, zd, variables   , param_int, param_real, varType, compact, flagibc, bcType, Gamma, Cv, MuS, Cs, Ts )
            elif (param_int[1] != 0 and flagibc == 1):
                connector.__setInterpTransfers(zones, zd, variablesIBC, param_int, param_real, varType, compact, flagibc, bcType, Gamma, Cv, MuS, Cs, Ts )
    return None
 

#==============================================================================
# Mise a plat (compactage) arbre donneur au niveau de la base
# fonctionne avec ___setInterpTransfer
#==============================================================================
def miseAPlatDonnorTree__(zones, tc, procDict=None, procList=None):
    size_int  = 0
    size_real = 0
    listproc  = []
    rac       = []
    rac_inst  = []
    sizeI     = []
    sizeR     = []
    sizeNbD   = []
    sizeType  = []
    nrac      = 0

    ordered_subRegions =[]
    No_zoneD =[]
    MeshTypeD=[]
    inst = {}
    numero_max =-100000000
    numero_min = 100000000

    zones_tc = Internal.getZones(tc)
    c        = 0
    for z in zones_tc:
      subRegions  =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
      meshtype   = 1
      zonetype   = Internal.getNodeFromType1(z, 'ZoneType_t')
      tmp        = Internal.getValue(zonetype)
      if tmp != "Structured": meshtype = 2
      for s in subRegions:
         #tri des pas de temps instationnaire
         #  1) les stationnaires
         #  2) les instationnaires regroupes par pas de temps
         if '#' not in s[0]: 
              ordered_subRegions.append(s)
              No_zoneD.append(c)
              MeshTypeD.append(meshtype)
         else: 
            numero_iter = int( s[0].split('#')[1].split('_')[0] )
            if numero_iter < numero_min : numero_min = numero_iter
            if numero_iter > numero_max : numero_max = numero_iter

            if numero_iter in inst.keys():
                sub = inst[ numero_iter ][0]
                sub = sub + [s]
                Noz = inst[ numero_iter ][1]
                Noz = Noz + [c]
                mesh= inst[ numero_iter ][2]
                mesh= mesh+ [meshtype]
                inst[ numero_iter ]=  [ sub , Noz , mesh ]
            else:
                inst[ numero_iter ]= [ [s],[c],[meshtype] ]

         TimeLevelNumber=  len(inst)
         if TimeLevelNumber != 1+numero_max-numero_min and len(inst)!= 0: 
              print "ERROR miseAPlatDonnorTree__: missing timestep in tc", numero_max,numero_min, TimeLevelNumber

         count_ID  = 0
         count_IBC = 0
         #alloc memoire
         pointlist     =  Internal.getNodeFromName1(s, 'PointList')
         pointlistD    =  Internal.getNodeFromName1(s, 'PointListDonor')
         InterpD       =  Internal.getNodeFromName1(s, 'InterpolantsDonor')
         Interptype    =  Internal.getNodeFromName1(s, 'InterpolantsType')
         RotationAngle =  Internal.getNodeFromName1(s, 'RotationAngle')
         RotationCenter=  Internal.getNodeFromName1(s, 'RotationCenter')

         Nbpts        =  numpy.shape(pointlist[ 1])[0]
         Nbpts_D      =  numpy.shape(pointlistD[1])[0]
         Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

         sname = s[0][0:2]
         utau =  Internal.getNodeFromName1(s, 'utau')
         if sname == 'IB': 
           if utau is None:
             # Creation des numpy d extraction
             pressNP = numpy.zeros((Nbpts_D),numpy.float64)
             densNP  = numpy.zeros((Nbpts_D),numpy.float64)
             utauNP  = numpy.zeros((Nbpts_D),numpy.float64)
             yplusNP = numpy.zeros((Nbpts_D),numpy.float64)

             Internal.createUniqueChild(s, 'Density'  , 'DataArray_t', densNP )
             Internal.createUniqueChild(s, 'Pressure' , 'DataArray_t', pressNP)
             Internal.createUniqueChild(s, 'utau'  , 'DataArray_t', utauNP )
             Internal.createUniqueChild(s, 'yplus' , 'DataArray_t', yplusNP )

             utau =  Internal.getNodeFromName1(s, 'utau')

         # on recupere le nombre de type different
         #typecell = Interptype[1][0]
         #Nbtype = [ typecell ]
         #for i in xrange(Nbpts_D):
         #  if Interptype[1][i] not in Nbtype: Nbtype.append(Interptype[1][i])
         #print 'nb type',  len(Nbtype), s[0],z[0], Nbtype
         nbType = numpy.unique(Interptype[1])
         nbTypeSize = nbType.size
         #print len(Nbtype), Nbtype2.size

         sizeIBC    =  0
         ntab_IBC   = 11
         if utau is not None: ntab_IBC = 13
         if sname == 'IB': 
            sizeIBC   = Nbpts_D*ntab_IBC
            count_IBC += 1
         else: 
            count_ID  += 1

         #periodicite azymuthal
         rotation = 0
         if RotationAngle is not None : rotation +=3 
         if RotationCenter is not None : rotation +=3 

         nrac  =  nrac + 1
         zRname = Internal.getValue(s)
         proc = 0
         if procDict is not None: proc = procDict[zRname]
         if proc not in listproc: 
                listproc.append(proc)
                rac.append(1)
                if '#' in s[0]: rac_inst.append(1)
                else          : rac_inst.append(0)
                sizeI.append(    Nbpts_D*2     + Nbpts   + nbTypeSize+1 )
                sizeR.append(    Nbpts_InterpD + sizeIBC + rotation      )
                sizeNbD.append(  Nbpts_D                                 )
                sizeType.append( Nbpts_D                 + nbTypeSize+1 )
         else:
                pos           = listproc.index(proc)
                rac[pos]      = rac[pos] + 1
                if '#' in s[0]: rac_inst[pos]= rac_inst[pos] + 1
                sizeI[pos]    = sizeI[pos]     + Nbpts_D*2     + Nbpts   + nbTypeSize+1
                sizeR[pos]    = sizeR[pos]     + Nbpts_InterpD + sizeIBC + rotation
                sizeNbD[pos]  = sizeNbD[pos]   + Nbpts_D
                sizeType[pos] = sizeType[pos]  + Nbpts_D                 + nbTypeSize+1

      c+=1

    base     = Internal.getNodeFromType1(tc, 'CGNSBase_t')  # noeud
    model    = 'NSLaminar'
    a        = Internal.getNodeFromName2(base, 'GoverningEquations')
    if a is not None: model = Internal.getValue(a)

    NbP2P     = len(listproc)
    sizeproc  = []
    ntab_int  = 16
    for i in xrange(NbP2P): sizeproc.append(5 + TimeLevelNumber*2 + ntab_int*rac[i] + sizeI[i])
           
    size_int =  1 + NbP2P + sum(sizeproc)
    size_real=  sum(sizeR)

    param_int  = numpy.empty(size_int , dtype=numpy.int32  )
    param_real = numpy.empty(size_real, dtype=numpy.float64)
    Internal.createUniqueChild(tc, 'Parameter_int' , 'DataArray_t', param_int )
    if size_real !=0 : 
        Internal.createUniqueChild(tc, 'Parameter_real', 'DataArray_t', param_real)

    # Dictionnaire pour optimisation
    znd = []
    for z in zones: znd.append(z[0])

    #
    #initialisation numpy
    #
    param_int[0] = NbP2P
    size_ptlist = []
    size_ptlistD= []
    size_ptType = []
    nb_rac      = []
    size_coef   = []
    adr_coef    = []   # pour cibler debut de echange dans param_real

    shift_coef  =0
    shift       =0
    for i in xrange(NbP2P):
       adr_coef.append(shift_coef)                    #adresse echange dans param_real
       shift_coef = shift_coef + sizeR[i]

       param_int[i+1] = 1 + NbP2P + shift              #adresse echange
       shift          =  shift  + sizeproc[i]
       size_ptlist.append(0)
       size_ptlistD.append(0)
       size_ptType.append(0)
       size_coef.append(0)
       nb_rac.append(0)

    for iter in range(numero_min,numero_max+1): 
        ordered_subRegions =  ordered_subRegions + inst[ iter ][0]
        No_zoneD           =  No_zoneD           + inst[ iter ][1]
        MeshTypeD          =  MeshTypeD          + inst[ iter ][2]

    #loop sur les raccord tries
    c        = 0
    Nbtot    = 0
    for s in ordered_subRegions:

       NozoneD  = No_zoneD[c]
       meshtype = MeshTypeD[c]

       zRname = Internal.getValue(s)
       proc = 0
       if procDict is not None: proc = procDict[zRname]
       pos  = listproc.index(proc)


       pt_ech = param_int[ pos+1 ]                 # adresse debut raccord pour l'echange pos
       pt_coef= adr_coef[pos] + size_coef[pos]     # adresse debut coef 

       pointlist     =  Internal.getNodeFromName1(s, 'PointList')
       pointlistD    =  Internal.getNodeFromName1(s, 'PointListDonor')
       Interptype    =  Internal.getNodeFromName1(s, 'InterpolantsType')
       InterpD       =  Internal.getNodeFromName1(s, 'InterpolantsDonor')
       RotationAngle =  Internal.getNodeFromName1(s, 'RotationAngle')
       RotationCenter=  Internal.getNodeFromName1(s, 'RotationCenter')
       Nbpts         =  numpy.shape(pointlist[ 1])[0]
       Nbpts_D       =  numpy.shape(pointlistD[1])[0]
       Nbpts_InterpD =  numpy.shape(InterpD[ 1  ])[0]

       param_int[ pt_ech    ] = proc
       param_int[ pt_ech +1 ] = rac[pos]
       param_int[ pt_ech +2 ] = rac_inst[pos]
       nrac_steady            = rac[pos] - rac_inst[pos]

       param_int[ pt_ech +3 ] = TimeLevelNumber
       nrac_inst_deb  =  nrac_steady
       for i in xrange(TimeLevelNumber):
            # len(inst[i][0])  = nb de raccord instationnaire pour le temps i
            nrac_inst_fin  = nrac_inst_deb + len(inst[i][0])

            param_int[ pt_ech +4 + i                  ] = nrac_inst_deb
            param_int[ pt_ech +4 + i + TimeLevelNumber] = nrac_inst_fin

            nrac_inst_deb  = nrac_inst_fin
 
       #iadr = pt_ech + 4 + nb_rac[pos]   # ptr echange + dest + nrac + norac
       iadr = pt_ech + 4 + TimeLevelNumber*2 + nb_rac[pos]   # ptr echange + dest + nrac + norac

       param_int[ iadr            ] = Nbpts
       param_int[ iadr + rac[pos] ] = Nbpts_InterpD

       #on recupere le nombre de type different
       typecell = Interptype[1][0]
       Nbtype= [ typecell ]
       for i in xrange(Nbpts_D):
          if Interptype[1][i] not in Nbtype: Nbtype.append(Interptype[1][i])

       #Si le type zero existe, on le place a la fin: sinon adressage openmp=boom dans donnorPts
       if 0 in Nbtype: Nbtype += [Nbtype.pop( Nbtype.index( 0 ) )]

       param_int[ iadr +rac[pos]*2 ] = len(Nbtype)

       param_int[ iadr +rac[pos]*3 ] = 0
       size_IBC = 0

       sname = s[0][0:2]
       xc=None;yc=None;zc=None; xi=None;yi=None;zi=None; xw=None;yw=None;zw=None;density=None;pressure=None;utau=None;yplus=None;
       ptxc=0;ptyc=0;ptzc=0;ptxi=0;ptyi=0;ptzi=0;ptxw=0;ptyw=0;ptzw=0;ptdensity=0;ptpressure=0; ptutau=0;ptyplus=0;
       if sname == 'IB': 

           param_int[ iadr +rac[pos]*3 ]  = 1

           xc        = Internal.getNodeFromName1(s , 'CoordinateX_PC')
           yc        = Internal.getNodeFromName1(s , 'CoordinateY_PC')
           zc        = Internal.getNodeFromName1(s , 'CoordinateZ_PC')
           xi        = Internal.getNodeFromName1(s , 'CoordinateX_PI')
           yi        = Internal.getNodeFromName1(s , 'CoordinateY_PI')
           zi        = Internal.getNodeFromName1(s , 'CoordinateZ_PI')
           xw        = Internal.getNodeFromName1(s , 'CoordinateX_PW')
           yw        = Internal.getNodeFromName1(s , 'CoordinateY_PW')
           zw        = Internal.getNodeFromName1(s , 'CoordinateZ_PW')
           density   = Internal.getNodeFromName1(s , 'Density')
           pressure  = Internal.getNodeFromName1(s , 'Pressure')
           utau      = Internal.getNodeFromName1(s , 'utau')
           yplus     = Internal.getNodeFromName1(s , 'yplus')

           ptxc      = pt_coef + Nbpts_InterpD
           ptyc      = pt_coef + Nbpts_InterpD + Nbpts_D
           ptzc      = pt_coef + Nbpts_InterpD + Nbpts_D*2
           ptxi      = pt_coef + Nbpts_InterpD + Nbpts_D*3
           ptyi      = pt_coef + Nbpts_InterpD + Nbpts_D*4
           ptzi      = pt_coef + Nbpts_InterpD + Nbpts_D*5
           ptxw      = pt_coef + Nbpts_InterpD + Nbpts_D*6
           ptyw      = pt_coef + Nbpts_InterpD + Nbpts_D*7
           ptzw      = pt_coef + Nbpts_InterpD + Nbpts_D*8
           ptdensity = pt_coef + Nbpts_InterpD + Nbpts_D*9
           ptpressure= pt_coef + Nbpts_InterpD + Nbpts_D*10

           size_IBC    = 11*Nbpts_D

           if utau is not None:
               ptutau    = pt_coef + Nbpts_InterpD + Nbpts_D*11
               ptyplus   = pt_coef + Nbpts_InterpD + Nbpts_D*12
               size_IBC  = 13*Nbpts_D
                   

       tmp = Internal.getNodeFromName1(s , 'ZoneRole')
       if tmp[1][0] == 'D': param_int[ iadr +rac[pos]*4 ] = 0   # role= Donor
       else               : param_int[ iadr +rac[pos]*4 ] = 1   # role= Receiver
            
       param_int[ iadr +rac[pos]*5 ] = NozoneD                    # No zone donneuse

       lst                              = pt_ech +5 +ntab_int*rac[pos] + sizeNbD[pos] + sizeType[pos] + size_ptlist[pos]  
       param_int[ iadr +rac[pos]*6    ] = lst                                                                      # PointlistAdr
       ptTy                             = pt_ech +5 +ntab_int*rac[pos] + sizeNbD[pos] + size_ptType[pos]
       param_int[ iadr +rac[pos]*7    ] = ptTy                                                                     # TypAdr
       lstD                             = pt_ech +5 +ntab_int*rac[pos] + size_ptlistD[pos]
       param_int[ iadr +rac[pos]*12 +1] = lstD                                                                     # PointlistDAdr

       Nbtot = Nbtot + Nbpts

       param_int[ ptTy  ] = len(Nbtype)
       noi       = 0
       nocoef    = 0
       sizecoef  = 0
       shift_typ = 1 + len(Nbtype)
       ctyp      = 0
       l0        = 0

       #recopie dans tableau a plat + tri par type
       if len(Nbtype)==1:
           tri_monoType(Nbpts_D, Nbpts,Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp, ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                        Interptype, pointlist, pointlistD, param_int,
                        ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, ptdensity,ptpressure,ptutau,ptyplus,
                        xc,yc,zc,xi,yi,zi,xw,yw,zw, density,pressure,utau,yplus,InterpD,param_real )
       else:
           tri_multiType(Nbpts_D,Nbpts,Nbpts_InterpD, meshtype, noi, lst,lstD,l0,ctyp, ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                         Interptype, pointlist, pointlistD, param_int,
                         ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, ptdensity,ptpressure,ptutau,ptyplus,
                         xc,yc,zc,xi,yi,zi,xw,yw,zw, density,pressure,utau,yplus,InterpD,param_real )

       pointlist[ 1] = param_int[ lst             : lst              + Nbpts         ]    # supression numpy initial pointlist
       Interptype[1] = param_int[ ptTy + shift_typ: ptTy + shift_typ + Nbpts_D       ]    # supression numpy initial interpolantType
       pointlistD[1] = param_int[ lstD            : lstD             + Nbpts_D       ]    # supression numpy initial pointlistDonnor
       InterpD[   1] = param_real[ pt_coef        : pt_coef          + Nbpts_InterpD ]    # supression numpy initial interpDonnor

       #if (s[0] == 'ID_cart3' and z[0]=='cart1'): print 'verif',  InterpD[   1][0], pt_coef,numpy.shape(InterpD[ 1  ])

       if sname == 'IB':
           xc[1]       = param_real[ ptxc: ptxc+ Nbpts_D ]
           yc[1]       = param_real[ ptyc: ptyc+ Nbpts_D ]
           zc[1]       = param_real[ ptzc: ptzc+ Nbpts_D ]
           xi[1]       = param_real[ ptxi: ptxi+ Nbpts_D ]
           yi[1]       = param_real[ ptyi: ptyi+ Nbpts_D ]
           zi[1]       = param_real[ ptzi: ptzi+ Nbpts_D ]                                      # supression numpy initial IBC
           xw[1]       = param_real[ ptxw: ptxw+ Nbpts_D ]
           yw[1]       = param_real[ ptyw: ptyw+ Nbpts_D ]
           zw[1]       = param_real[ ptzw: ptzw+ Nbpts_D ]
           density[1]  = param_real[ ptdensity : ptdensity + Nbpts_D ]
           pressure[1] = param_real[ ptpressure: ptpressure+ Nbpts_D ]
           if utau is not None:
               utau[1]  = param_real[ ptutau : ptutau  + Nbpts_D ]
               yplus[1] = param_real[ ptyplus: ptyplus + Nbpts_D ]

       param_int[ iadr +rac[pos]*8 ] = adr_coef[pos] + size_coef[pos]          # PtcoefAdr
          
       iadr = iadr +1
       param_int[ iadr +rac[pos]*8 ] = rac[pos]                  # nrac pour mpi

       tmp = Internal.getNodeFromName1(s , 'GridLocation')
       if tmp[1][4] == 'C': param_int[ iadr +rac[pos]*9 ] = 1   # location= CellCenter
       else               : param_int[ iadr +rac[pos]*9 ] = 0   # location= CellVertex

       param_int[ iadr +rac[pos]*10 ] = Nbpts_D

       #chercher No zone receveuse grace a tc ou dico (si mpi)
       if procDict is None:
         param_int[ iadr +rac[pos]*11 ] = znd.index( zRname )         # No zone receveuse
       else:
         param_int[ iadr +rac[pos]*11  ]= procList[proc].index( zRname )  # No zone raccord
         #no_zone =0
         #for cle in procDict:
         #   if (cle == zRname ): 
         #      param_int[ iadr +rac[pos]*11  ]= no_zone  # No zone raccord
         #print "AV dic", param_int[ iadr +rac[pos]*11  ] , proc
         #   if (procDict[cle] == proc ): no_zone = no_zone +1

       #print 'rac', s[0], 'zoneR=', zRname, param_int[ iadr +rac[pos]*11 ], 'NozoneD=', zones_tc[No_zoneD[c]][0]

       # ATTENTION STEPHANIE A TOUT COMMENTE
       #a = Internal.getNodeFromName2(zones_tc[ param_int[ iadr +rac[pos]*11  ]  ], 'GoverningEquations')
       #if a is not None: model = Internal.getValue(a)

       #print 'model=',model,'zoneR',zones_tc[param_int[ iadr +rac[pos]*11  ]][0], 'NoR=', param_int[ iadr +rac[pos]*11  ], 'NoD=', c
       if model=='NSTurbulent': neq_loc = 6
       else                   : neq_loc = 5
            
       tmp =  Internal.getNodeFromName1(s , 'RANSLES')
       if tmp is not None: param_int[ iadr +rac[pos]*13  ] = min (5, neq_loc)   # RANSLES
       else:               param_int[ iadr +rac[pos]*13  ] = neq_loc
       
       # raccord periodique avec rotation
       if RotationAngle is not None: 
             param_int[ iadr +rac[pos]*14  ] = 1
             shiftRotation                   = 6
             ptdeb =   pt_coef + Nbpts_InterpD   
             param_real[ ptdeb   : ptdeb+3 ] = RotationAngle[1][0:3]
             param_real[ ptdeb+3 : ptdeb+6 ] = RotationCenter[1][0:3]
             RotationAngle[1]  =    param_real[ ptdeb   : ptdeb+3]                                    
             RotationCenter[1] =    param_real[ ptdeb+3 : ptdeb+6]                                    

       else: 
             param_int[ iadr +rac[pos]*14  ] = 0
             shiftRotation                   = 0


       # raccord instationnaire
       param_int[ iadr +rac[pos]*15  ] =-1
       if '#' in s[0]:
           numero_iter = int( s[0].split('#')[1].split('_')[0] )
           param_int[ iadr +rac[pos]*15  ] = numero_iter
          
       #print 'model=', model, 'tmp=', tmp, 'neq_loc=', neq_loc

       ### Verifier position et chois entre Nbpts et NbptsD 
       size_ptlistD[pos] = size_ptlistD[pos] +  Nbpts_D
       size_ptlist[pos]  = size_ptlist[pos]  +  Nbpts
       size_ptType[pos]  = size_ptType[pos]  +  Nbpts_D + len(Nbtype)+1


       size_coef[pos] = size_coef[pos] + Nbpts_InterpD + size_IBC + shiftRotation
       nb_rac[pos]    = nb_rac[pos] + 1

       c += 1

    return None

#==============================================================================
# tri multitype
#==============================================================================
def tri_multiType(Nbpts_D, Nbpts, Nbpts_InterpD,meshtype, noi, lst,lstD,l0,ctyp,ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                  Interptype, pointlist, pointlistD, param_int,
                  ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, ptdensity,ptpressure,ptutau,ptyplus,
                  xc,yc,zc,xi,yi,zi,xw,yw,zw, density,pressure,utau,yplus,InterpD,param_real ):

  for ntype in Nbtype:
    noi_old   = 0
    nocoef_old= 0
    l         = 0
    for i in xrange(Nbpts_D):

       ltype = Interptype[1][i]
       if meshtype == 1:
         if ltype == 1: sizecoef=1
         elif ltype == 2: sizecoef=8
         elif ltype == 3: sizecoef=9
         elif ltype == 4: sizecoef=8
         elif ltype == 5: sizecoef=15
         elif ltype == 22: sizecoef=4
       else:
         if ltype == 1: sizecoef=1
         elif ltype == 4: sizecoef=4

       if ltype == ntype:
            #recopieinterpolantType
            param_int[ ptTy + shift_typ + l + l0 ] = ltype

            # recopie pointlist
            if ntype != 0:
              param_int[  lst + noi] = pointlist[ 1][ noi_old ]
              noi     = noi     + 1
              noi_old = noi_old + 1
            else:
              ncfLoc   = pointlist[ 1][ noi_old ]
              sizecoef = ncfLoc
              param_int[  lst + noi] = ncfLoc
              param_int[ lst+noi+1: lst+noi+1+ncfLoc] = pointlist[1][ noi_old+1: noi_old+1+ncfLoc]
              noi     = noi     + 1 + ncfLoc
              noi_old = noi_old + 1 + ncfLoc

            #recopie pointListDonnor
            param_int[ lstD +  l + l0] = pointlistD[1][i]
            #recopie Maillage IBC
            if sname == 'IB':
               param_real[ ptxc      + l + l0 ]= xc[1][i]
               param_real[ ptyc      + l + l0 ]= yc[1][i]
               param_real[ ptzc      + l + l0 ]= zc[1][i]
               param_real[ ptxi      + l + l0 ]= xi[1][i]
               param_real[ ptyi      + l + l0 ]= yi[1][i]
               param_real[ ptzi      + l + l0 ]= zi[1][i]
               param_real[ ptxw      + l + l0 ]= xw[1][i]
               param_real[ ptyw      + l + l0 ]= yw[1][i]
               param_real[ ptzw      + l + l0 ]= zw[1][i]
               param_real[ ptdensity + l + l0 ]= density[1][i]
               param_real[ ptpressure+ l + l0 ]= pressure[1][i]
               if utau is not None:
                   param_real[ ptutau    + l + l0 ]= utau[1][i]
                   param_real[ ptyplus   + l + l0 ]= yplus[1][i]

            #recopie  InterpD
            param_real[ pt_coef + nocoef: pt_coef + nocoef+sizecoef] = InterpD[1][ nocoef_old: nocoef_old+sizecoef]
            nocoef     = nocoef     + sizecoef
            nocoef_old = nocoef_old + sizecoef
            l += 1

       else:
            if ntype != 0:
                   noi_old = noi_old + 1
            else:           
                   ncfLoc  = pointlist[ 1][ noi_old ]
                   noi_old = noi_old + 1 + ncfLoc

            nocoef_old += sizecoef

    l0 = l0 + l

    param_int[ ptTy + ctyp +1 ] = l
    ctyp                        = ctyp +1

  return None
#==============================================================================
# tri monotype
#==============================================================================
def tri_monoType(Nbpts_D, Nbpts, Nbpts_InterpD,meshtype, noi, lst,lstD,l0,ctyp,ptTy,shift_typ,pt_coef,nocoef,sname,Nbtype,
                  Interptype, pointlist, pointlistD, param_int,
                  ptxc,ptyc,ptzc,ptxi,ptyi,ptzi,ptxw,ptyw,ptzw, ptdensity,ptpressure,ptutau,ptyplus,
                  xc,yc,zc,xi,yi,zi,xw,yw,zw, density,pressure,utau,yplus,InterpD,param_real ):

  ntype     = Nbtype[0]
  noi_old   = 0
  nocoef_old= 0
  l         = 0
  ltype     = Interptype[1][0]

  #recopieinterpolantType
  ideb =  ptTy + shift_typ
  val = float(ltype)
  connector.initNuma( None, param_int, ideb, Nbpts_D, 1, val )

  # recopie pointlist
  ideb =lst
  connector.initNuma( pointlist[1] , param_int, ideb, Nbpts , 1, val   )
  #recopie pointListDonnor
  ideb =lstD
  connector.initNuma( pointlistD[1] , param_int, ideb, Nbpts_D , 1, val   )
  #recopie Maillage IBC
  if sname == 'IB':
       connector.initNuma( xc[1] , param_real, ptxc, Nbpts_D , 0, val   )
       connector.initNuma( yc[1] , param_real, ptyc, Nbpts_D , 0, val   )
       connector.initNuma( zc[1] , param_real, ptzc, Nbpts_D , 0, val   )
       connector.initNuma( xi[1] , param_real, ptxi, Nbpts_D , 0, val   )
       connector.initNuma( yi[1] , param_real, ptyi, Nbpts_D , 0, val   )
       connector.initNuma( zi[1] , param_real, ptzi, Nbpts_D , 0, val   )
       connector.initNuma( xw[1] , param_real, ptxw, Nbpts_D , 0, val   )
       connector.initNuma( yw[1] , param_real, ptyw, Nbpts_D , 0, val   )
       connector.initNuma( zw[1] , param_real, ptzw, Nbpts_D , 0, val   )
       
       connector.initNuma( density[1] , param_real, ptdensity , Nbpts_D , 0, val   )
       connector.initNuma(pressure[1] , param_real, ptpressure, Nbpts_D , 0, val   )
       if utau is not None:
           connector.initNuma( utau[1] , param_real, ptutau , Nbpts_D , 0, val   )
           connector.initNuma(yplus[1] , param_real, ptyplus, Nbpts_D , 0, val   )

  #recopie  InterpD
  connector.initNuma( InterpD[1] , param_real, pt_coef , Nbpts_InterpD , 0, val   )

  param_int[ ptTy + ctyp +1 ] = Nbpts_D

  return None
#==============================================================================
# Mise a plat (compactage) arbre donneur au niveau de la zone donneuse
# fonctionne avec __setInterpTransfer
#==============================================================================
def miseAPlatDonnorZone__(zones, tc, procDict):
       zones_tc = Internal.getZones(tc)
       for z in zones_tc:
         racs      =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
         size_int  = 0
         size_real = 0
         count_ID  = 0
         count_IBC = 0
         #alloc memoire
         for rac in racs:
            pointlist    =  Internal.getNodeFromName1(rac, 'PointList')
            pointlistD   =  Internal.getNodeFromName1(rac, 'PointListDonor')
            InterpD      =  Internal.getNodeFromName1(rac, 'InterpolantsDonor')
            utau         =  Internal.getNodeFromName1(rac, 'utau')

            sizeIBC    =  0
            ntab_IBC   = 11
            if utau is not None: ntab_IBC = 13

            Nbpts        =  numpy.shape(pointlist[ 1])[0]
            Nbpts_D      =  numpy.shape(pointlistD[1])[0]
            Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

            size_int   =  size_int + 7 + Nbpts_D*2 + Nbpts
            size_real  =  size_real+   Nbpts_InterpD
            sname = rac[0][0:2]
            if sname == 'IB': 
               size_real  =  size_real +   Nbpts_D*ntab_IBC
               count_IBC = count_IBC +1;
            else: 
               count_ID  = count_ID  +1;
            #print 'nbpt, nbpt_donor', sname,Nbpts,Nbpts_InterpD

         size_int = size_int + 2 + (count_IBC + count_ID)*2  # 2: nbr rac ID et IBC, stockage adresse debut raccord et coef
         param_int  = numpy.empty(size_int , dtype=numpy.int32  )
         param_real = numpy.empty(size_real, dtype=numpy.float64)
	 Internal.createUniqueChild(z, 'Parameter_int' , 'DataArray_t', param_int )
         if size_real !=0 :
	    Internal.createUniqueChild(z, 'Parameter_real', 'DataArray_t', param_real)# recopie pointlis

         #print 'size int et real', size_int, size_real
         #print 'ID et IBC', count_ID, count_IBC

         param_int[0] = count_ID
         param_int[1] = count_IBC
         #initialisation
         c        = 0
         size_rac = 0
         size_coef= 0
         for rac in racs:
            pt_rac = 2 + (count_ID + count_IBC)*2 + size_rac # adresse debut raccord 
            pt_coef= size_coef                               # adresse debut coef 
            #print 'indice', pt_rac , Nbpts,count_ID,count_IBC, size_rac

            pointlist    =  Internal.getNodeFromName1(rac, 'PointList')
            pointlistD   =  Internal.getNodeFromName1(rac, 'PointListDonor')
            Interptype   =  Internal.getNodeFromName1(rac, 'InterpolantsType')
            InterpD      =  Internal.getNodeFromName1(rac, 'InterpolantsDonor')
            Nbpts        =  numpy.shape(pointlist[ 1])[0]
            Nbpts_D      =  numpy.shape(pointlistD[1])[0]
            Nbpts_InterpD=  numpy.shape(InterpD[ 1  ])[0]

            param_int[ 2+c                    ] = pt_rac 
            param_int[ 2+c+count_ID+count_IBC ] = pt_coef
            param_int[ pt_rac                 ] = Nbpts
            param_int[ pt_rac +1              ] = Nbpts_D
            param_int[ pt_rac +2              ] = Nbpts_InterpD

            typecell = Interptype[1][0]
            #on cree un tableau qui contient les elements non egaux au premier element
            b = Interptype[1] [ Interptype[1] != typecell ]
            if (len(b) == 0):   param_int[ pt_rac +3 ] = 0    # type homogene
            else:                param_int[ pt_rac +3 ] = 1    # type melange

            tmp =  Internal.getNodeFromName1(rac, 'ZoneRole')
            if (tmp[1][0] == 'D'): param_int[ pt_rac +4 ] = 0           # role= Donor
            else                 : param_int[ pt_rac +4 ] = 1           # role= Receiver
            tmp =  Internal.getNodeFromName1(rac, 'GridLocation')
            if (tmp[1][4] == 'C'): param_int[ pt_rac +5 ] = 1           # location= CellCenter
            else                 : param_int[ pt_rac +5 ] = 0           # location= CellVertex
            
            zrcvname = Internal.getValue(rac)
            no_zone =0
            for z in zones:
               if (z[0] == zrcvname ): param_int[ pt_rac +6 ]= no_zone  # No zone raccord                    
               no_zone = no_zone +1

            ideb =  pt_rac +7
            param_int[ ideb:ideb + Nbpts   ] = pointlist[ 1][0:Nbpts  ]           # recopie pointlist
            pointlist[ 1]                    = param_int[ ideb : ideb + Nbpts ]   # supression numpy initial

            ideb =  ideb  + Nbpts
            param_int[ ideb:ideb+ Nbpts_D  ] = pointlistD[1][0:Nbpts_D]           # recopie pointlistdonor
            pointlistD[ 1]                   = param_int[ ideb : ideb + Nbpts_D ] # supression numpy initial

            ideb =  ideb  + Nbpts_D
            param_int[ ideb:ideb + Nbpts_D ] = Interptype[1][0:Nbpts_D]           #recopieinterpolantType 
            Interptype[ 1]                   = param_int[ ideb : ideb + Nbpts_D ] # supression numpy initial

            size_rac   =  size_rac + 7 + Nbpts_D*2 + Nbpts

            param_real[ pt_coef:pt_coef + Nbpts_InterpD   ] = InterpD[1][0:Nbpts_InterpD]
            ### supression numpy initial
            InterpD[1]  =  param_real[ pt_coef:pt_coef + Nbpts_InterpD ]
            
            size_coef = size_coef + Nbpts_InterpD

            sname = rac[0][0:2]
            if (sname == 'IB'): 
                if (utau is not None):
                   var_ibc=['CoordinateX_PC','CoordinateY_PC','CoordinateZ_PC','CoordinateX_PI','CoordinateY_PI','CoordinateZ_PI','CoordinateX_PW','CoordinateY_PW','CoordinateZ_PW', 'Density','Pressure','utau','yplus']
                else:
                   var_ibc=['CoordinateX_PC','CoordinateY_PC','CoordinateZ_PC','CoordinateX_PI','CoordinateY_PI','CoordinateZ_PI','CoordinateX_PW','CoordinateY_PW','CoordinateZ_PW', 'Density','Pressure']

                count_ibc = 0
                ideb      = pt_coef + Nbpts_InterpD
                for v_ibc in var_ibc:
                   tmp                            = Internal.getNodeFromName1(rac, v_ibc)
                   param_real[ ideb:ideb+ Nbpts_D]= tmp[1][0:Nbpts_D]
                   tmp[1]                         = param_real[ ideb:ideb+ Nbpts_D ]
                   ideb                           = ideb + Nbpts_D
                   count_ibc   = count_ibc +1

                size_coef = size_coef + count_ibc*Nbpts_D

            c = c+1

#===============================================================================
# General transfers: Chimera + IBC - inplace version optimiser par arbre tc compacte par zone donneuse
# Interpolation is applied to aR 
# Beware: variables must be defined in topTreeD at nodes in order to be 
# consistent with the computation of connectivity by setInterpData and 
# setIBCData 
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variablesI =['var1','var2',...]: variables to be used in Chimera transfers
#                = None: the whole FlowSolutionNodes variables in topTreeD are transferred 
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers 
# IN: bcType (IBC only) 0: glissement
#                       1: adherence
#                       2: loi de paroi log
#                       3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE)+ronultideSA
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t)+nultideSA
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p)+nultideSA
# IN: storage=-1/0/1: unknown/direct/inverse
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
#===============================================================================
def ___setInterpTransfers(aR, topTreeD, 
                         variables=[], 
                         variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'], 
                         bcType=0, varType=1, storage=-1, compact=0,
                         Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08, 
                         Cs=0.3831337844872463, Ts=1.0):

    # Recup des donnees a partir des zones receveuses    
    if (storage != 1 or compact==0):
        print 'erreur __setInterpTransfers: Mode receveur a coder. Mode compact obligatoire: compact=1'
                            
    # Recup des donnees a partir des zones donneuses
    zones     = Internal.getZones(aR)
    zonesD    = Internal.getZones(topTreeD)
    param_int = Internal.getNodeFromName2(topTreeD, 'Parameter_int' )[1]
    param_real= Internal.getNodeFromName2(topTreeD, 'Parameter_real')[1]

    connector.___setInterpTransfers(zones, zonesD, variables , param_int, param_real, varType, bcType, Gamma, Cv, MuS, Cs, Ts )

    return None
    

#===============================================================================
# IN : subregion sr
#def _setIBCTransfers__(zsr, zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
#                       bcType, loc, varType, Gamma, Cv, MuS, Cs, Ts, extract):
#    xPC = Internal.getNodeFromName1(zsr,'CoordinateX_PC')[1]
#    yPC = Internal.getNodeFromName1(zsr,'CoordinateY_PC')[1]
#    zPC = Internal.getNodeFromName1(zsr,'CoordinateZ_PC')[1]
#    xPW = Internal.getNodeFromName1(zsr,'CoordinateX_PW')[1]
#    yPW = Internal.getNodeFromName1(zsr,'CoordinateY_PW')[1]
#    zPW = Internal.getNodeFromName1(zsr,'CoordinateZ_PW')[1]
#    xPI = Internal.getNodeFromName1(zsr,'CoordinateX_PI')[1]
#    yPI = Internal.getNodeFromName1(zsr,'CoordinateY_PI')[1]
#    zPI = Internal.getNodeFromName1(zsr,'CoordinateZ_PI')[1]
#
#    # Creation des numpy d extraction
#    pressNP = None; utauNP = None; yplusNP = None; densNP = None
#    if extract == 1:
#        nIBC = xPC.shape[0]
#        pressNode = Internal.getNodeFromName1(zsr,__PRESSURE__)
#        if pressNode is None: pressNP = numpy.zeros((nIBC),numpy.float64)
#        else: pressNP = Internal.getValue(pressNode)
#        densNode = Internal.getNodeFromName1(zsr,__DENSITY__)
#        if densNode is None: densNP = numpy.zeros((nIBC),numpy.float64)
#        else: densNP = Internal.getValue(densNode)
#
#        if bcType > 1:
#            utauNode = Internal.getNodeFromName1(zsr,__UTAU__)
#            if utauNode is None: utauNP = numpy.zeros((nIBC),numpy.float64)
#            else: utauNP = Internal.getValue(utauNode)
#            yplusNode = Internal.getNodeFromName1(zsr,__YPLUS__)
#            if yplusNode is None: yplusNP = numpy.zeros((nIBC),numpy.float64)
#            else: yplusNP = Internal.getValue(yplusNode)
#
#    # Transferts
#    connector._setIBCTransfers(zr, zd, variablesIBC, ListRcv, ListDonor,
#                               DonorType, Coefs,
#                               xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
#                               pressNP, utauNP, yplusNP, densNP,
#                               bcType, loc, varType, Gamma, Cv, MuS, Cs, Ts)
#
#    if extract == 1:
#        Internal._createUniqueChild(zsr,__PRESSURE__,'DataArray_t',value=pressNP)
#        Internal._createUniqueChild(zsr,__DENSITY__,'DataArray_t',value=densNP)
#
#        if bcType > 1:
#            Internal._createUniqueChild(zsr,__UTAU__, 'DataArray_t', value=utauNP)
#            Internal._createUniqueChild(zsr,__YPLUS__, 'DataArray_t', value=yplusNP)
#
#    return None

#===============================================================================
# General transfers: Chimera + IBC parallele - getFromZone
# Beware: variables must be defined in topTreeD at nodes in order to be consistent with
# the computation of connectivity by setInterpData and setIBCData
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variables=['var1','var2',...]: variables to be used in Chimera transfers
#              =[]: the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']:
# IN: bcType  0: glissement    (IBC only) p
#             1: adherence
#             2: loi de paroi log
#             3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE(,ronultideSA))
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t(,nutildeSA))
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p(,nutildeSA))
# Adim: KCore.adim1 for Minf=0.1 (IBC only)
#===============================================================================
def setInterpTransfersD(topTreeD, variables=[],
                        variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                        bcType=0, varType=1, compact=0,
                        Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                        Cs=0.3831337844872463, Ts=1.0, extract=0):

    # Recup des donnees a partir des zones donneuses
    zonesD = Internal.getZones(topTreeD)
    infos = []
    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:

          sname = s[0][0:2]
          #test pour eviter parcours arbre inutile 
          if ( (sname=='ID' or sname=='AD') and variables is not None) or (sname == 'IB' and variablesIBC is not None):
            dname = Internal.getValue(s) 
            idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
            if idn is not None: # la subRegion decrit des interpolations
                zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                zoneRole = Internal.getValue(zoneRole)
                if zoneRole == 'Donor':
                    location = Internal.getNodeFromName1(s, 'GridLocation') # localisation des donnees des receveurs
                    if location is not None:
                        location = Internal.getValue(location)
                    if location == 'CellCenter': loc = 'centers'
                    else: loc = 'nodes'
                    Coefs     = idn[1]
                    DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                    ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                    ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                    if (sname=='ID' or sname=='AD'):
                        #print 'transfert ID: zd ', zd[0]
                        arrayT = connector._setInterpTransfersD(zd, variables, ListDonor, DonorType, Coefs, varType, compact,                                                                
                                                                Internal.__GridCoordinates__, 
                                                                Internal.__FlowSolutionNodes__, 
                                                                Internal.__FlowSolutionCenters__)    
                        infos.append([dname,arrayT,ListRcv,loc])

                    elif sname == 'IB':
                        xPC = Internal.getNodeFromName1(s,'CoordinateX_PC')[1]
                        yPC = Internal.getNodeFromName1(s,'CoordinateY_PC')[1]
                        zPC = Internal.getNodeFromName1(s,'CoordinateZ_PC')[1]
                        xPW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                        yPW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                        zPW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                        xPI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                        yPI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                        zPI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                        Density = Internal.getNodeFromName1(s,'Density')[1]
                        Pressure= Internal.getNodeFromName1(s,'Pressure')[1]
                        utau    = Internal.getNodeFromName1(s,'utau')[1]
                        yplus   = Internal.getNodeFromName1(s,'yplus')[1]

                        #print 'transfert IBC : zd ', zd[0]
                        arrayT = connector._setIBCTransfersD(zd, variablesIBC, ListDonor, DonorType, Coefs, 
                                                             xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, Density, Pressure, utau, yplus,
                                                             bcType, varType, compact, Gamma, Cv, MuS, Cs, Ts,                                                             
                                                             Internal.__GridCoordinates__, 
                                                             Internal.__FlowSolutionNodes__, 
                                                             Internal.__FlowSolutionCenters__)
                        infos.append([dname,arrayT,ListRcv,loc])
    # Sortie
    return infos

#=============================================================================
# Calcul des coefficients d'interpolation et stockage dans l'arbre CGNS/Python
# ----------------------------------------------------------------------------
# setSeqInterpolation: calcul en sequentiel (stockage direct ou inverse
# setDistInterpolation: calcul dans un environnement distribue (stockage de type inverse)
# -------------------
# le celln doit etre defini dans l arbre et valoir :
# 0 pour les pts masques, 2 pour les pts interpoles
# -------------------
# IN: t: arbre sur lequel les interpolations sont calculees
# IN: loc: localisation des interpolations ('cell':cellules ou 'face':faces)
# IN: double_wall: activation de la technique double wall
# IN: storage: type de stockage (direct: sur le bloc interpole, inverse: sur le bloc d'interpolation)
# IN: prefixFile: prefixe pour le nom du fichier elsA
# IN: sameBase: si sameBase=1, les zones d'une meme base peuvent etre des zones d'interpolation
# IN: parallelDatas: liste de donnees a fournir dans un contexte parallele (=[] en sequentiel)
#     la liste contient: la liste de voisinage, le rang et la liste des cellules ou points EX interpoles.
# -------------------
# prefixFile == '': stockage dans l'arbre CGNS/Python uniquement
# prefixFile != '': stockage dans l'arbre CGNS/Python et ecriture dans un
#                   fichier relisible par elsA
#==============================================================================
def setInterpolations(t, loc='cell', double_wall=0,
                      storage='direct', prefixFile='',
                      sameBase=0, solver='elsA', nGhostCells=2, parallelDatas=[], check=True):

    # Solveur :
    if solver == 'elsA': Solver = 1
    elif solver == 'Cassiopee': Solver = 2

    # Localisation
    if loc =='cell': depth=2
    elif loc =='face': depth=1
    else: raise TypeError("setInterpolations: bad value for attribute 'loc'.")
    if parallelDatas == []: # mode sequentiel
        t = setSeqInterpolations(t, depth, double_wall, storage, prefixFile, sameBase, Solver, nGhostCells, check)
    else: # mode distribue
        t = setDistInterpolations(t, parallelDatas, depth, double_wall, sameBase, Solver, check)
    return t

#==============================================================================
# Calcul des coefficients d'interpolation et stockage dans l'arbre CGNS/Python
# dans un environnement distribue
#==============================================================================
def setDistInterpolations(t, parallelDatas=[], depth=2, double_wall=0, sameBase=0, solver=1, check=True):
    if double_wall == 1: import DoubleWall
    a = Internal.copyRef(t)

    if parallelDatas == []: return a
    else:
        if len(parallelDatas) != 3 :
            raise TypeError("setDistInterpolations: missing datas in parallelDatas.")
        else:
            graph = parallelDatas[0]
            rank = parallelDatas[1]
            interpDatas = parallelDatas[2]

    localGraph = graph[rank]
    for oppNode in localGraph.keys():
        oppZones = graph[oppNode][rank]
        for s in xrange(len(interpDatas[oppNode])):
            if depth == 2:
                interpCells=interpDatas[oppNode][s]
            else:
                EXPts=interpDatas[oppNode][s]
            oppZone = oppZones[s]
            # ----------------------------------------
            # Initialisation
            # ----------------------------------------
            # creation et initialisation a 1 du champ cellN s'il n'est pas deja present dans l'arbre
            a = addCellN__(a); a = addCellN__(a, loc='centers')
            bases = Internal.getBases(a); nbases = len(bases)

            # initialisation de liste de parois et de surfaces pour le traitement 'double_wall'
            firstWallCenters = []; surfacesExtC = []; walls1 = []

            # ----------------------------------------------------------------------------------------
            # Calcul des cellules et coefficients d'interpolation et ecriture dans un arbre CGNS/Python
            # ----------------------------------------------------------------------------------------
            # Calcul des donnees sur les domaines d interpolations (coordonnees, noms des zones, cellN, surfaces)
            if double_wall == 1:
                dwInfo = DoubleWall.extractDoubleWallInfo__(a)
                firstWallCenters = dwInfo[0]; surfacesExtC = dwInfo[1]
            interpolationZones=[]; interpolationZonesName=[]; interpolationCelln=[]; surfi=[]
            for nob in xrange(nbases):
                zones = Internal.getNodesFromType1(bases[nob], 'Zone_t')
                nzones = len(zones)
                for noz in xrange(nzones):
                    z = zones[noz]
                    zn = C.getFields(Internal.__GridCoordinates__,z)[0]
                    cellN = C.getField('centers:cellN', z)[0]
                    interpolationZonesName.append(z[0])
                    interpolationZones.append(zn)
                    interpolationCelln.append(cellN)
                    if (surfacesExtC != []): surfi.append(surfacesExtC[nob][noz])

            nBlksI = len(interpolationZones)

            # recuperation des domaines d interpolations
            # Calcul des cellules d interpolations et des coefficients associes
            # -----------------------------------------------------------------
            if nBlksI > 0:
                if depth == 2:
                    listOfInterpCells = []
                    # traitement double_wall
                    for i in xrange(nBlksI):
                        # changewall des centres et EX si paroi existe pour le domaine d interpolation
                        if surfi ==[] or walls1 == []:# attention: walls1 est tjs [] ici, pourquoi ?
                            listOfInterpCells.append(interpCells)
                        else:
                            if surfi[i] == []:
                                listOfInterpCells.append(interpCells)
                            else: # surface existe
                                ### DBG TO DO ...
                                #zc2 = Connector.changeWall__(zc, firstWallCenters1, surfi[i])
                                #interpCells = Connector.getInterpolatedPoints__(zc2)
                                #del zc2
                                listOfInterpCells.append(interpCells)
                    resInterp = Connector.setInterpolations__(oppZone, 1, 1, listOfInterpCells,
                                                              interpolationZones, interpolationCelln, isEX=0, check=check)
                    del listOfInterpCells
                # traitement des pts EX
                elif (depth == 1):
                    listOfEXPts = [];
                    for i in xrange(nBlksI):
                        if surfi == [] or walls1 == []:
                            listOfEXPts.append(EXPts);
                        else:
                            if surfi[i] == []:
                                listOfEXPts.append(EXPts)
                            else:
                                ### DBG TO DO ...
                                #expts = Connector.changeWallEX__(EXPts, zc, zn, walls1, surfi[i],\
                                #                                 doubleWallTol)
                                #listOfEXPts.append(expts);
                                listOfEXPts.append(EXPts);
                    del zn
                    resInterp = Connector.setInterpolations__(oppZone, 1, 1, listOfEXPts, interpolationZones, \
                                                              interpolationCelln, isEX=1, check=check)
                    del listOfEXPts

                for nozd in xrange(nBlksI):
                    zdonorname = interpolationZonesName[nozd]
                    zdonor = Internal.getNodesFromName(a, zdonorname)[0]
                    interpInverseStorage(oppZone,zdonor,nozd,resInterp,depth)
    return a

#=============================================================================
# Calcul des coefficients d'interpolation et stockage dans l'arbre CGNS/Python
# dans un environnement sequentiel
#==============================================================================
def setSeqInterpolations(t, depth=2, double_wall=0, storage='direct',
                         prefixFile='', sameBase=0, solver=1, nGhostCells=2, 
                         check=True):
    import Generator.PyTree as G
    if double_wall == 1: import DoubleWall

    a = Internal.copyRef(t)
    if storage != 'direct' and storage != 'inverse':
        raise TypeError("setInterpolations: bad value for attribute 'storage'.")
    # Ecriture ou non dans des fichiers elsA
    writeInFile=0
    if prefixFile != '': writeInFile=1
    # ----------------------------------------
    # Initialisation
    # ----------------------------------------
    # creation et initialisation a 1 du champ cellN s'il n'est pas deja present dans l'arbre
    a = addCellN__(a, loc='centers') # aux centres
    bases = Internal.getBases(a); nbases = len(bases)

    # Correspondance entre une zone et son Id globale pour l'ecriture dans les fichiers elsA
    ZonesId={}; Id = 0
    if writeInFile == 1:
        # declaration de dictionnaire pour le stockage de donnees par bloc donneur pour l'ecriture en fichier elsA
        listRcvId = {}; nbInterpCellsForDonor = {}; listIndicesRcv = {}; listEXdir={}; listIndicesDonor = {}; listInterpolantsDonor = {}; listCellN={}; listZoneDim={}; listInterpTypes = {}
    for nob in xrange(nbases):
        zones = Internal.getNodesFromType1(bases[nob], 'Zone_t')
        for noz in xrange(len(zones)):
            ZonesId[zones[noz][0]]=Id; Id = Id+1
    ntotZones = len(ZonesId)

    #=======================================================================
    # 2 - Recherche des periodicites :
    #     duplication des blocs periodiques dans les bases associees
    #     creation d un noeud fils au niveau de la zone dupliquee de nom 'TemporaryPeriodicZone'
    #=======================================================================
    for nob in xrange(nbases):
        parentb,db = Internal.getParentOfNode(a,bases[nob])
        bases[nob] = C.duplicatePeriodicZones__(bases[nob])
        parentb[2][db] = bases[nob]

    # initialisation de liste de parois et de surfaces pour le traitement 'double_wall'
    wallBndIndicesN=[]; surfacesExtC=[]; surfacesN=[]

    # -------------------------------------------------------------------------
    # Calcul des cellules et coefficients d interpolation et ecriture dans un
    # arbre CGNS/Python
    # -------------------------------------------------------------------------
    if writeInFile == 1:
        nbInterpCellsForDonor = [0]*ntotZones
    if double_wall == 1:
        dwInfo = DoubleWall.extractDoubleWallInfo__(a)
        firstWallCenters = dwInfo[0]; surfacesExtC = dwInfo[1]

    intersectionDict = getIntersectingDomains(a, method='AABB')

    for nob in xrange(nbases):
        if (bases[nob][0] == 'CARTESIAN'): sameBase = 1 # peut etre a enlever ??
        zones = Internal.getNodesFromType1(bases[nob], 'Zone_t')
        nzones = len(zones)
        for noz in xrange(nzones):
            z = zones[noz]
            zname = z[0]
            if Internal.getNodeFromName(z,'TempPeriodicZone') is not None: break

            # calcul d'arrays contenant les coefficients d'interpolation et les indices des cellules donneuses
            zn = C.getFields(Internal.__GridCoordinates__, z)[0]
            nirc = Internal.getZoneDim(z)[1]-1; njrc = Internal.getZoneDim(z)[2]-1; nkrc = Internal.getZoneDim(z)[3]-1
            if writeInFile == 1: listZoneDim[ZonesId[z[0]]]=[nirc,njrc,nkrc]
            # calcul des cellules d interpolations et des coefficients d interpolations
            zc = Converter.node2Center(zn)
            cellN = C.getField('centers:cellN',z)[0] # celln aux centres des cellules
            if writeInFile == 1: listCellN[ZonesId[z[0]]]=cellN[1]
            zc = Converter.addVars([zc,cellN]); del cellN
            firstWallCenters1 = []
            if double_wall == 1: firstWallCenters1 = firstWallCenters[nob][noz]

            # traitement des centres
            interpCells = Connector.getInterpolatedPoints__(zc)
            interpCells = Converter.extractVars(interpCells,['CoordinateX','CoordinateY','CoordinateZ','indcell'])
            # traitement des pts EX
            if depth == 1: EXPts = Connector.getEXPoints__(zn, zc)
            if (interpCells[1].size != 0): # pts interpoles
                # recuperation des domaines d interpolations
                res = getInterpolationDomains__(bases,nob,noz,zn,surfacest=surfacesExtC,sameBase=sameBase, intersectionDict=intersectionDict[zname])
                interpolationZonesName = res[0]; interpolationZones = res[1]; interpolationCelln = res[2]; surfi = res[3]; periodicZones = res[4]
                nBlksI = len(interpolationZones)
                # Calcul des cellules d interpolations et des coefficients associes
                # -------------------------------------------------------------
                if nBlksI > 0:
                    if depth == 2:
                        listOfInterpCells = []
                        # traitement double_wall
                        for i in xrange(nBlksI):
                            # changewall des centres et EX si paroi existe pour le domaine d interpolation
                            if (firstWallCenters1 == [] or surfi == []):
                                listOfInterpCells.append(interpCells)
                            else:
                                if (surfi[i] == []):
                                    listOfInterpCells.append(interpCells)
                                else: # surface existe
                                    zc2 = Connector.changeWall__(zc, firstWallCenters1, surfi[i])
                                    interpCells = Connector.getInterpolatedPoints__(zc2)
                                    interpCells = Converter.extractVars(interpCells,['CoordinateX','CoordinateY','CoordinateZ','indcell'])
                                    del zc2
                                    listOfInterpCells.append(interpCells)
                        resInterp = Connector.setInterpolations__(z[0],nirc, njrc, listOfInterpCells, interpolationZones, \
                                                                  interpolationCelln, isEX=0, zoneId=ZonesId[z[0]]+1, check=check)
                        del listOfInterpCells

                    # traitement des pts EX
                    elif depth == 1:
                        listOfEXPts = []
                        for i in xrange(nBlksI):
                            zdonorId = ZonesId[interpolationZonesName[i]]
                            if (surfi == [] or firstWallCenters1 == []): listOfEXPts.append(EXPts)
                            else:
                                if surfi[i] == []: listOfEXPts.append(EXPts)
                                else:
                                    expts = Connector.changeWallEX__(EXPts, zc, zn, firstWallCenters1, surfi[i])
                                    listOfEXPts.append(expts)
                        del zn
                        resInterp = Connector.setInterpolations__(z[0],nirc, njrc, listOfEXPts, interpolationZones, \
                                                                  interpolationCelln, isEX=1, zoneId=ZonesId[z[0]]+1, check=check)
                        del listOfEXPts
                    for nozd in xrange(nBlksI):
                        indicesRcv = resInterp[0][nozd]
                        indicesDonor=resInterp[1][nozd]
                        indicesDonorExtC = resInterp[2][nozd]
                        interpCoef=resInterp[3][nozd];
                        interpVol=resInterp[4][nozd]; interpType=resInterp[5][nozd]
                        indicesExtrap = resInterp[6][nozd]; indicesOrphan = resInterp[7][nozd]
                        if depth == 1: EXdir = resInterp[8][nozd]
                        zdonorname = interpolationZonesName[nozd]
                        isperiodic = periodicZones[nozd]
                        if isperiodic == 2:
                            print 'Periodic interpolation from +theta: ', interpType.shape[0]
                            interpType = 102*numpy.ones((interpType.shape[0]),numpy.int32)
                        elif isperiodic == 3:
                            print 'Periodic interpolation from -theta: ', interpType.shape[0]
                            interpType = 103*numpy.ones((interpType.shape[0]),numpy.int32)
                        resInterp[5][nozd] = interpType
                        zdonor = Internal.getNodesFromName(a,zdonorname)[0]
                        # ----------------------------------------
                        # stockage direct (sur le bloc interpole)
                        # ----------------------------------------
                        if storage == 'direct':
                            interpDirectStorage(z,zdonorname,nozd,resInterp,depth)
                        # -------------------------------------
                        # stockage inverse (sur le bloc donneur)
                        # -------------------------------------
                        else:
                            interpInverseStorage(z[0],zdonor,nozd,resInterp,depth)
                        # -----------------------------------------------------
                        # Stockage de donnees pour l ecriture des fichiers elsA
                        # -----------------------------------------------------
                        # - nombre total de points interpoles par le bloc donneur zdonorname
                        # - liste des donnees a ecrire par bloc donneur (indices cell. interpolees, indices cell. donneuses, coefficients d interpolation)
                        if (writeInFile == 1):
                            donorId = ZonesId[zdonorname]
                            # dimensions de la grille d interpolation
                            ni = Internal.getZoneDim(zdonor)[1]; nj = Internal.getZoneDim(zdonor)[2]; nk = Internal.getZoneDim(zdonor)[3]
                            nbInterpCellsForDonor[donorId] = nbInterpCellsForDonor[donorId] + len(indicesDonor)
                            if donorId in listRcvId:
                                listRcvId[donorId].append(ZonesId[z[0]])
                                listIndicesRcv[donorId].append(indicesRcv)
                                listIndicesDonor[donorId].append(indicesDonorExtC)
                                listInterpolantsDonor[donorId].append(interpCoef)
                                listInterpTypes[donorId].append(interpType)
                            else:
                                listRcvId[donorId]=[ZonesId[z[0]]]
                                listIndicesRcv[donorId]=[indicesRcv]
                                listIndicesDonor[donorId]=[indicesDonorExtC]
                                listInterpolantsDonor[donorId]=[interpCoef]
                                listInterpTypes[donorId]=[interpType]

                            if depth ==1:
                                if donorId in listEXdir: listEXdir[donorId].append(EXdir)
                                else: listEXdir[donorId]=[EXdir]

                else: # no interpolation domain found
                    if depth == 2: indicesOrphan = numpy.array(interpCells[1][3],dtype='int32')
                    elif depth == 1: indicesOrphan = numpy.array(EXPts[1][3],dtype='int32')
                    # on cree une zone subregion avec les pts orphelins
                    nameSubRegion='Orphan_'+z[0]
                    z[2].append([nameSubRegion, None, [],'ZoneSubRegion_t'])
                    info = z[2][len(z[2])-1]
                    v = numpy.fromstring('Receiver', 'c'); info[2].append(['ZoneRole',v , [], 'DataArray_t'])
                    if depth == 1: info[2].append(['OrphanFaceList',numpy.unique(indicesOrphan) , [], 'DataArray_t'])
                    else: info[2].append(['OrphanPointList',numpy.unique(indicesOrphan) , [], 'DataArray_t'])
    # Delete duplicated periodic zones
    a = C.removeDuplicatedPeriodicZones__(a)

    # -------------------------------
    # Ecriture dans un fichier "elsA"
    # -------------------------------
    # Parcours des zones. On regarde si ce sont des zones donneuses. Si oui, on ecrit le(s) fichier(s) elsA correspondant(s).
    if writeInFile == 1:
        if depth == 2: EX=0
        elif depth == 1: EX=1
        Connector.writeCoefs(ntotZones, listRcvId,listIndicesRcv,listEXdir,listIndicesDonor,listInterpolantsDonor,
                             listInterpTypes,listCellN,listZoneDim,nbInterpCellsForDonor,
                             prefixFile, isEX=EX, solver=solver, nGhostCells=nGhostCells)
    return a

#------------------------------------------------------------------------------
# Calcul des transferts Chimere
# IN: vars: liste des noms des champs sur lesquels le transfert est applique
#------------------------------------------------------------------------------
def chimeraTransfer(t, storage='direct', variables=[], loc='cell',mesh='standard'):
    vars2 = []
    for v in variables:
        v2 = v.split(':')
        if len(v2) == 2 and v2[0] == 'centers': vars2.append(v)
        # DBG else: print 'Warning: chimeraTransfer: only variables located at centers taken into account.'
    if vars2 == []:
        raise ValueError("chimeraTransfer: no variable to transfer.")
    if loc != 'cell': raise ValueError("chimeraTransfer: only valid for cell centers.")
    if storage == 'direct': return directChimeraTransfer(t, vars2, loc,mesh)
    else: return inverseChimeraTransfer(t, vars2, loc,mesh)

#------------------------------------------------------------------------------
# Calcul des transferts Chimere a partir d'un stockage direct
# Les champs sur lesquels on applique les transferts Chimere sont
# entres par l'utilisateur
# IN: vars: liste des noms des champs sur lesquels le transfert est applique
# locinterp = 'cell' ou 'face'
# Pour l instant on suppose que les interpolations sont realisees aux centres
# donc les transferts aussi
#------------------------------------------------------------------------------
def directChimeraTransfer(t, variables, locinterp,mesh='standard'):
    loc = 'centers'; depth = 2
    if locinterp == 'face': depth = 1

    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp); nzones = len(zones)
    for noz in xrange(nzones):
        z = zones[noz]
        parent1,d1 = Internal.getParentOfNode(tp,z)
        lRcvArrays = []
        for i in variables:
            f = C.getField(i,z)[0]
            if f != []: lRcvArrays.append(f)
        if lRcvArrays != []: rcvArray = Converter.addVars(lRcvArrays)
        # Les interpolations sont stockees dans les noeuds ZoneSubRegion_t
        subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        subRegions=[]
        for s in subRegions2:
            sname = s[0]
            sname0 = sname.split('_')[0]
            if sname0=='ID' or sname0=='AD': subRegions.append(s)
        for s in subRegions:
            idn = Internal.getNodeFromName1(s,'InterpolantsDonor')
            if idn is not None: # la subRegion decrit des interpolations Chimere
                interpDonor = idn[1]
                if depth == 1:
                    cellDonorEX = Internal.getNodeFromName1(s,'FaceListDonor')[1]
                    interpDonorEX = Internal.getNodeFromName1(s,'FaceInterpolantsDonor')[1]
                    interpDonorEXType = Internal.getNodeFromName1(s,'FaceInterpolantsType')[1]
                interpDonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                cellRcv   = Internal.getNodeFromName1(s,'PointList')[1]
                cellDonor = Internal.getNodeFromName1(s,'PointListDonor')[1]
                zdonorname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(tp, zdonorname)
                zdonor = Internal.getNodesFromType1(zdonors,'Zone_t')[0]
                lDonorArrays = []
                for i in variables:
                    f = C.getField(i,zdonor)[0]
                    if f !=[]: lDonorArrays.append(f)
                if lDonorArrays != []: donorArray = Converter.addVars(lDonorArrays)
                # dimensions de la grille d interpolation
                ni = Internal.getZoneDim(zdonor)[1]; nj = Internal.getZoneDim(zdonor)[2]; nk = Internal.getZoneDim(zdonor)[3]
                # dimensions de la grille interpolee
                nir = Internal.getZoneDim(z)[1]; njr = Internal.getZoneDim(z)[2]
                if loc == 'centers':
                    ni = ni-1; nj = nj-1; nk = nk-1; nir = nir-1; njr = njr-1
                rcvArray = Connector.chimeraTransfer(cellRcv, cellDonor, interpDonorType, interpDonor,
                                                     rcvArray, donorArray)
##                 # traitement des pts EX
##                 if depth == 1:
##                     rcvArray=Connector.chimeraTransfer(cellRcv, cellDonorEX, interpDonorEXType,interpDonorEX,
##                rcvArray,donorArray)
                C.setFields([rcvArray], z, loc)
                parent1[2][d1] = z
    return tp

#------------------------------------------------------------------------------
# Calcul des transferts Chimere a partir d'un stockage inverse
# Les noms des champs sur lesquels on applique les transferts Chimere sont
# entres par l'utilisateur
# IN: variables: liste des noms des champs sur lesquels le transfert est applique
# IN: locinterp='cell' ou 'face', ne marche que pour cell actuellement
# IN: mesh='standard' or 'extended'
#------------------------------------------------------------------------------
def inverseChimeraTransfer(t, variables, locinterp= 'centers',mesh='standard'):
    depth = 2 # valeur par defaut pour des interpolations 'centers'
    if locinterp == 'face': depth = 1

    tp = Internal.copyRef(t)
    # On regarde si t est une zone ou une liste de zones
    isZoneOrListZone = False
    if Internal.isStdNode(tp) == -1: # noeud
        if tp[3] == 'Zone_t': isZoneOrListZone = True; tp = [tp]
    elif Internal.isStdNode(tp) == 0: #liste
        for zt in tp:
            if zt[3] == 'Zone_t': isZoneOrListZone = True
            else: isZoneOrListZone = False; break

    if isZoneOrListZone: zones = tp
    else: zones = Internal.getZones(tp)

    # liste des zones locales (sur le processeur courant)
    localZonesName = []
    for z in zones: localZonesName.append(z[0])
     
    # Definition de dictionnaires pour les champs interpoles : 
    # Cle : nom zone receveuse. Valeurs : champs interpole + indices pts receveurs + volume cellule donneuse    
    rcvFields = {} 

    if depth == 1: rcvLocalFields = {}

    # Parcours des zones locales (qui sont zones donneuses)
    nzones = len(zones)
    for noz in xrange(nzones):
        zdonor = zones[noz]
        # Recuperation des champs des zones donneuses pour les grandeurs a transferer
        # Ces champs sont stockes dans un array unique donorArray
        lDonorArrays = []
        for i in variables:
            f = C.getField(i,zdonor)[0]
            if f !=[]: 
                # Passage du maillage donneur en centres etendus
                if mesh == 'extended': lDonorArrays.append(Converter.center2ExtCenter(f))
                else: lDonorArrays.append(f)
        if lDonorArrays != []: donorArray = Converter.addVars(lDonorArrays)
       # Parcours des ZoneSubRegion pour effectuer les transferts
       # (Les interpolations sont stockees dans les noeuds ZoneSubRegion_t)
        subRegions = Internal.getNodesFromType1(zdonor,'ZoneSubRegion_t')
        for s in subRegions:
            if s[0].split('_')[0] == 'ID': # la ZoneSubRegion decrit des interpolations Chimere
                idn = Internal.getNodeFromName1(s,'InterpolantsDonor')
                idEX = Internal.getNodeFromName1(s,'FaceInterpolantsDonor')
                if idn is not None or idEX is not None: # la ZoneSubRegion decrit des interpolations Chimere
                    rcvZoneName = Internal.getValue(s) # nom de la zone interpolee
                    #A.  traitement locinterp = centers
                    if idn is not None and depth == 2:
                        # A.1 recuperation des champs pour l interpolation
                        if mesh == 'extended': cellDonorList = Internal.getNodeFromName1(s,'PointListExtC')
                        else: cellDonorList = Internal.getNodeFromName1(s,'PointList')
                        if cellDonorList is not None and cellDonorList[1] != []: # la region traite des interpolations depth = 2
                            cellDonor = cellDonorList[1]
                            cellRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                            interpDonor = idn[1]
                            cellVol   = Internal.getNodeFromName1(s,'InterpolantsVol')[1]
                            interpDonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                            # A.2 creation du array pour stocker la solution interpolee
                            if rcvZoneName not in localZonesName: # zone receveuse pas sur le processeur locale
                                vars2 = ''
                                for v in variables:
                                    v2 = v.split(':')[1];
                                    if variables.index(v) == len(variables)-1: vars2=vars2+v2
                                    else: vars2=vars2+v2+','
                                cellRcvPara = numpy.arange(cellRcv.shape[0],dtype=numpy.int32)
                                rcvField = numpy.empty((donorArray[1].shape[0],cellRcv.shape[0])) # DBG 
                                # DBG rcvField = numpy.empty((donorArray[1].shape[0],cellRcv.shape[0]), order='Fortran') # DBG FORTRAN
                                rcvArray = [vars2,rcvField,cellRcv.shape[0],1,1] # tableau monodimensionnel
                            else: # zone receveuse sur le processeur locale
                                cellRcvPara = cellRcv
                                zs = Internal.getNodesFromName2(tp,rcvZoneName)
                                for zr in zs:
                                    if zr[3] == 'Zone_t':
                                        z = zr; break
                                lRcvArrays = []
                                for i in variables:
                                    f = C.getField(i,z)[0]
                                    if f !=[]: lRcvArrays.append(f)
                                if lRcvArrays != []: rcvArray = Converter.addVars(lRcvArrays)
                                else:
                                    raise ValueError("inverseChimeraTransfer: field to interpolate not in zone %d."%rcvZoneName)      
                            # A.3 interpolation Chimere
                            rcvField = Connector.chimeraTransfer(cellRcvPara, cellDonor, interpDonorType, 
                                                                  interpDonor, rcvArray, donorArray)[1]
                            if rcvZoneName not in localZonesName:
                                # redimensionnement du tableau au nombre de cellules interpolees
                                rcvField2 = numpy.empty((donorArray[1].shape[0]+2,cellRcv.shape[0]), order='Fortran')
                                rcvField2[:rcvField.shape[0]]=rcvField
                                rcvField2[donorArray[1].shape[0]]=cellRcv
                                rcvField2[donorArray[1].shape[0]+1]=cellVol
                                if rcvZoneName in rcvFields.keys():
                                    rcvFields[rcvZoneName].append(rcvField2)
                                else:
                                    rcvFields[rcvZoneName]=[rcvField2]
                            else: # solution mise directement dans l arbre local
                                rcvArray[1]=rcvField
                                C.setFields([rcvArray], z, locinterp)
                    # traitement locinterp = face
                    elif idEX != None and depth == 1:
                        if mesh == 'extended': cellDonorListEX = Internal.getNodeFromName1(s,'FacePointListExtC')
                        else : cellDonorListEX = Internal.getNodeFromName1(s,'FacePointList')
                        if cellDonorListEX is not None: # la region traite des interpolations depth = 1
                            cellDonorEX = cellDonorListEX[1]
                            if cellDonorEX != []:
                                cellRcvEX = Internal.getNodeFromName1(s,'FaceListDonor')[1]
                                # DBG rcvFieldEX = numpy.empty((donorArray[1].shape[0],cellRcvEX.shape[0]), order='Fortran')
                                EXdir = Internal.getNodeFromName1(s,'FaceDirection')[1]
                                interpDonorEX = idEX[1]
                                cellVolEX   = Internal.getNodeFromName1(s,'FaceInterpolantsVol')[1]
                                interpDonorEXType = Internal.getNodeFromName1(s,'FaceInterpolantsType')[1]
                                rcvArrayEX = Converter.array(donorArray[0],cellRcvEX.shape[0], 1, 1)
                                rcvFieldEX = Connector.chimeraTransfer(cellRcvEX,cellDonorEX,interpDonorEXType,
                                                                       interpDonorEX, rcvArrayEX, donorArray)[1]
                                rcvField2 = numpy.empty((donorArray[1].shape[0]+3,cellRcvEX.shape[0]), order='Fortran')
                                rcvField2[:rcvFieldEX.shape[0]]=rcvFieldEX
                                rcvField2[donorArray[1].shape[0]]=cellRcvEX
                                rcvField2[donorArray[1].shape[0]+1]=cellVolEX
                                rcvField2[donorArray[1].shape[0]+2]=EXdir
                                if rcvZoneName not in localZonesName:
                                    if rcvZoneName in rcvFields.keys():
                                        rcvFields[rcvZoneName].append(rcvField2)
                                    else:
                                        rcvFields[rcvZoneName]=[rcvField2]
                                else:
                                    if rcvZoneName in rcvLocalFields.keys():
                                        rcvLocalFields[rcvZoneName].append(rcvField2)
                                    else:
                                        rcvLocalFields[rcvZoneName]=[rcvField2]

    if depth==1: return rcvLocalFields, rcvFields
    return tp, rcvFields

#-----------------------------------------------------------------------------
# IN: zname: name of interpolated zone
# IN: zdonor: donor zone
# OUT: modified donor zone
#-----------------------------------------------------------------------------
def interpInverseStorage(zname, zdonor, nozd, resInterp, depth):
    if depth == 1: loc = 'faces'
    else: loc = 'centers'
    indicesDonor=resInterp[1][nozd]
    indicesDonorExtC=resInterp[2][nozd]
    indicesOrphan = resInterp[7][nozd]
    if indicesDonor.size > 0 or indicesOrphan.size > 0:
        indRcv = resInterp[0][nozd]; interpCoef=resInterp[3][nozd];
        interpVol=resInterp[4][nozd]; interpType=resInterp[5][nozd]
        indicesExtrap = resInterp[6][nozd];
        if depth == 2: EXdir=numpy.array([],numpy.int32)
        elif depth == 1: EXdir = resInterp[8][nozd]
        cellIndExtrap = indicesExtrap; cellIndOrphan = indicesOrphan
        coef=interpCoef; vol=interpVol; interptype=interpType              
        # creation des noeuds ZoneSubRegion_t contenant les coefficients d'interpolation et les indices des cellules donneuses
        # Il y a un noeud ZoneSubRegion_t par paire de blocs interpole / bloc donneur
        createInterpRegion__(zdonor,zname, indicesDonor, indRcv, coef,
                             interptype, vol, cellIndExtrap,cellIndOrphan,
                             EXdir,indicesDonorExtC, 'Donor',loc=loc)
    return

#-----------------------------------------------------------------------------
# IN: zdonorname: name of donor zone
# IN: z: interpolated zone
# OUT: modified interpolated zone
#-----------------------------------------------------------------------------
def interpDirectStorage(z, zdonorname, nozd, resInterp, depth):
    if depth == 1: loc = 'faces'
    else: loc = 'centers'
    indicesDonor=resInterp[1][nozd]
    indicesDonorExtC=resInterp[2][nozd]
    indicesOrphan = resInterp[7][nozd]
    if indicesDonor.size > 0 or indicesOrphan.size >0:
        indRcv = resInterp[0][nozd]; interpCoef=resInterp[3][nozd]
        interpVol=resInterp[4][nozd]; interpType=resInterp[5][nozd]
        indicesExtrap = resInterp[6][nozd];
        if depth == 2: EXdir=numpy.array([],numpy.int32)
        elif depth == 1: EXdir = resInterp[8][nozd]
        cellIndExtrap = indicesExtrap; cellIndOrphan = indicesOrphan
        coef = interpCoef; vol = interpVol; interptype = interpType
        # creation des noeuds ZoneSubRegion_t contenant les coefficients d'interpolation et
        # les indices des cellules donneuses
        # Il y a un noeud ZoneSubRegion_t par paire de blocs interpole / bloc donneur
        createInterpRegion__(z, zdonorname, indRcv, indicesDonor, coef,
                             interptype, vol, cellIndExtrap, cellIndOrphan,
                             EXdir,indicesDonorExtC, 'Receiver', loc=loc)
    return

#------------------------------------------------------------------------------
# Creation or completion (if it already exists) of a node ZoneSubRegion_t for interpolation datas
# IN/OUT: z zone to be modified
# IN: zname: nom de la zone (donneuse si direct ou receveuse si inverse)
# IN: pointlist: liste des points ds PointList (points donneurs si stockage inverse, receveurs sinon)
# IN: pointlistdonor: liste des pts (points donneurs si stockage direct, receveurs sinon)
# IN: interpCoef: coefs d interpolation
# IN: interpType: type d interpolation
# IN: indicesExtrap: indices des pts extrapoles
# IN: indicesOrphan: indices de ts les pts orphelins (dupliques pour ttes les zones donneuses)
# IN: tag:'Receiver' or 'Donor'
# IN: loc='centers', 'nodes' ou 'faces'
# IN: EXDir: direction pour les pts EX (Chimere depth=1)
# IN: itype: 'abutting' pour transferts multidomaine, 'chimera' pour les interp chimere, 
# 'ibc' pour les IBC
#------------------------------------------------------------------------------
def createInterpRegion__(z, zname, pointlist, pointlistdonor, interpCoef, interpType, interpVol, indicesExtrap, indicesOrphan, \
                             EXDir = [],pointlistdonorExtC = None, tag='Receiver', loc='centers', itype='chimera',
                             prefix=None, RotationAngle=None, RotationCenter=None):
    if prefix is None:
        if itype == 'chimera': nameSubRegion = 'ID_'+zname
        elif itype == 'abutting': nameSubRegion = 'ID_'+zname
        elif itype=='ibc': nameSubRegion='IBCD_'+zname
    else: nameSubRegion=prefix+'_'+zname

    zsr = Internal.getNodesFromName1(z,nameSubRegion)
    # create new subregion for interpolations
    if zsr == []:
        v = numpy.fromstring(zname, 'c')
        z[2].append([nameSubRegion, v, [],'ZoneSubRegion_t'])
        info = z[2][len(z[2])-1]
        v = numpy.fromstring(tag, 'c'); info[2].append(['ZoneRole',v , [], 'DataArray_t'])
        if loc == 'faces':
            info[2].append(['FaceList',      pointlist, [], 'IndexArray_t'])
            info[2].append(['FaceListDonor', pointlistdonor, [], 'IndexArray_t'])
            if pointlistdonorExtC is not None:
                if tag == 'Receiver':
                    info[2].append(['FaceListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                else:
                    info[2].append(['FaceListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
            info[2].append(['FaceDirection', EXDir, [], 'DataArray_t'])
            info[2].append(['FaceInterpolantsDonor',interpCoef , [], 'DataArray_t'])
            if interpVol.size > 0: info[2].append(['FaceInterpolantsVol',interpVol , [], 'DataArray_t'])
            info[2].append(['FaceInterpolantsType',interpType , [], 'DataArray_t'])
            if indicesOrphan.size>0: info[2].append(['OrphanFaceList',numpy.unique(indicesOrphan) , [], 'DataArray_t'])
            if indicesExtrap.size>0: info[2].append(['ExtrapFaceList',indicesExtrap , [], 'DataArray_t'])
        else:
            if loc == 'centers':
                v = numpy.fromstring('CellCenter', 'c')
                info[2].append(['GridLocation', v, [], 'GridLocation_t'])
            info[2].append(['PointList',      pointlist , [], 'IndexArray_t'])
            info[2].append(['PointListDonor', pointlistdonor , [], 'IndexArray_t'])
            if pointlistdonorExtC != None:
                if tag == 'Receiver':
                    info[2].append(['PointListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                else:
                    info[2].append(['PointListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
            info[2].append(['InterpolantsDonor',interpCoef , [], 'DataArray_t'])
            if interpVol.size>0: info[2].append(['InterpolantsVol',interpVol , [], 'DataArray_t'])
            info[2].append(['InterpolantsType',interpType , [], 'DataArray_t'])
            if indicesOrphan.size>0: info[2].append(['OrphanPointList',numpy.unique(indicesOrphan) , [], 'DataArray_t'])
            if indicesExtrap.size>0: info[2].append(['ExtrapPointList',indicesExtrap, [], 'DataArray_t'])
        if RotationAngle is not None: info[2].append(RotationAngle)
        if RotationCenter is not None: info[2].append(RotationCenter)
    else:
        if loc == 'faces':
            intd = Internal.getNodesFromName1(zsr[0],'FaceInterpolantsDonor')
            if intd == []:
                l = len(zsr[0][2]); info = zsr[0][2][l-1]
                zsr[0][2].append(['FaceList',      pointlist, [], 'IndexArray_t'])
                zsr[0][2].append(['FaceListDonor', pointlistdonor, [], 'IndexArray_t'])
                if pointlistdonorExtC != None:
                    if tag == 'Receiver':
                        zsr[0][2].append(['FaceListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                    else:
                        zsr[0][2].append(['FaceListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                zsr[0][2].append(['FaceInterpolantsDonor',interpCoef, [], 'DataArray_t'])
                zsr[0][2].append(['FaceDirection',EXDir, [], 'DataArray_t'])
                if interpVol.size>0: zsr[0][2].append(['FaceInterpolantsVol',interpVol, [], 'DataArray_t'])
                zsr[0][2].append(['FaceInterpolantsType',interpType, [], 'DataArray_t'])
                if indicesOrphan.size>0: zsr[0][2].append(['OrphanFaceList',numpy.unique(indicesOrphan),[],'DataArray_t'])
                if indicesExtrap.size>0: zsr[0][2].append(['ExtrapFaceList', indicesExtrap,[],'DataArray_t'])
                if RotationAngle is not None: zsr[0][2].append(RotationAngle)
                if RotationCenter is not None: zsr[0][2].append(RotationCenter)
            else:
                intd[0][1] = numpy.concatenate((intd[0][1],interpCoefs))
                intv = Internal.getNodesFromName1(zsr[0],'FaceInterpolantsVol')
                if intv!= [] and interpVol.size>0: intv[0][1]=numpy.concatenate((intv[0][1],interpVol))
                inttype = Internal.getNodesFromName1(zsr[0],'FaceInterpolantsType')
                if inttype != []: inttype[0][1] = numpy.concatenate((inttype[0][1],interpType))
                node = Internal.getNodesFromName1(zsr[0],'FaceList')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlist))
                node = Internal.getNodesFromName1(zsr[0],'FaceListDonor')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonor))
                node = Internal.getNodesFromName1(zsr[0],'FaceDirection')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],EXDir))
                node = Internal.getNodesFromName1(zsr[0],'OrphanFaceList')
                if node!=[]: node[0][1]=numpy.unique(indicesOrphan)
                node = Internal.getNodesFromName1(zsr[0],'ExtrapFaceList')
                if node!=[]: node[0][1]=numpy.concatenate((node[0][1],indicesExtrap))
        else:
            intd = Internal.getNodesFromName1(zsr[0],'InterpolantsDonor')
            if intd == []:
                l = len(zsr[0][2]); info = zsr[0][2][l-1]
                if loc == 'centers':
                    v = numpy.fromstring('CellCenter', 'c')
                    zsr[0][2].append(['GridLocation', v, [], 'GridLocation_t'])
                zsr[0][2].append(['PointList',      pointlist, [], 'IndexArray_t'])
                zsr[0][2].append(['PointListDonor', pointlistdonor, [], 'IndexArray_t'])
                if pointlistdonorExtC != None:
                    if tag == 'Receiver':
                        zsr[0][2].append(['PointListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                    else:
                        zsr[0][2].append(['PointListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                zsr[0][2].append(['InterpolantsDonor',interpCoef, [], 'DataArray_t'])
                if interpVol.size>0: zsr[0][2].append(['InterpolantsVol',interpVol, [], 'DataArray_t'])
                zsr[0][2].append(['InterpolantsType',interpType, [], 'DataArray_t'])
                if indicesOrphan.size>0: zsr[0][2].append(['OrphanPointList',numpy.unique(indicesOrphan),[], 'DataArray_t'])
                if indicesExtrap.size>0: zsr[0][2].append(['ExtrapPointList',indicesExtrap, [], 'DataArray_t'])
                if RotationAngle is not None: zsr[0][2].append(RotationAngle)
                if RotationCenter is not None: zsr[0][2].append(RotationCenter)
            else:
                intd[0][1] = numpy.concatenate((intd[0][1],interpCoef))
                intv = Internal.getNodesFromName1(zsr[0],'InterpolantsVol')
                if intv!= [] and interpVol.size>0: intv[0][1] = numpy.concatenate((intv[0][1],interpVol))
                inttype = Internal.getNodesFromName1(zsr[0],'InterpolantsType')
                if inttype != []: inttype[0][1] = numpy.concatenate((inttype[0][1],interpType))
                node = Internal.getNodesFromName1(zsr[0],'PointList')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlist))
                node = Internal.getNodesFromName1(zsr[0],'PointListDonor')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonor))
                node = Internal.getNodesFromName1(zsr[0],'PointListExtC')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonorExtC))
                node = Internal.getNodesFromName1(zsr[0],'PointListDonorExtC')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonorExtC))
                node = Internal.getNodesFromName1(zsr[0],'OrphanPointList')
                if node!=[]: node[0][1]=numpy.unique(indicesOrphan)
                node = Internal.getNodesFromName1(zsr[0],'ExtrapPointList')
                if node!=[]: node[0][1]=numpy.concatenate((node[0][1],indicesExtrap))

    return None

#=============================================================================
# Get overset info : if a cell is interpolated, orphan, extrapolated
#                    or the quality of the donor cell
# 'orphan' (0,1), 'extrapolated'(sum |cf|), 'interpolated'(0,1)
# 'cellRatio' : max(volD/volR,volR/volD)
# 'donorAspect' : edgeRatio of donor cell
# Compliant with setInterpolations storage
#=============================================================================
def chimeraInfo(a, type='interpolated'):
    try: import OversetInfo
    except: raise ImportError("chimeraInfo requires OversetInfo module.")
    t = Internal.copyRef(a)
    if type == 'interpolated':
        C._initVars(t,'{centers:interpolated}=0.')
        t = OversetInfo.chimeraNatureOfCells__(t,'interpolated')
    elif type == 'extrapolated':
        C._initVars(t,'{centers:extrapolated}=0.')
        t = OversetInfo.chimeraNatureOfCells__(t,'extrapolated')
    elif type =='orphan':
        C._initVars(t,'{centers:orphan}=0.')
        t = OversetInfo.chimeraNatureOfCells__(t,'orphan')
    elif type =='cellRatio':
        C._initVars(t,'{centers:cellRatio}=0.')
        t = OversetInfo.chimeraCellRatio__(t)
    elif type == 'donorAspect':
        C._initVars(t,'{centers:donorAspect}=0.')
        t = OversetInfo.chimeraDonorAspect__(t)
    else: raise ValueError("chimeraInfo: type is not valid")
    return t

#=============================================================================
# Get overset info: if a cell is interpolated, orphan, extrapolated
#                   or the quality of the donor cell
# 'orphan' (0,1), 'extrapolated'(sum |cf|), 'interpolated'(0,1)
# 'cellRatio': max(volD/volR,volR/volD)
# 'donorAspect': edgeRatio of donor cell
# Compliant with setInterpData storage
#=============================================================================
def getOversetInfo(aR, topTreeD, loc='nodes', type='interpolated'):
    try: import OversetInfo
    except: raise ImportError("getOversetInfo requires OversetInfo module.")
    tR = Internal.copyRef(aR)

    if type == 'interpolated':
        if loc == 'nodes': C._initVars(tR,'{interpolated}=0.')
        else: C._initVars(tR,'{centers:interpolated}=0.')
        tR = OversetInfo.oversetNatureOfCells__(tR,topTreeD,'interpolated')
    elif type == 'extrapolated':
        if loc == 'nodes': C._initVars(tR,'{extrapolated}=0.')
        else: C._initVars(tR,'{centers:extrapolated}=0.')
        tR = OversetInfo.oversetNatureOfCells__(tR,topTreeD,'extrapolated')
    elif type =='orphan':
        if loc == 'nodes': C._initVars(tR,'{orphan}=0.')
        else: C._initVars(tR,'{centers:orphan}=0.')
        tR = OversetInfo.oversetNatureOfCells__(tR,topTreeD,'orphan')
    elif type =='cellRatio':
        if loc == 'nodes': C._initVars(tR,'{cellRatio}=0.')
        else: C._initVars(tR,'{centers:cellRatio}=0.')
        tR = OversetInfo.oversetCellRatio__(tR,topTreeD)
    elif type == 'donorAspect':
        if loc == 'nodes': C._initVars(tR,'{donorAspect}=0.')
        else: C._initVars(tR,'{centers:donorAspect}=0.')
        tR = OversetInfo.oversetDonorAspect__(tR,topTreeD)

    else: raise ValueError("getOversetInfo: type is invalid.")
    return tR
