# -OversetInfo -
# Module for internal functions used for overset info
import Connector
import numpy
__version__ = Connector.__version__

try:
    import Converter.Internal as Internal
    import Converter.PyTree as C
    import Converter
except:
    raise ImportError("Connector.OversetInfo requires Converter module.")

#===============================================================================
# Information on overset quality : ratio of volumes between donors and receivers
# Compliant with setInterpolations
#===============================================================================
def chimeraCellRatio__(t):
    try: import Generator.PyTree as G
    except: raise ImportError("chimeraInfo: requires Generator module.")
    t = G.getVolumeMap(t)    
    zones = Internal.getZones(t)

    # First pass : direct storage
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        parentr,dr = Internal.getParentOfNode(t,zr)
        volR = C.getField('centers:vol',zr)[0][1]
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok                  
                ListVolD = Internal.getNodesFromName2(s,"InterpolantsVol")
                if ListVolD != []:
                    ListVolD = ListVolD[0][1]
                    ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]                
                    nindI = ListRcv.size                             
                    if nindI > 0:
                        field = Converter.array('cellRatio',nindI,1,1)
                        ListIndex = []                                 
                        for noind in xrange(nindI):
                            index = ListRcv[noind]
                            volr = volR[0,index]
                            vold = ListVolD[noind]
                            cr = max(volr/vold,vold/volr)
                            field[1][0,noind] = cr
                        #
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc='centers')
        #
        parentr[2][dr] = zr

    # 2nd pass: inverse storage
    zones = Internal.getZones(t)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(t,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []: 
                    zr = zr[0]                    
                    volR = C.getField('centers:vol',zr)[0][1]
                    parentr,dr = Internal.getParentOfNode(t,zr)
                    ListVolD = Internal.getNodesFromName2(s,"InterpolantsVol")
                    if ListVolD != []:
                        ListVolD = ListVolD[0][1]
                        ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                        nindI = ListRcv.size         
                        if nindI > 0:
                            field = Converter.array('cellRatio',nindI,1,1)
                            for noind in xrange(nindI):
                                index = ListRcv[noind]
                                volr = volR[0,index]
                                vold = ListVolD[noind]
                                cr = max(volr/vold,vold/volr)
                                field[1][0,noind] = cr

                            zr = C.setPartialFields(zr, [field], [ListRcv],loc='centers')                        
                            parentr[2][dr] = zr

    t = C.rmVars(t,'centers:vol') # faut il la detruire ou non ? pas de test leger pour savoir si c etait ds l arbre avant
    return t

#=============================================================================
# Information on overset quality : aspect ratio of donor cells
# Compliant with setInterpolations
#==============================================================================
def chimeraDonorAspect__(t):
    try: import Generator.PyTree as G
    except: raise ImportError("chimeraInfo: requires Generator module.")
    t = G.getEdgeRatio(t)  
    zones = Internal.getZones(t)
    
    # First pass : direct storage
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        parentr,dr = Internal.getParentOfNode(t,zr)
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                zdnrname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(t,zdnrname)
                zd = Internal.getNodesFromType1(zdonors,'Zone_t')
                if zd == []: raise ValueError("chimeraInfo: donor zone %s not found."%zdnrname)
                else: zd = zd[0]                
                #
                ListDnr = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                #
                ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]                
                nindI = ListRcv.size
                #
                if nindI > 0:
                    edgeRatio = C.getField('centers:EdgeRatio',zd)[0][1]
                    field = Converter.array('donorAspect',nindI,1,1)
                    for noind in xrange(nindI):     
                        indD = ListDnr[noind]          
                        field[1][0,noind] = edgeRatio[0,indD]
                    zr = C.setPartialFields(zr, [field], [ListRcv],loc='centers')
        parentr[2][dr] = zr
    
    # Inverse storage
    zones = Internal.getZones(t)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        dimsD = Internal.getZoneDim(zd)
        nidc = max(dimsD[1]-1,1)
        njdc = max(dimsD[2]-1,1)
        nidcnjdc = nidc*njdc
        edgeRatio = C.getField('centers:EdgeRatio',zd)[0][1]
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(t,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if (zr != []): 
                    zr = zr[0]                                       
                    parentr,dr = Internal.getParentOfNode(t,zr)
                    #
                    ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                    nindI = ListRcv.size         
                    #
                    if nindI > 0:
                        ListDnr = Internal.getNodesFromName1(s,'PointList')[0][1]
                        #
                        field = Converter.array('donorAspect',nindI,1,1)
                        for noind in xrange(nindI):
                            indD = ListDnr[noind]       
                            field[1][0,noind] = edgeRatio[0,indD]              
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc='centers')                        
                        parentr[2][dr] = zr

    t = C.rmVars(t,'centers:EdgeRatio')   
    return t

# =============================================================================
# extract info for interpolated/extrapolated/orphan cells
# Compliant with setInterpolations
# =============================================================================
def chimeraNatureOfCells__(t,nature):
    zones = Internal.getZones(t)
    # First pass : direct storage
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        parentr,dr = Internal.getParentOfNode(t,zr)
            
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                if nature == 'interpolated':
                    ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                    if ListRcv.size > 0:
                        field = Converter.array('interpolated',ListRcv.size,1,1)
                        field = Converter.initVars(field,'interpolated',1.)
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc='centers')

                elif nature == 'extrapolated':
                    ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                    if ListExtrap != []:
                        ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                        ListExtrap = ListExtrap[0][1]                       
                        DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1] 
                        Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                        # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                        field = Connector.connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                        zr = C.setPartialFields(zr, [field], [ListExtrap],loc='centers')       

                elif nature == 'orphan':
                    orphans = Internal.getNodesFromName(zr, 'OrphanPointList')
                    if orphans != []:
                        ListOrphan = orphans[0][1]
                        field = Converter.array('orphan',ListOrphan.size,1,1)
                        field = Converter.initVars(field,'orphan',1.)
                        zr = C.setPartialFields(zr, [field], [ListOrphan],loc='centers')          
        parentr[2][dr] = zr

    # 2nd pass: storage in donor zones
    zones = Internal.getZones(t)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(t,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if (zr != []): 
                    zr = zr[0]                    
                    parentr,dr = Internal.getParentOfNode(t,zr)
                    if nature == 'interpolated':
                        ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                        if ListRcv.size > 0:
                            field = Converter.array('interpolated',ListRcv.size,1,1)
                            field = Converter.initVars(field,'interpolated',1.)
                            zr = C.setPartialFields(zr, [field], [ListRcv],loc='centers')                        
                    elif nature == 'extrapolated':
                         ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                         if ListExtrap != []:
                             ListExtrap = ListExtrap[0][1]
                             DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1] 
                             Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                             ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                             # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                             field = Connector.connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                             zr = C.setPartialFields(zr, [field], [ListExtrap],loc='centers')       
                    elif nature == 'orphan':
                        orphans = Internal.getNodesFromName(zd, 'OrphanPointList')
                        if orphans != []:
                            ListOrphan = orphans[0][1]
                            field = Converter.array('orphan',ListOrphan.size,1,1)
                            field = Converter.initVars(field,'orphan',1.)
                            zr = C.setPartialFields(zr, [field], [ListOrphan],loc='centers')           
                    parentr[2][dr] = zr
    return t
# =============================================================================
# extract info for interpolated/extrapolated/orphan cells
# Compliant with setInterpData
# =============================================================================
def oversetNatureOfCells__(aR,topTreeD,nature):
    tR = Internal.copyRef(aR)
    zonesR = Internal.getZones(tR)
    # First pass : direct storage
    for zr in zonesR:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        parentr,dr = Internal.getParentOfNode(tR,zr)
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'
                if location == 'CellCenter': locr = 'centers'
                
                if nature == 'interpolated':
                    ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]                
                    if ListRcv.size > 0:
                        field = Converter.array('interpolated',ListRcv.size,1,1)
                        field = Converter.initVars(field,'interpolated',1.)
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)

                elif nature == 'extrapolated':
                    ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                    if ListExtrap != []:
                        ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                        ListExtrap = ListExtrap[0][1]
                        DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1] 
                        Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                        # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                        field = Connector.connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                        zr = C.setPartialFields(zr, [field], [ListExtrap],loc=locr)       
                elif nature == 'orphan':
                    orphans = Internal.getNodesFromName(zr, 'OrphanPointList')
                    if orphans != []:
                        ListOrphan = orphans[0][1]
                        field = Converter.array('orphan',ListOrphan.size,1,1)
                        field = Converter.initVars(field,'orphan',1.)
                        zr = C.setPartialFields(zr, [field], [ListOrphan],loc=locr)          
        parentr[2][dr] = zr

    # 2nd pass: storage in donor zones
    zones = Internal.getZones(topTreeD)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'
                if location == 'CellCenter': locr = 'centers'

                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(tR,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if (zr != []): 
                    zr = zr[0]                    
                    parentr,dr = Internal.getParentOfNode(tR,zr)
                    if nature == 'interpolated':
                        ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                        if ListRcv.size > 0:
                            field = Converter.array('interpolated',ListRcv.size,1,1)
                            field = Converter.initVars(field,'interpolated',1.)
                            zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)                        
                    elif nature == 'extrapolated':
                        ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                        if ListExtrap != []:
                            ListExtrap = ListExtrap[0][1]
                            DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1] 
                            Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                            ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                            # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                            field = Connector.connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                            zr = C.setPartialFields(zr, [field], [ListExtrap],loc=locr)       
                    elif nature == 'orphan':
                        orphans = Internal.getNodesFromName(zd, 'OrphanPointList')
                        if orphans != []:
                            ListOrphan = orphans[0][1]
                            field = Converter.array('orphan',ListOrphan.size,1,1)
                            field = Converter.initVars(field,'orphan',1.)
                            zr = C.setPartialFields(zr, [field], [ListOrphan],loc=locr)           
                    parentr[2][dr] = zr
    return tR

#===============================================================================
# Information on overset quality : ratio of volumes between donors and receivers
# Compliant with setInterpolations
#===============================================================================
def oversetCellRatio__(aR, topTreeD):
    try: import Generator.PyTree as G
    except: raise ImportError("oversetInfo: requires Generator module.")

    tR = Internal.copyRef(aR)
    tR = G.getVolumeMap(tR)    
    tR = C.center2Node(tR,'centers:vol')

    tD = Internal.copyRef(topTreeD)
    tD = G.getVolumeMap(tD)    
    tD = C.center2Node(tD,'centers:vol')
    tD = C.rmVars(tD,['centers:vol'])
    
    # First pass : direct storage
    zones = Internal.getZones(tR)
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        parentr,dr = Internal.getParentOfNode(tR,zr)
        volNR = C.getField('vol',zr)[0][1]
        volCR = C.getField('centers:vol',zr)[0][1]
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok      
                # interpolated node or cell ?
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'; volRcv = []
                if location == 'CellCenter': locr = 'centers'; volRcv = volCR
                else: volRcv = volNR
                ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]       
                # donor zone
                zdnrname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(tD,zdnrname)
                zd = Internal.getNodesFromType1(zdonors,'Zone_t')
                if zd == []: raise ValueError("oversetInfo: donor zone %s not found."%zdnrname)
                else: zd = zd[0]                
                volDnr = C.getField('vol',zd)[0][1]
                ListDnr = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
  
                nindI = ListDnr.size
                if nindI > 0:
                    field = Converter.array('cellRatio',nindI,1,1)
                    for noind in xrange(nindI):
                        voldl = volDnr[0,ListDnr[noind]]
                        volrl = volRcv[0,ListRcv[noind]]
                        cr = max(voldl/volrl,volrl/voldl)
                        field[1][0,noind] = cr
                    #
                    zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)   
        #
        parentr[2][dr] = zr

    # 2nd pass : inverse storage 
    zones = Internal.getZones(tD)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        volDnr = C.getField('vol',zd)[0][1]
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(tR,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if (zr != []): 
                    zr = zr[0]                   
                    parentr,dr = Internal.getParentOfNode(tR,zr) 
                    # interpolated node or cell ?
                    location = Internal.getNodesFromName1(s, 'GridLocation')
                    if location != []: location = Internal.getValue(location[0])
                    locr = 'nodes'; volRcv = []
                    if location == 'CellCenter': locr = 'centers'; volRcv = C.getField('centers:vol',zr)[0][1]
                    else: volRcv = C.getField('vol',zr)[0][1]
                    ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]       
                    ListDnr = Internal.getNodesFromName1(s,'PointList')[0][1]
                    #
                    nindI = ListDnr.size
                    if nindI > 0:
                        field = Converter.array('cellRatio',nindI,1,1)
                        for noind in xrange(nindI):
                            voldl = volDnr[0,ListDnr[noind]]
                            volrl = volRcv[0,ListRcv[noind]]
                            cr = max(voldl/volrl,volrl/voldl)
                            field[1][0,noind] = cr
                        #
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)  
                        parentr[2][dr] = zr
    #
    tR = C.rmVars(tR,'centers:vol') # faut il la detruire ou non ? pas de test leger pour savoir si c etait ds l arbre avant
    tR = C.rmVars(tR,'vol')
    tD = C.rmVars(tD,'vol')
    return tR
#===============================================================================
# Information on overset quality : ratio of volumes between donors and receivers
# Compliant with setInterpolations
#===============================================================================
def oversetDonorAspect__(aR, topTreeD):
    try: import Generator.PyTree as G
    except: raise ImportError("oversetInfo: requires Generator module.")
    # 
    tR = Internal.copyRef(aR)
    #
    tD = Internal.copyRef(topTreeD)
    tD = G.getEdgeRatio(tD)    
    tD = C.center2Node(tD,'centers:EdgeRatio')
    tD = C.rmVars(tD,['centers:EdgeRatio'])
    #
    # First pass : direct storage
    zones = Internal.getZones(tR)
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        parentr,dr = Internal.getParentOfNode(tR,zr)
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok      
                # interpolated node or cell ?
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes';
                if location == 'CellCenter': locr = 'centers'; 
              
                # donor zone
                zdnrname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(tD,zdnrname)
                zd = Internal.getNodesFromType1(zdonors,'Zone_t')
                if zd == []: raise ValueError("oversetInfo: donor zone %s not found."%zdnrname)      
                else: zd = zd[0]                
                ER = C.getField('EdgeRatio', zd)[0][1]
                ListDnr = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]         
                nindI = ListDnr.size
                if nindI > 0:
                    field = Converter.array('donorAspect',nindI,1,1)
                    for noind in xrange(nindI):             
                        field[1][0,noind] = ER[0,ListDnr[noind]]
                    #
                    zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)   
        #
        parentr[2][dr] = zr

    # 2nd pass : inverse storage 
    zones = Internal.getZones(tD)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0]
            if sname.split('_')[0] == 'ID': 
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        ER = C.getField('EdgeRatio',zd)[0][1]
        for s in subRegions:                
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(tR,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []: 
                    zr = zr[0]                   
                    parentr,dr = Internal.getParentOfNode(tR,zr) 
                    # interpolated node or cell ?
                    location = Internal.getNodesFromName1(s, 'GridLocation')
                    if location != []: location = Internal.getValue(location[0])
                    locr = 'nodes';
                    if location == 'CellCenter': locr = 'centers'
                    ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]       
                    ListDnr = Internal.getNodesFromName1(s,'PointList')[0][1]
                    #
                    nindI = ListDnr.size
                    if nindI > 0:
                        field = Converter.array('donorAspect',nindI,1,1)
                        for noind in xrange(nindI):
                            field[1][0,noind] = ER[0,ListDnr[noind]]
                        #
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)  
                        parentr[2][dr] = zr
    #
    tD = C.rmVars(tD,'EdgeRatio')
    return tR
