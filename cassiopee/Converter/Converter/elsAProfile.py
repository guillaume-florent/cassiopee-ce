import PyTree as C
import Internal
import Connector.PyTree as X
import Converter
import numpy
import math
from fnmatch import fnmatch

#############################################################################
# Gestion des recouvrements dans Cassiopee
# ----------------------------------------
# Les recouvrements sont associes a des connectivites de type 'Overset'.
# Dans l'arbre CGNS, les donnees sont donc stockees dans un noeud
# ZoneGridConnectivity_t avec une connectivite de type
# GridConnectivityType_t == 'Overset'.
# Les zones donneuses ne sont pas prefixees par le nom de leur base.
#
#
# Gestion des recouvrements dans elsAxdt
# ----------------------------------------
# Les recouvrements sont associes a des conditions aux limites de type
# 'BCOverlap'.
# Dans l'arbre CGNS, les donnees sont donc stockees dans un noeud ZoneBC_t avec
# un noeud BC_t de type FamilySpecified => Ces conditions aux limites sont donc
# obligatoirement definies via une famille associee.
# La famille associee possede un noeud FamilyBC_t == 'BCOverlap'.
# Les zones donneuses sont contenues dans un noeud NeighbourList, lui-meme
# contenu dans un noeud de type UserDefinedData appele .Solver#Overlap.
# Les zones donneuses sont prefixees par le nom de leur base
#
#
# Conversion proposee
# -------------------
# Creation d'une famille (ou deux s'il y a des conditions doubly_defined) par
# zone.
# Dans un premier temps, les ZoneBC_t de type BCOverlap sont construites.
# Dans un second temps, on remplit les NeighbourList en tenant compte des
# intersections des "bounding box" des domaines.
#############################################################################
#
#                      QUELQUES LIMITATIONS
# ---------------------------------------------------------------------------
# La conversion ne gere que le passage de donnees creees par Cassiopee.
# Par exemple, la transition ou les conditions aux limites demandant des
# parametres "utilisateurs" ne sont pas gerees
#
##############################################################################

# Traduction : keyword elsA -> keyword CGNS
keyselsA2CGNS = {\
'config'         :'EquationDimension'             , \
'fluid'          :'GasModel'                      , \
''               :'GasModelType'                  , \
'pg'             :'Ideal'                         , \
#'pg'            :'CalloricallyPerfect'            , \
'gamma'          :'SpecificHeatRatio'             , \
'cv'             :'SpecificHeatVolume'            , \
'_'              :'IdealGasConstant'              , \
'__'             :'ThermalConductivityModel'      , \
'___'            :'TurbulenceClosure'             , \
'prandtl'        :'ConstantPrandtl'               , \
'prandtltb'      :'PrandtlTurbulent'              , \
'phymod'         :"GoverningEquations"            , \
'turbmod'        :'TurbulenceModel'               , \
'euler'          :'Euler'                         , \
'nstur'          :'NSTurbulent'                   , \
'nslam'          :'NSLaminar'                     , \
'spalart'        :'OneEquation_SpalartAllmaras'   , \
'kepsjl'         :'TwoEquation_JonesLaunder'      , \
'komega_menter'  :'TwoEquation_MenterSST'         , \
'komega_wilcox'  :'TwoEquation_Wilcox'            , \
'komega_kok'     :'UserDefined'                   , \
'smith'          :'UserDefined'                   , \
'visclaw'        :'ViscosityModel'                , \
'sutherland'     :'Sutherland'                    , \
'suth_const'     :'SutherlandLawConstant'         , \
'suth_muref'     :'ViscosityMolecularReference'   , \
'suth_tref'      :'TemperatureReference'          , \
'walladia'       :'BCWall'                        , \
'wallslip'       :'BCWallInviscid'                , \
'cell'           :'CellCenter'                    , \
'node'           :'Vertex'                        , \
'x'              :'CoordinateX'                   , \
'y'              :'CoordinateY'                   , \
'z'              :'CoordinateZ'                   , \
'X'              :'CoordinateX'                   , \
'Y'              :'CoordinateY'                   , \
'Z'              :'CoordinateZ'                   , \
'ro'             :'Density'                       , \
'rou'            :'MomentumX'                     , \
'rov'            :'MomentumY'                     , \
'row'            :'MomentumZ'                     , \
'rovx'           :'MomentumX'                     , \
'rovy'           :'MomentumY'                     , \
'rovz'           :'MomentumZ'                     , \
'roe'            :'EnergyStagnationDensity'       , \
'roE'            :'EnergyStagnationDensity'       , \
'rok'            :'TurbulentEnergyKineticDensity' , \
'roeps'          :'TurbulentDissipationDensity'   , \
'mach'           :'Mach'                          , \
'psta'           :'Pressure'                      , \
'tsta'           :'Temperature'                   , \
'viscrapp'       :'Viscosity_EddyMolecularRatio'  , \
'walldistance'   :'TurbulentDistance'             , \
'wallglobalindex':'TurbulentDistanceIndex'        , \
}

# Traduction : keyword CGNS -> keyword elsA
# keysCGNS2elsA = dict((v, k) for k, v in keyselsA2CGNS.iteritems())
keysCGNS2elsA={
'EquationDimension'             :'config'       , \
'GasModel'                      :'fluid'        , \
'GasModelType'                  :''             , \
'Ideal'                         :'pg'           , \
#'CalloricallyPerfect'           :'pg'           , \
'SpecificHeatRatio'             :'gamma'        , \
'SpecificHeatVolume'            :'cv'           , \
'IdealGasConstant'              :'_'            , \
'ThermalConductivityModel'      :'__'           , \
'TurbulenceClosure'             :'___'          , \
'ConstantPrandtl'               :'prandtl'      , \
'PrandtlTurbulent'              :'prandtltb'    , \
"GoverningEquations"            :'phymod'       , \
'TurbulenceModel'               :'turbmod'      , \
'Euler'                         :'euler'        , \
'NSTurbulent'                   :'nstur'        , \
'NSLaminar'                     :'nslam'        , \
'OneEquation_SpalartAllmaras'   :'spalart'      , \
'TwoEquation_JonesLaunder'      :'kepsjl'       , \
'TwoEquation_MenterSST'         :'komega_menter', \
'TwoEquation_Wilcox'            :'komega_wilcox', \
'UserDefined'                   :'komega_kok'   , \
'UserDefined'                   :'smith'        , \
'ViscosityModel'                :'visclaw'      , \
'Sutherland'                    :'sutherland'   ,\
'SutherlandLawConstant'         :'suth_const'   , \
'ViscosityMolecularReference'   :'suth_muref'   , \
'TemperatureReference'          :'suth_tref'    , \
'BCWall'                        :'walladia'     , \
'BCWallInviscid'                :'wallslip'     , \
'CellCenter'                    :'cell'         , \
'Vertex'                        :'node'         , \
'CoordinateX'                   :'x'            , \
'CoordinateY'                   :'y'            , \
'CoordinateZ'                   :'z'            , \
'Density'                       :'ro'           , \
'MomentumX'                     :'rou'          , \
'MomentumY'                     :'rov'          , \
'MomentumZ'                     :'row'          , \
'EnergyStagnationDensity'       :'roe'          , \
"TurbulentEnergyKineticDensity" :'rok'          , \
"TurbulentDissipationDensity"   :'roeps'        , \
'Mach'                          :'mach'         , \
'Pressure'                      :'psta'         , \
'Temperature'                   :'tsta'         , \
'Viscosity_EddyMolecularRatio'  :'viscrapp'     , \
}

# -----------------------------------------------------------------------------
def getCGNSkeys(key, verbose=True):
    """Return the CGNS key, if it exists, equivalent to the elsA key.
    Usage: getCGNSkeys(elsAKey)"""
    if key in keyselsA2CGNS.keys(): return keyselsA2CGNS[key]
    elif key in keyselsA2CGNS.values(): return key
    else:
        if verbose: print 'Warning: getCGNSkeys: the given key '" + key + "' could not be translated in a CGNS key.'
        return key

# -----------------------------------------------------------------------------
# Cree un noeud output contenant les informations pour les sorties elsA
def createOutput(Dict, name):
    """Create an output node containing informations for elsA extraction.
       - `Dict` - a python dictionary containing all output informations. For example::
            >>> Dict={}
            >>> Dict["var"]="convflux_rou convflux_rov convflux_row"
            >>> Dict["loc"]=1
            >>> Dict["fluxcoef"]=1.
            >>> Dict["period"]=10

       - `name` - a string to complete the node name.
    """
    node = Internal.createNode(".Solver#Output"+name, 'UserDefinedData_t')
    try:
        for each in Dict.keys():
            Internal.createChild(node, each, "DataArray_t", Dict[each])
        return node
    except: raise TypeError("createOutput: the first argument is not a python dictionary.")

# -----------------------------------------------------------------------------
def _addOutput(node, Dict, name):
    """Add an output node to extract required value from a node.
       - `node` - node pyTree concerned with extraction.
       - `Dict` - a python dictionary containing all output informations. For example::
            >>> Dict={}
            >>> Dict["var"]="convflux_rou convflux_rov convflux_row"
            >>> Dict["loc"]=1
            >>> Dict["fluxcoef"]=1.
            >>> Dict["period"]=10
       - `name` - a string to complete the node name.
    Returns
       - `outputNode` - the newly created output node
    """
    outputNode = createOutput(Dict, name)
    Internal._addChild(node, outputNode)
    return None

def addOutput(node, Dict, name):
  """Add an output node to extract required value from a node.
  """
  nodep = Internal.copyRef(node)
  _addOutput(nodep, Dict, name)
  return nodep

# -----------------------------------------------------------------------------
def addOutputForces(node, name="", var=None, loc=4, writingmode=1,
                    period=None, pinf=None, fluxcoef=None,
                    torquecoef=None, xyztorque=None, frame=None,
                    governingEquations="NSTurbulent", xtorque=None,
                    ytorque=None, ztorque=None):
    """Add an output node for forces.
       - `node` - a pyTree node concerned with extraction.
       - `name` - an optional string to complete the output node name.
       - `var` - an optional list of string of the variable to extract.
       - `loc` - an optional integer for the location of the extraction.
       - `writingmode` - an optional integer of the writing mode.
       - `period` - an optional integer of file writing frequency.
       - `pinf` - an optional real of the infinite stagnation pressure.
       - `fluxcoef` - an optional real to correct the extract value.
       - `torquecoef` - an optional real to correct the extract value.
       - `xyztorque` - an optional list of real of the torque origin.
       - `frame` - an optional string of the writing frame. ['relative','absolute']. Default: elsA default value.
       - `governingEquations` - an optional string of the governing equations of the model.
    """
    nodep = Internal.copyRef(node)
    _addOutputForces(nodep, name, var, loc, writingmode,
                     period, pinf, fluxcoef,
                     torquecoef, xyztorque, frame,
                     governingEquations, xtorque, ytorque, ztorque)
    return nodep

def _addOutputForces(node, name="", var=None, loc=4, writingmode=1,
                     period=None, pinf=None, fluxcoef=None,
                     torquecoef=None, xyztorque=None, frame=None,
                     governingEquations="NSTurbulent", xtorque=None,
                     ytorque=None, ztorque=None):
    Dict = {}
    Dict['var'] = "flux_rou flux_rov flux_row convflux_rou convflux_rov convflux_row diffflux_rou diffflux_rov diffflux_row torque_rou torque_rov torque_row convtorque_rou convtorque_rov convtorque_row difftorque_rou difftorque_rov difftorque_row"
    Dict["writingmode"] = writingmode
    Dict["loc"] = loc
    if var is not None: Dict["var"] = var
    elif governingEquations == "Euler": Dict["var"] = "convflux_rou convflux_rov convflux_row convtorque_rou convtorque_rov convtorque_row"
    if period is not None: Dict["period"] = period
    if pinf is not None: Dict["pinf"] = pinf
    if fluxcoef is not None: Dict["fluxcoeff"] = fluxcoef
    if torquecoef is not None: Dict["torquecoeff"] = torquecoef
    if xtorque is not None:
        Dict["xtorque"] = xtorque
        deprecation("'xtorque' is obsolete. Please use 'xyztorque' argument instead.", stacklevel=3)
    if ytorque is not None:
        Dict["ytorque"] = ytorque
        deprecation("'ytorque' is obsolete. Please use 'xyztorque' argument instead.", stacklevel=3)
    if ztorque is not None:
        Dict["ztorque"] = ztorque
        deprecation("'ztorque' is obsolete. Please use 'xyztorque' argument instead.", stacklevel=3)
    if xyztorque is not None:
        Dict["xtorque"] = xyztorque[0]
        Dict["ytorque"] = xyztorque[1]
        Dict["ztorque"] = xyztorque[2]
    if writingmode is not None: Dict["writingmode"] = writingmode
    if frame is not None:
        if frame not in ['relative','absolute']:
            raise AttributeError('Frame should be in %s')%(['relative','absolute'])
        Dict["writingframe"] = frame
    _addOutput(node, Dict, ":Forces"+name)
    return None

#==============================================================================
def addOutputFriction(node, name="", var=None, loc=4, writingmode=1,
                      period=None, fluxcoef=None,
                      torquecoef=None, writingframe=None):
    """Add an output node for frictions.
       - `node` - a pyTree node concerned with extraction.
       - `name` - an optional string to complete the output node name.
       - `var` - an optional list of string of the variable to extract.
                 Default='SkinFrictionX SkinFrictionY SkinFrictionZ SkinFrictionMagnitude WallCellSize'
       - `loc` - an optional integer for the location of the extraction.
       - `writingmode` - an optional integer of the writing mode.
       - `period` - an optional integer of file writing frequency.
       - `fluxcoef` - an optional real to correct the extract value.
       - `writingframe` - an optional integer of the writing frame.
    """
    nodep = Internal.copyRef(node)
    _addOutputFriction(node, name, var, loc, writingmode,
                       period, fluxcoef,
                       torquecoef, writingframe)
    return nodep

def _addOutputFriction(node, name="", var=None, loc=4, writingmode=1,
                       period=None, fluxcoef=None,
                       torquecoef=None, writingframe=None):
    Dict = {}
    Dict['var'] = "SkinFrictionX SkinFrictionY SkinFrictionZ SkinFrictionMagnitude WallCellSize"
    Dict["writingmode"] = writingmode
    Dict["loc"] = loc
    if var is not None: Dict["var"] = var
    if period is not None: Dict["period"] = period
    if fluxcoef is not None: Dict["fluxcoeff"] = fluxcoef
    if writingmode is not None: Dict["writingmode"] = writingmode
    if writingframe is not None: Dict["writingframe"] = writingframe
    _addOutput(node, Dict, ":Friction"+name)
    return None

#=============================================================================
def addGlobalConvergenceHistory(t, normValue=0):
   """Create a node for global convergence history storage for each base.
   Usage: addGlobalConvergenceHistory(t, normValue)
   IN: pyTree t
   IN: normValue - an optional integer, specifying the type of norm.
   """
   tp = Internal.copyRef(t)
   _addGlobalConvergenceHistory(tp, normValue)
   return tp

def _addGlobalConvergenceHistory(t, normValue=0):
   bases = Internal.getBases(t)
   for base in bases:
       child = Internal.createNode("NormDefinitions", "Descriptor_t", "ConvergenceHistory")
       Internal._addChild(base,Internal.createNode("GlobalConvergenceHistory","ConvergenceHistory_t", normValue, [child]))
   return None

#=============================================================================
def addReferenceState(t, conservative=None, temp=None, turbmod='spalart',
                      name='ReferenceState',
                      comments=None, k=None, eps=None):
   """Add a reference state for each base.
   Usage: addreferenceState(t, conservative, temp, turbmod, name, comments, k, eps)
   IN: t: pyTree containing CGNS bases
   IN: conservative: list of all conservative variables [ro,rovx,rovy,rovz,roe,rok,roeps] in this order. Only the first five are required.
   IN: temp: optional real, temperature.
   IN: turbmod: - turbulence model in order to set the correct names of the turbulent variables
         - `komega` - k-omega type models
         - `spalart` - Spalart-Allmaras model
         - `keps` - k-epsilon type models
         - `rsm` - DRSM type models
   IN: name: string. The name of the ReferenceState node (default: 'ReferenceState' )
   IN: comments: optional string, comments to add to the flow state description.
   """
   tp = Internal.copyRef(t)
   _addReferenceState(tp, conservative, temp, turbmod,
                       name, comments, k, eps)
   return tp

def _addReferenceState(t, conservative=None, temp=None, turbmod='spalart',
                       name='ReferenceState',
                       comments=None, k=None, eps=None):
   bases = Internal.getBases(t)
   for b in bases:
       NoeudReferenceState = Internal.createNode(name, 'ReferenceState_t')
       if conservative is not None and len(conservative) >= 5:
           Internal.createChild(NoeudReferenceState,'Density','DataArray_t', conservative[0])
           Internal.createChild(NoeudReferenceState,'MomentumX','DataArray_t', conservative[1])
           Internal.createChild(NoeudReferenceState,'MomentumY','DataArray_t', conservative[2])
           Internal.createChild(NoeudReferenceState,'MomentumZ','DataArray_t', conservative[3])
           Internal.createChild(NoeudReferenceState,'EnergyStagnationDensity','DataArray_t', conservative[4])
           if len(conservative) == 5: pass
           elif len(conservative) == 6:
               if turbmod == 'spalart':
                   Internal.createChild(NoeudReferenceState,'TurbulentSANuTildeDensity','DataArray_t', conservative[5])
               else:
                   raise Exception('Number of conservative variables: conservative='+str(conservative)+' is not compatible with choice of turbulence model: turbmod=\''+turbmod+'\'.')
           elif len(conservative) == 7:
               if turbmod[0:6] =='komega':
                  Internal.createChild(NoeudReferenceState,'TurbulentEnergyKineticDensity','DataArray_t', conservative[5])
                  Internal.createChild(NoeudReferenceState,'TurbulentDissipationRateDensity','DataArray_t', conservative[6])

               elif turbmod[0:4] =='keps' or turbmod == 'chien' or turbmod == 'asm':
                   Internal.createChild(NoeudReferenceState,'TurbulentEnergyKineticDensity','DataArray_t', conservative[5])
                   Internal.createChild(NoeudReferenceState,'TurbulentDissipationDensity','DataArray_t', conservative[6])

               elif turbmod == 'smith':
                   Internal.createChild(NoeudReferenceState,'TurbulentEnergyKineticDensity','DataArray_t', conservative[5])
                   Internal.createChild(NoeudReferenceState,'TurbulentLengthScaleDensity','DataArray_t', conservative[6])

               elif turbmod == 'kkl' or turbmod[0:5] == 'earsm':
                   Internal.createChild(NoeudReferenceState,'TurbulentEnergyKineticDensity','DataArray_t', conservative[5])
                   Internal.createChild(NoeudReferenceState,'TurbulentKineticPLSDensity','DataArray_t', conservative[6])

               else:
                   raise Exception('Number of conservative variables: conservative='+str(conservative)+' is not compatible with choice of turbulence model: turbmod=\''+turbmod+'\'.')
           elif len(conservative) == 12:
               if turbmod == 'rsm':
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressXX','DataArray_t', conservative[5])
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressXY','DataArray_t', conservative[6])
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressXZ','DataArray_t', conservative[7])
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressYY','DataArray_t', conservative[8])
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressYZ','DataArray_t', conservative[9])
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressZZ','DataArray_t', conservative[10])
                   Internal.createChild(NoeudReferenceState,'ReynoldsStressDissipationScale','DataArray_t', conservative[11])
               else:
                   raise Exception('Number of conservative variables: conservative='+str(conservative)+' is not compatible with choice of turbulence model: turbmod=\''+turbmod+'\'.')
           else:
               raise Exception('Number of conservative variables: conservative='+str(conservative)+' is not supported.')
       else:
           raise Exception('Number of conservative variables: conservative='+str(conservative)+' is not supported.')
       if temp is not None: Internal.createChild(NoeudReferenceState, 'Temperature','DataArray_t', temp)
       Internal.createChild(NoeudReferenceState,'ReferenceStateDescription','Descriptor_t', str(comments))
       Internal._addChild(b, NoeudReferenceState)
   return None

#===============================================================================
def createFlowSolution(name='', loc='CellCenter', variables=None,
                       governingEquations=None, writingMode=None,
                       writingFrame=None, period=None, output=None):
    """Create a node to extract flow solution.
       - `name` - an optional string of the suffix node name to add to 'FlowSolution'.
       - `loc` - an optional string of the location of extraction (``CellCenter``, ``Vertex`` or ``cellfict``).
       - `variables` - an optional list of variables.
       - `governingEquations` - an optional string to precise what kind of variables (``Euler``, ``NSTurbulent``, ...). If given, ``Conservative`` and ``Commons`` will be automatically added to `Variables`.
       - `writingMode` - an optional integer.
       - `writingFrame` - an optional string in ['absolute','relative']
       - `period` - an optional integer.
       - `output` - an optional dictionary of node to add to '.Solver#Output' if necessary.

    :Keywords:
       - ``Conservative`` for all conservative variables.
       - ``Turbulent`` for turbulent variables.
       - ``Commons`` for *Pressure*, *Mach* and *Temperature*, plus ``CommonsNS`` for NSTurbulent model if given.
       - ``CommonsNS`` for *Viscosity_EddyMolecularRatio* ,
       - ``xyz`` for coordinates.
       - ``WallDistance`` for *TurbulentDistance* and *TurbulentDistanceIndex*
    """
    def addVariable(var, List, ind=None):
        if ind is None: ind = len(List)
        if var not in List: List.insert(ind, getCGNSkeys(var, verbose=False))
    # ----------------------
    def addVariables(Vars, List):
        Conservative = ["Density", "MomentumX", "MomentumY", "MomentumZ", "EnergyStagnationDensity"]
        Turbulent = ["TurbulentEnergyKineticDensity", "TurbulentDissipationDensity"]
        Commons = ["Pressure", "Mach", "Temperature"]
        CommonsNS = ["Viscosity_EddyMolecularRatio"]
        WallDistance = ["TurbulentDistance","TurbulentDistanceIndex"]
        xyz = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
        ind = List.index(Vars)
        List.pop(ind)
        VARS = eval(Vars)
        for var in VARS:
            addVariable(var, List, ind + VARS.index(var))

    if variables is None: variables = []
    elif isinstance(variables, str): variables=variables.split()
    if governingEquations is not None:
        addVariable("Conservative", variables, 0)
        addVariable("Commons", variables, 1)
    else: governingEquations = 'GoverningEquations'
    if "xyz" in variables:
        addVariables("xyz", Variables)
    if "Conservatives" in variables:
        deprecation("'Conservatives' keywords is obsolete. Please use 'Conservative' instead => Forced to 'Conservative'")
        ind = variables.index('Conservatives')
        variables.pop(ind)
        variables.insert(ind, 'Conservative')
    if "Turbulents" in variables:
        deprecation("'Turbulents' keywords is obsolete. Please use 'Turbulent' instead => Forced to 'Turbulent'")
        ind = variables.index('Turbulents')
        variables.pop(ind)
        variables.insert(ind, 'Turbulent')
    if "Conservative" in variables:
        addVariables("Conservative", variables)
    if "Turbulent" in variables:
        addVariables("Turbulent", variables)
    if "Commons" in variables:
        if fnmatch(governingEquations, "NS*"): addVariable("CommonsNS", variables, variables.index("Commons"))
        addVariables("Commons", variables)
    if "CommonsNS" in variables:
        addVariables("CommonsNS", variables)
    if "WallDistance" in variables:
        addVariables("WallDistance", variables)

    if output is None: output = {}
    if loc == 'cellfict': output['loc'] = 2
    elif getCGNSkeys(loc) == 'CellCenter' or getCGNSkeys(loc) == 'Vertex': pass
    else: raise AttributeError, "'loc' attribute should be 'CellCenter','Vertex' or 'cellfict' => loc='%s'" % (loc)
    if writingMode is not None: output['writingmode'] = writingMode
    if period is not None: output['period'] = period
    if writingFrame is not None: output['writingframe'] = writingFrame
    FlowSolution = Internal.createNode("FlowSolution"+name,"FlowSolution_t")
    Internal.createChild(FlowSolution, "GridLocation", "GridLocation_t", loc)
    if output != {}: _addOutput(FlowSolution, output, '')
    for v in variables:
        Internal.createChild(FlowSolution, v, "UserDefinedData_t")
    return FlowSolution

#===============================================================================
def addFlowSolution(t, name='', loc='CellCenter', variables=None,
                    governingEquations=None, writingMode=None,
                    writingFrame='relative', period=None, output=None,
                    addBCExtract=False, protocol="end"):
    """Add a node to extract flow solution for each zone of the pyTree.
       - `name` - an optional string of the suffix node name to add to 'FlowSolution'.
       - `loc` - an optional string of the location of extraction (``CellCenter``, ``Vertex`` or ``cellfict``).
       - `variables` - an optional list of variables.
       - `governingEquations` - an optional string to precise what kind of variables (``Euler``, ``NSTurbulent``, ...). If given, ``Conservative`` and ``Commons`` will be automatically added to `Variables`.
       - `writingMode` - an optional integer.
       - `writingFrame` - an optional string in ['absolute','relative']
       - `period` - an optional integer.
       - `output` - an optional dictionary of node to add to '.Solver#Output' if necessary.
       - `addBCExtract` - an optional boolean to add extract from BC windows (pseudo CellFict).
       - `protocol` - an optional string of the value of 'Protocol' node ('iteration','end','after').

    :Keywords:
       - ``Conservative`` for all conservative variables.
       - ``Turbulent`` for turbulent variables.
       - ``Commons`` for *Pressure*, *Mach* and *Temperature*, plus ``CommonsNS`` for NSTurbulent model if given.
       - ``CommonsNS`` for *Viscosity_EddyMolecularRatio*.
       - ``xyz`` for coordinates.
    """
    tp = Internal.copyRef(t)
    _addFlowSolution(tp, name, loc, variables,
                     governingEquations, writingMode,
                     writingFrame, period, output,
                     addBCExtract, protocol)
    return tp

def _addFlowSolution(t, name='', loc='CellCenter', variables=None,
                     governingEquations=None, writingMode=None,
                     writingFrame='relative', period=None, output=None,
                     addBCExtract=False, protocol="end"):
    import string
    for zone in Internal.getZones(t):
        FlowSolution = createFlowSolution(name=name, loc=loc, variables=variables, governingEquations=governingEquations, writingMode=writingMode, writingFrame=writingFrame, period=period, output=output)
        if addBCExtract:
            child = Internal.createNode(".Solver#SolutionSubSetInBC", "UserDefinedData_t")
            for bc in Internal.getNodesFromType2(zone, "BC_t"):
                subChild = Internal.createNode(bc[0], "UserDefinedData_t")
                Internal.createChild(subChild, "PointRange", "DataArray_t", Internal.getNodesFromName1(bc, "PointRange")[0][1])
                Internal.createChild(subChild, "Protocol", "Descriptor_t", protocol)
                Internal.createChild(subChild, "Path", "DataArray_t", string.join(['SurfaceSolution/NeumannData/'+variable for variable in variables],'\n'))
                Internal._addChild(child, subChild)

            Internal._addChild(FlowSolution, child)
        Internal._addChild(zone, FlowSolution)
    return None

#==============================================================================
def addFlowSolutionEoR(t, name='', variables=None, governingEquations=None,
                       writingFrame='relative', addBCExtract=False,
                       protocol="end"):
    """Add a node to extract flow solution at End of Run for each zone of the pyTree.
       - `name` - an optional string of the suffix node name to add to 'FlowSolution#EndOfRun'.
       - `variables` - an optional list of variables.
       - `writingFrame` - an optional string in ['absolute','relative']
       - `governingEquations` - an optional string to precise what kind of variables (``Euler``, ``NSTurbulent``, ...). If given, ``Conservative`` and ``Commons`` will be automatically added to `variables`.
       - `addBCExtract` - an optional boolean to add extract from BC windows (pseudo CellFict).
       - `protocol` - an optional string of the value of 'Protocol' node ('iteration','end','after').

    :Keywords:
       - ``Conservative`` for all conservative variables.
       - ``Turbulent`` for turbulent variables.
       - ``Commons`` for *Pressure*, *Mach* and *Temperature*, plus ``CommonsNS`` for NSTurbulent model if given.
       - ``CommonsNS`` for *Viscosity_EddyMolecularRatio*.
       - ``xyz`` for coordinates.
    """
    tp = Internal.copyRef(t)
    _addFlowSolutionEoR(tp, name, variables, governingEquations,
                        writingFrame, addBCExtract, protocol)
    return tp

def _addFlowSolutionEoR(t, name='', variables=None, governingEquations=None,
                        writingFrame='relative', addBCExtract=False,
                        protocol="end"):
    _addFlowSolution(t, name='#EndOfRun'+name, loc='CellCenter', governingEquations=governingEquations, writingFrame=writingFrame, variables=variables, addBCExtract=addBCExtract, protocol=protocol)
    return None

#==============================================================================
def buildBCOverlap__(t):
    """ Convert the Overlap condition from GC (Grid Connectivity condition) to BC (Boundary Condition) for elsA solver.
    """
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)

    c = 0 # compteur pour le nommage des conditions BCOverlap
    for i in xrange(len(bases)):
        zones = Internal.getNodesFromType1(bases[i], 'Zone_t')
        for j in xrange(len(zones)):
            (parentBase, numZone) = Internal.getParentOfNode(tp,zones[j])
            # Creation d'une famille par zone"
            familyName='F_'+zones[j][0]
            hasOverset = 0
            zonegct = Internal.getNodesFromType(zones[j],'GridConnectivityType_t')
            for zgct in zonegct:
                valz = Internal.getValue(zgct)
                if valz == 'Overset': hasOverset = 1
            if hasOverset == 0: # pas de connectivite Overset dans la zone => famille "vide"
                # Creation de la famille
                bases[i] = C.addFamily2Base(bases[i], familyName, 'UserDefined')
                # Ajout pour la zone d'un noeud fils FamilyName, pour faire reference a la famille creee
                famNameArray = numpy.fromstring(familyName, 'c')
                familyNameList = Internal.getNodesFromType1(zones[j], 'FamilyName_t')
                if familyNameList == []:
                    zones[j][2].append(['FamilyName', famNameArray, [], 'FamilyName_t'])
                else:
                    addFamilyNameList = Internal.getNodesFromType1(zones[j],'AdditionalFamilyName_t')
                    if addFamilyNameList == []:
                        zones[j][2].append(['AddFamilyName', famNameArray, [], 'AdditionalFamilyName_t'])
            else:
                bases[i] = C.addFamily2Base(bases[i], familyName, 'BCOverlap')
                # Creation d'un noeud NeighbourList dans le noeud .Solver#Overlap de la famille creee
                F = Internal.getNodesFromName1(bases[i], familyName)[0]
                (parentBase, numFamily) = Internal.getParentOfNode(bases[i],F)
                Ovlp = Internal.getNodesFromName(F,'.Solver#Overlap')[0]
                (parentFamily, numOvlp) = Internal.getParentOfNode(F,Ovlp)
                F[2][numOvlp][2].append(['NeighbourList', None, [], 'DataArray_t'])
                bases[i][2][numFamily] = F
                # Creation des noeuds ZoneBC_t de type BCOverlap
                gc = Internal.getNodesFromType2(zones[j],'GridConnectivity_t')
                for k in xrange(len(gc)):
                    # search in parent GridConnectivity if a PointRange is present
                    # if no, connectivity is linked to blanking and not to a BCOverlap
                    prange = Internal.getNodesFromName1(gc[k],'PointRange')
                    if prange != []: # corresponds to a BCOverlap
                        gct = Internal.getNodesFromType3(gc[k],'GridConnectivityType_t')
                        for o in gct:
                            if o != []: val=o[1]
                            if isinstance(val, numpy.ndarray):
                                val = val.tostring()
                            if val == 'Overset':
                                # Recuperation des donnees necessaires a la creation d'un noeud ZoneBC_t
                                #   range
                                lPointRange = Internal.getNodesFromName1(gc[k], 'PointRange')
                                r = Internal.range2Window(lPointRange[0][1])
                                i1=int(r[0]); j1=int(r[2]); k1=int(r[4])
                                i2=int(r[1]); j2=int(r[3]); k2=int(r[5])
                                range = [i1,i2,j1,j2,k1,k2]
                                #   nom de la ZoneBC_t
                                overlapName='overlapBC'+str(c); c=c+1
                                #   doubly_defined
                                userDef = Internal.getNodesFromName(gc[k], 'UserDefinedData')
                                if userDef != []:
                                    if (len(userDef[0]) == 4):
                                        info = userDef[0][2][0]
                                        if (info[0] == 'doubly_defined'):
                                            # Creation d'une famille doubly_defined pour la zone"
                                            familyNameDD='FDD_'+zones[j][0]+'_'+gc[k][0]
                                            ListFDD = (Internal.getNodesFromName1(bases[i],familyNameDD))
                                            if ListFDD == []:
                                                bases[i] = C.addFamily2Base(bases[i], familyNameDD, 'BCOverlap')
                                            zones[j] = C.addBC2Zone(zones[j], overlapName, 'FamilySpecified:'+familyNameDD,range, rangeDonor='doubly_defined')
                                            zones[j] = C.tagWithFamily(zones[j],familyNameDD)
                                            FDD = Internal.getNodesFromName1(bases[i],familyNameDD)[0]
                                            (parentBaseDD, numFamilyDD) = Internal.getParentOfNode(bases[i],FDD)
                                            #    Creation d'un noeud doubly_defined
                                            OvlpDD  = Internal.getNodesFromName(FDD,'.Solver#Overlap')[0]
                                            (parentFamilyDD, numOvlpDD) = Internal.getParentOfNode(FDD,OvlpDD)
                                            FDD[2][numOvlpDD][2].append(['NeighbourList', None, [], 'DataArray_t'])
                                            dd = numpy.fromstring('active', 'c')
                                            FDD[2][numOvlpDD][2].append(['doubly_defined', dd, [], 'DataArray_t'])
                                            bases[i][2][numFamilyDD] = FDD
                                else:
                                    zones[j] = C.addBC2Zone(zones[j], overlapName, 'FamilySpecified:'+familyName,range)
                                    familyNameList = Internal.getNodesFromType1(zones[j],'FamilyName_t')
                                    if familyNameList == []: zones[j] = C.tagWithFamily(zones[j],familyName)

            bases[i][2][numZone] = zones[j]

    cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv != []:
        (parent, pos) = Internal.getParentOfNode(tp,cgnsv[0])
        bases = [tp[2][pos]]+bases
    tp[2] = bases

    return tp

#==============================================================================
def buildBCOverlap(t):
    """ Convert the Overlap condition from GC (Grid Connectivity condition) to BC (Boundary Condition) for elsA solver.
    """
    tp = buildBCOverlap__(t)
    return tp

#==============================================================================
def rmGCOverlap__(t):
    """ Remove the Overlap condition described as Grid Connectivity.
    """
    tp = Internal.copyRef(t)
    zgc = Internal.getNodesFromType3(tp,'ZoneGridConnectivity_t')
    for z in zgc:
      gc = Internal.getNodesFromType1(z,'GridConnectivity_t')
      for node in gc:
        gct = Internal.getNodeFromType1(node,'GridConnectivityType_t')
        if gct is not None: Internal._rmNode(z,node)
    return tp

#==============================================================================
def rmGCOverlap(t):
    """ Remove the Overlap condition described as Grid Connectivity.
    """
    tp = rmGCOverlap__(t)
    return tp

#==============================================================================
def addNeighbours__(t, sameBase=0):
    """ Fill the NeighbourList nodes with bounding-box domains intersection.
    """
    # On regarde si les zones donneuses existent deja dans Cassiopee
    #   - Si oui, on les recopie
    #   - Si non, on les construit via getBBIntersectingDomainsForBase
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)
    for i in xrange(len(bases)):
        fams = []
        doms = X.getCEBBIntersectingDomains(bases[i] , bases, sameBase)
        if doms != [[]]:
            for j in xrange(len(doms)):
                famsByZone=[]
                if doms[j] !=[]:
                    for donorZoneName in doms[j]:
                        z = Internal.getNodesFromName2(tp,donorZoneName)
                        if z != []:
                            (donorBase, num) = Internal.getParentOfNode(tp,z[0])
                            baseName=donorBase[0]
                            donorFamilyName='F_'+ donorZoneName
                            donorF = Internal.getNodesFromName2(tp,donorFamilyName)[0]
                            famsByZone.append(baseName+'/'+donorF[0])
                fams.append(famsByZone)
        zones = Internal.getNodesFromType1(bases[i],'Zone_t')
        for zi in xrange(len(zones)):
            familyName='F_'+ zones[zi][0]
            F = Internal.getNodesFromName1(bases[i],familyName)[0]
            lOvlp = Internal.getNodesFromName(F,'.Solver#Overlap')
            if lOvlp != []:
                Ovlp = lOvlp[0]
                N = Internal.getNodesFromName(F,'NeighbourList')[0]
                (parentBase, numFamily) = Internal.getParentOfNode(bases[i],F)
                (parentFamily, numOvlp) = Internal.getParentOfNode(F,Ovlp)
                (parentOvlp, numN) = Internal.getParentOfNode(Ovlp,N)
                listFamily = ''
                for l in xrange(len(fams[zi])):
                    listFamily = listFamily+fams[zi][l]
                    if l != len(fams[zi])-1:
                        listFamily = listFamily + ' '
                v = numpy.fromstring(listFamily, 'c')
                N[1]=v
                Ovlp[2][numN]=N
                F[2][numOvlp]=Ovlp
                bases[i][2][numFamily] = F
            # famille doubly-defined
            gc = Internal.getNodesFromType2(zones[zi],'GridConnectivity_t')
            for k in xrange(len(gc)):
                familyNameDD='FDD_'+ zones[zi][0]+'_'+gc[k][0]
                FDD = (Internal.getNodesFromName1(bases[i],familyNameDD))
                if FDD != []:
                    listVal="";domsDD=""
                    if isinstance(gc[k][1], numpy.ndarray):
                        val = gc[k][1].tostring()
                        listVal = val.split(",")
                    for vi in xrange(len(listVal)):
                        nvi = (Internal.getNodesFromName(tp,listVal[vi]))[0]
                        (pvi, nvi) = Internal.getParentOfNode(tp,nvi)
                        namevi = pvi[0]+'/F_'+listVal[vi]
                        if vi != len(listVal)-1:
                            namevi = namevi + ' '
                        domsDD = domsDD + namevi
                    OvlpDD  = Internal.getNodesFromName(FDD[0],'.Solver#Overlap')[0]
                    NDD  = Internal.getNodesFromName(FDD[0],'NeighbourList')[0]
                    (parentBaseDD, numFamilyDD) = Internal.getParentOfNode(bases[i],FDD[0])
                    (parentFamilyDD, numOvlpDD) = Internal.getParentOfNode(FDD[0],OvlpDD)
                    (parentOvlpDD, numNDD) = Internal.getParentOfNode(OvlpDD,NDD)
                    NDD[1]=numpy.fromstring(domsDD, 'c')
                    OvlpDD[2][numNDD]=NDD
                    FDD[0][2][numOvlpDD]=OvlpDD
                    bases[i][2][numFamilyDD] = FDD[0]

    cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv != []:
        (parent, pos) = Internal.getParentOfNode(tp,cgnsv[0])
        bases = [tp[2][pos]]+bases
    tp[2] = bases
    return tp

#==============================================================================
def addNeighbours(t, sameBase=0):
    """ Fill the NeighbourList nodes with bounding-box domains intersection.
    """
    tp = addNeighbours__(t)
    return tp

#==============================================================================
def addFamilyBCNode__(t):
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)
    for i in xrange(len(bases)):
        families = Internal.getNodesFromType1(bases[i], 'Family_t')
        for f in families:
            fbc = Internal.getNodesFromType1(f, 'FamilyBC_t')
            (parentBase, numFamily) = Internal.getParentOfNode(bases[i],f)
            if fbc == []:
                ud = numpy.fromstring('UserDefined', 'c')
                f[2].append(['FamilyBC', ud, [], 'FamilyBC_t'])
            else:
                fbc[0][0]='FamilyBC'
            bases[i][2][numFamily] = f

    cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv != []:
        (parent, pos) = Internal.getParentOfNode(tp, cgnsv[0])
        bases = [tp[2][pos]]+bases
    tp[2] = bases
    return tp

#==============================================================================
def addTurbulentDistanceIndex__(t):
    """Add the TurbulentDistance index node for elsA solver.
    """
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for j in xrange(len(zones)):
            (parentBase, numZone) = Internal.getParentOfNode(tp,zones[j])
            sol = Internal.getNodesFromName(zones[j], 'FlowSolution#Init')
            for k in xrange(len(sol)):
                (parentSol, numSol) = Internal.getParentOfNode(zones[j],sol[k])
                dist = Internal.getNodesFromName(zones[j], 'TurbulentDistance')
                if dist != []:
                    (parentDist, numDist) = Internal.getParentOfNode(zones[j],dist[0])
                    size = zones[j][2][numSol][2][numDist][1].shape
                    tdi = numpy.empty(size, numpy.int32); tdi.fill(-1)
                    zones[j][2][numSol][2].append(['TurbulentDistanceIndex', tdi, [], 'DataArray_t'])
            b[2][numZone] = zones[j]

    cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv != []:
        (parent, pos) = Internal.getParentOfNode(tp,cgnsv[0])
        bases = [tp[2][pos]]+bases
    tp[2] = bases
    return tp

#==============================================================================
def addTurbulentDistanceIndex(t):
    """Add the TurbulentDistance index node for elsA solver.
    """
    tp = addTurbulentDistanceIndex__(t)
    return tp

#==============================================================================
# def _addTurbulentDistanceIndex__(t):
#     zones = Internal.getZones(t)
#     for z in zones:
#         dist = Internal.getNodeFromName(z, 'TurbulentDistance')
#         if dist is not None:
#             size = dist[1].shape
#             tdi = numpy.empty(size, numpy.int32); tdi.fill(-1)
#             Internal.createChild(['TurbulentDistanceIndex', tdi, [], 'DataArray_t'])
#             b[2][numZone] = zones[j]

#     cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
#     if cgnsv != []:
#         (parent, pos) = Internal.getParentOfNode(tp,cgnsv[0])
#         bases = [tp[2][pos]]+bases
#     tp[2] = bases
#     return None


#==============================================================================
def buildMaskFiles__(t, maskInTree=True, keepOversetHoles=True, fileDir='.'):
    """ Build the mask files for elsA solver.
    """
    tp = Internal.copyRef(t)

    # dump des noeuds OversetHoles => construction des fichiers de masquage
    zones = Internal.getZones(tp)
    c = 0
    for z in zones:
        ho = Internal.getNodesFromType(z, 'OversetHoles_t')
        if ho != []:
            h = ho[0][2][1][1]
            dim = Internal.getZoneDim(z)
            h = Internal.convertIJKArray21DArray(h, dim[1]-1,dim[2]-1,dim[3]-1)
            hf = numpy.empty(h.shape, dtype=numpy.float64)
            hf[:] = h[:]
            array = ['cell_index', hf, hf.size, 1, 1]
            Converter.convertArrays2File([array], fileDir+'/hole_'+z[0]+'.v3d',
                                         'bin_v3d', dataFormat='%14.7e')
    if not keepOversetHoles: Internal._rmNodesByName(tp,'OversetHoles')
    return tp

#==============================================================================
def buildMaskFiles(t, maskInTree=True, keepOversetHoles=True, fileDir='.'):
    """ Build the mask files for elsA solver.
    """
    tp = buildMaskFiles__(t, maskInTree, keepOversetHoles, fileDir)
    return tp

#==============================================================================
def adaptNearmatch__(t):
    """Convert the nearmatch condition for use with elsA aerodynamic solver.
    """ 
    tp = Internal.copyRef(t)

    bases = Internal.getBases(tp)
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            (parentBase, numZone) = Internal.getParentOfNode(tp, z)
            gc = Internal.getNodesFromType(z, 'GridConnectivity_t')
            for l in xrange(len(gc)):
                nm = Internal.getNodesFromName2(gc[l], 'NMRatio')
                if nm != []:
                    gc[l][3]= 'GridConnectivity1to1_t'
                    transfo = Internal.getNodesFromName2(gc[l],'Transform')
                    prangedonor = Internal.getNodesFromName2(gc[l],'PointRangeDonor')
                    # copie des noeuds Transform et PointRangeDonor dans le noeud GridConnectivity1to1_t
                    if prangedonor != []:
                        gc[l][2].append(['PointRangeDonor', prangedonor[0][1], [], 'IndexRange_t'])
                    else:
                        print "Warning: adaptNearmatch__: PointRangeDonor missing in nearmatch join ",gc[l][0]
                    if transfo != []:
                        gc[l][2].append(['Transform', transfo[0][1], [], '\"int[IndexDimension]\"'])
                    else:
                        print "Warning: adaptNearmatch__: Transform missing in nearmatch join ",gc[l][0]
                    # on regarde si on est sur le raccord fin ou grossier (requis par elsAxdt)
                    fine = 0
                    iratio = 1; jratio = 1; kratio = 1
                    nmratioFact = nm[0][1][0]* nm[0][1][1]* nm[0][1][2]
                    if nmratioFact < 1:
                        fine = 1
                        iratio = int(math.ceil(1./nm[0][1][0]))
                        jratio = int(math.ceil(1./nm[0][1][1]))
                        kratio = int(math.ceil(1./nm[0][1][2]))
                    else:
                        iratio = int(nm[0][1][0])
                        jratio = int(nm[0][1][1])
                        kratio = int(nm[0][1][2])

                    gct = Internal.getNodesFromType1(gc[l],'GridConnectivityType_t')[0]
                    pr1,r1=Internal.getParentOfNode(gc,gct)
                    del pr1[2][r1]
                    gct = Internal.getNodesFromName1(gc[l],'PointListDonor')[0]
                    pr1,r1=Internal.getParentOfNode(gc,gct)
                    del pr1[2][r1]
                    gct = Internal.getNodesFromName1(gc[l],'UserDefinedData')[0]
                    pr1,r1=Internal.getParentOfNode(gc,gct)
                    del pr1[2][r1]
                    # Creation d'un noeud .Solver#Property contenant les noeuds fils :
                    #   jtype, type, matchside et nmratio
                    gc[l][2].append(['.Solver#Property', None, [], 'UserDefinedData_t'])
                    pos = len(gc[l][2])-1
                    jot = numpy.fromstring('nearmatch', 'c')
                    jo = numpy.fromstring('join', 'c')
                    if fine == 1: ms = numpy.fromstring('fine', 'c')
                    else: ms = numpy.fromstring('coarse', 'c')
                    gc[l][2][pos][2].append(['jtype', jot, [], 'DataArray_t'])
                    gc[l][2][pos][2].append(['type', jo, [], 'DataArray_t'])
                    gc[l][2][pos][2].append(['matchside', ms, [], 'DataArray_t'])
                    gc[l][2][pos][2].append(['i_ratio', iratio, [], 'DataArray_t'])
                    gc[l][2][pos][2].append(['j_ratio', jratio, [], 'DataArray_t'])
                    gc[l][2][pos][2].append(['k_ratio', kratio, [], 'DataArray_t'])

            b[2][numZone] = z

    cgnsv = Internal.getNodeFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv is not None:
        (parent, pos) = Internal.getParentOfNode(tp, cgnsv)
        bases = [tp[2][pos]]+bases
    tp[2] = bases

    return tp

#==============================================================================
def adaptNearmatch(t):
    """Convert the nearmatch condition for use with elsA aerodynamic solver.
    """ 
    tp = adaptNearmatch__(t)
    return tp

#===============================================================================
# clean='yes': delete obsolete Cassiopee nodes
# the current window to the
# In the CGNS Std (connectMatchPeriodic is compliant) :
# RotationAngle defines the angle from the current interface to the connecting interface
# Translation defines the vector from the current interface to the connecting interface
# In elsA: it is the opposite
#===============================================================================
def adaptPeriodicMatch__(t, clean='no'):
    """Convert the periodic match condition for use with elsA aerodynamic solver.
    """ 
    tp = Internal.copyRef(t)

    jointype = numpy.fromstring('join', 'c')
    jtopo = numpy.fromstring('periodic', 'c')
    ptype_rot = numpy.fromstring('rot', 'c')
    ptype_tra = numpy.fromstring('tra', 'c')
    jtype = numpy.fromstring('match', 'c')

    bases = Internal.getBases(tp)
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            dim = Internal.getZoneDim(z)
            ni = dim[1]; nj = dim[2]; nk = dim[3]; ninj = ni*nj
            connect = Internal.getNodesFromType(z, 'GridConnectivity1to1_t')
            for c in connect:
                donorName = Internal.getValue(c)
                zopp = Internal.getNodeFromName(tp, donorName)
                dimopp = Internal.getZoneDim(zopp)
                niopp = dimopp[1]; njopp = dimopp[2]; nkopp = dimopp[3]; ninjopp = niopp*njopp
                periodic = Internal.getNodesFromType(c, 'Periodic_t')
                if periodic != []:
                    p = periodic[0]
                    rotAngle = Internal.getNodesFromName(p, 'RotationAngle')
                    if rotAngle != []: rotAngle=rotAngle[0][1]
                    else: rotAngle = numpy.empty(0)
                    rotCenter = Internal.getNodesFromName(p, 'RotationCenter')
                    if rotCenter != []: rotCenter=rotCenter[0][1]
                    else: rotCenter = numpy.zeros(3)
                    translation = Internal.getNodesFromName(p, 'Translation')
                    if translation != []: translation=translation[0][1]
                    else: translation = numpy.empty(0)
                    Internal.createChild(c,'.Solver#Property','UserDefinedData_t',value=None, children=[])
                    solverProperty = c[2][len(c[2])-1]
                    Internal.createChild(solverProperty,'type'  ,'DataArray_t',value=jointype, children=[])
                    Internal.createChild(solverProperty,'jtopo' ,'DataArray_t',value=jtopo, children=[])
                    Internal.createChild(solverProperty,'jtype' ,'DataArray_t',value=jtype, children=[])
                    if rotAngle.any(): # at least one rotation angle not nul => periodicity by rotation
                        if len(numpy.where(rotAngle)[0]) != 1:
                            print "Warning: adaptPeriodicMatch__: rotation angle must have one non-zero component."
                            continue
                        axis = numpy.zeros(3)
                        axis[numpy.where(rotAngle)[0][0]]=1.0
                        angle = -rotAngle[numpy.where(rotAngle)[0][0]] # angle of the periodicity by rotation
                        angle1 = int(round(numpy.pi/angle))           # number of angular sectors of the component which includes the zone in radian
                        angle2 = 1                                    # number of channels defined by the component which includes the zone. Has to be set by user afterward
                        pointRange = Internal.getNodesFromName1(c,'PointRange')
                        pointRangeDonor = Internal.getNodesFromName1(c,'PointRangeDonor')
                        if pointRange != [] and pointRangeDonor != []:
                            pointRange = pointRange[0]; pointRangeDonor = pointRangeDonor[0]
                        else:
                            print "Warning: adaptPeriodicMatch__: missing PointRange(Donor) for join ",c[0],"."
                            continue
                        [wincur_imin, wincur_imax, wincur_jmin, wincur_jmax, wincur_kmin, wincur_kmax] = Internal.range2Window(pointRange[1])
                        [winopp_imin, winopp_imax, winopp_jmin, winopp_jmax, winopp_kmin, winopp_kmax] = Internal.range2Window(pointRangeDonor[1])
                        index_cur = ninj*(wincur_kmin-1)    + ni*(wincur_jmin-1)    + (wincur_imin-1) # Index of point Pcur
                        index_opp = ninjopp*(winopp_kmin-1) + niopp*(winopp_jmin-1) + (winopp_imin-1) # Index of point Popp
                        # Sign of the periodicity angle wrt the trigonometrical convention
                        if angle > 0.: pangle=1
                        else: pangle=-1

                        Internal.createChild(solverProperty,'ptype' ,'DataArray_t', value=ptype_rot)
                        Internal.createChild(solverProperty,'pangle','DataArray_t', value=pangle)
                        # angle and center are managed by node .Solver#Param, child of the zone
                        solverParam = Internal.getNodesFromName(z, '.Solver#Param')
                        if solverParam == []:
                            Internal.createChild(z,'.Solver#Param', 'UserDefinedData_t', value=None)
                            solverParam = z[2][len(z[2])-1]
                            Internal.createChild(solverParam,'axis_ang_1','DataArray_t',value=angle1)
                            Internal.createChild(solverParam,'axis_ang_2','DataArray_t',value=angle2)
                            Internal.createChild(solverParam,'axis_pnt_x','DataArray_t',value=rotCenter[0])
                            Internal.createChild(solverParam,'axis_pnt_y','DataArray_t',value=rotCenter[1])
                            Internal.createChild(solverParam,'axis_pnt_z','DataArray_t',value=rotCenter[2])
                            Internal.createChild(solverParam,'axis_vct_x','DataArray_t',value=axis[0])
                            Internal.createChild(solverParam,'axis_vct_y','DataArray_t',value=axis[1])
                            Internal.createChild(solverParam,'axis_vct_z','DataArray_t',value=axis[2])
                    else:
                        if translation.any(): # translation vector not nul => periodicity by translation
                            Internal.createChild(solverProperty,'ptype','DataArray_t',value=ptype_tra)
                            Internal.createChild(solverProperty,'xtran','DataArray_t',value= -translation[0])
                            Internal.createChild(solverProperty,'ytran','DataArray_t',value= -translation[1])
                            Internal.createChild(solverProperty,'ztran','DataArray_t',value= -translation[2])
                    if clean == 'yes':
                        c = Internal.rmNodes(c, 'Periodic_t')
    return tp

#==============================================================================
def adaptPeriodicMatch(t):
    """Convert the periodic match condition for use with elsA aerodynamic solver.
    """ 
    tp = adaptPeriodicMatch__(t)
    return tp

#==============================================================================
def addBaseToDonorZone__(t):
    tp = Internal.copyRef(t)

    zones = Internal.getZones(tp)
    for z in zones:
        subRegions= Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        interpSubRegions=[]
        for s in subRegions:
            sname = s[0]
            if sname.split('_')[0] == 'ID': interpSubRegions.append(s)
        for s in interpSubRegions:
            donorname = s[1].tostring()
            donorzone = Internal.getNodeFromName(tp, donorname)
            base,pos = Internal.getParentOfNode(tp, donorzone)
            s[1] = numpy.fromstring(base[0]+"/"+donorname,'c')

    return tp

#===============================================================================
# Remove useless families
def removeExtraFamily__(t, listOfNodes):
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)
    toDelete = []
    for i in xrange(len(bases)):
        baseName = bases[i][0]
        families = Internal.getNodesFromType1(bases[i], 'Family_t')
        for j in xrange(len(families)):
            familyName = families[j][0]
            path="%s/%s"%(baseName,familyName)
            # Checking if the family is not needed  if it is a default family with just a familyBC_t node of type UserDefined
            if path not in listOfNodes and len(families[j][2]) == 1 and families[j][2][0][2] == [] and families[j][2][0][1].tostring() == 'UserDefined':
                Internal._rmNode(tp,families[j])
                toDelete.append(path)

        zones = Internal.getNodesFromType1(bases[i], 'Zone_t')
        for j in xrange(len(zones)):
            familyNameNodes = Internal.getNodesFromType1(zones[j],'FamilyName_t')
            additionalFamilyNameNodes = Internal.getNodesFromType1(zones[j],'AdditionalFamilyName_t')

            for k in xrange(len(familyNameNodes)):
                familyName = familyNameNodes[k][1]
                if type(familyName) == 'numpy.ndarray':
                    familyName = familyName.tostring()
                familyPath = "%s/%s"%(baseName,familyName)
                if familyPath in toDelete:
                    Internal._rmNode(tp, familyNameNodes[k])

            for k in xrange(len(additionalFamilyNameNodes)):
                familyName = additionalFamilyNameNodes[k][1].tostring()
                familyPath = "%s/%s"%(baseName,familyName)
                if familyPath in toDelete:
                    Internal._rmNode(tp, additionalFamilyNameNodes[k])
    return tp

# Add new merged family and update family names
# - Add new family names
# - remove the old ones and get the .Solver#Overlap nodes
# - Replace value of familyName_t nodes with the new merged families for Zone_t and BC_t
# - Update NeighbourList node
def addMergedFamily__(t, equivalenceNodes):
    tp = Internal.copyRef(t)
    # Looking for overlap BC FamilySpecified!
    bases = Internal.getBases(tp)
    for i in xrange(len(bases)):
        baseName = bases[i][0]
        families = Internal.getNodesFromType1(bases[i], 'Family_t')

        # Adding new family
        for e in sorted(equivalenceNodes.keys()):
            if e.split("/")[0] == baseName:
                Fnew = [e.split("/")[1], None, [], 'Family_t']
                Fnew[2].append(['FamilyBC', 'UserDefined', [], 'FamilyBC_t'])
                ovlp = None
                for familyName in equivalenceNodes[e]:
                    fname = familyName.split("/")[1]
                    F = Internal.getNodeFromName(bases[i],fname)
                    if F is not None:
                        # Getting the .Solver#Overlap node
                        if len(Internal.getNodesFromName(F, '.Solver#Overlap')):
                            ovlp = Internal.getNodesFromName(F, '.Solver#Overlap')[0]
                    Internal._rmNode(tp, F)
                if ovlp:
                    Fnew[2].append(ovlp)
                bases[i][2].append(Fnew)

        # Updating NeighbourList
        for j in xrange(len(families)):
            familyName = families[j][0]
            path = "%s/%s"%(baseName,familyName)

            Ovlp  = Internal.getNodesFromName(families[j],'.Solver#Overlap')
            for k in xrange(len(Ovlp)):
                NeighbourList = Internal.getNodesFromName(Ovlp[k], 'NeighbourList')
                nl =  NeighbourList[0][1]
                if nl is not None:# can be None for doubly defined BCs
                    nl = nl.tostring().split()
                    for e in sorted(equivalenceNodes.keys()):
                        if equivalenceNodes[e] <= set(nl):
                            familyNodes = Internal.getNodesFromName(bases[i],e)
                            snl = set(nl) - equivalenceNodes[e]
                            snl.add(e)
                            NeighbourList[0][1] = numpy.fromstring(" ".join(list(snl)),'c')

        zones = Internal.getNodesFromType1(bases[i], 'Zone_t')
        for j in xrange(len(zones)):
            zoneName = zones[j][0]
            path = "%s/%s"%(baseName,zoneName)
            familyNameNodes = Internal.getNodesFromType1(zones[j],'FamilyName_t')
            additionalFamilyNameNodes = Internal.getNodesFromType1(zones[j],'AdditionalFamilyName_t')

            # Updating FamilyName_t for Zone_t
            for k in xrange(len(familyNameNodes)):
                familyName = familyNameNodes[k][1]
                if type(familyName) == 'numpy.ndarray':
                    familyName = familyName.tostring()
                familyPath = "%s/%s"%(baseName,familyName)
                for e in sorted(equivalenceNodes.keys()):
                    if familyPath in equivalenceNodes[e]:
                        familyNameNodes[k][1] = numpy.fromstring(e.split("/")[1],'c')

            # Updating AdditionalFamilyName_t for Zone_t
            for k in xrange(len(additionalFamilyNameNodes)):
                familyName = additionalFamilyNameNodes[k][1].tostring()
                familyPath = "%s/%s"%(baseName,familyName)
                for e in sorted(equivalenceNodes.keys()):
                    if familyPath in equivalenceNodes[e]:
                        additionalFamilyNameNodes[k][1] = numpy.fromstring(e.split("/")[1],'c')

            # Updating familyName_t for familySpeciefied BC_t
            zonebcs = Internal.getNodesFromType1(zones[j], 'ZoneBC_t')
            for k in xrange(len(zonebcs)):
                bcs = Internal.getNodesFromType1(zonebcs[k], 'BC_t')
                for l in xrange(len(bcs)):
                    if bcs[l][1].tostring() == 'FamilySpecified':
                        familyNameNodes = Internal.getNodesFromType1(bcs[l], 'FamilyName_t')
                        additionalFamilyNameNodes = Internal.getNodesFromType1(bcs[l], 'AdditionalFamilyName_t')
                        familyName = familyNameNodes[0][1].tostring()
                        familyPath = "%s/%s"%(baseName,familyName)
                        for e in sorted(equivalenceNodes.keys()):
                            if familyPath in equivalenceNodes[e]:
                                familyNameNodes[0][1] = numpy.fromstring(e.split("/")[1],'c')

                        # Updating AdditionalFamilyName_t for Zone_t
                        for k in xrange(len(additionalFamilyNameNodes)):
                            familyName = additionalFamilyNameNodes[k][1].tostring()
                            familyPath = "%s/%s"%(baseName,familyName)
                            for e in sorted(equivalenceNodes.keys()):
                                if familyPath in equivalenceNodes[e]:
                                    additionalFamilyNameNodes[k][1] = numpy.fromstring(e.split("/")[1],'c')


    return tp


# Building graph of neighbours zones
def buildGraph__(t):
    bases = Internal.getBases(t)
    g = {}
    for i in xrange(len(bases)):
        baseName = bases[i][0]
        families = Internal.getNodesFromType1(bases[i], 'Family_t')
        for j in xrange(len(families)):
            familyName = families[j][0]
            path = "%s/%s"%(baseName, familyName)
            Ovlp = Internal.getNodesFromName(families[j],'.Solver#Overlap')
            for k in xrange(len(Ovlp)):
                if path not in g.keys():
                    g[path] = set()
                NeighbourList = Internal.getNodesFromName(Ovlp[k], 'NeighbourList')
                nl = NeighbourList[0][1]
                if nl is not None: # can be None for doubly defined BCs
                    nl = nl.tostring().split()
                    g[path].update(nl)
                    for elt in nl:
                        if elt not in g.keys(): g[elt] = set()
                        g[elt].add(path)
    return g

# Compute equivalent families that can be merged.
# Families are equivalent and can be merged if they have the same neighbours
def buildPart__(g):
    equivalenceValues = {}
    equivalenceNodes = {}
    i = 0
    for k1 in sorted(g.keys()):
        for k2 in sorted(g.keys()):
            if k1 < k2:
                if g[k1] == g[k2]:
                    if g[k1] not in equivalenceValues.values():
                        i = i+1
                        equivalenceValues[i] = g[k1]
                    if i not in equivalenceNodes.keys():
                        equivalenceNodes[i] = set([])
                    equivalenceNodes[i].update([k1,k2])

    equivalenceNames = {}
    i = 0
    for e in equivalenceNodes:
        i += 1
        baseName = list(equivalenceNodes[e])[0].split("/")[0]
        newFamilyName = "F_Merged.%d"%(i)
        newPath = "%s/%s"%(baseName,newFamilyName)
        equivalenceNames[newPath] = equivalenceNodes[e]

    return equivalenceNames

# Nettoie l'arbre 
def cleanTree__(t):
    tp = Internal.copyRef(t)
    g = buildGraph__(tp)
    equivalenceNodes = buildPart__(g)
    neededFamily = g.keys()
    # Remove useless families
    tp = removeExtraFamily__(tp, g.keys())
    # Useless call!
    # Add new merged family and update family names
    tp = addMergedFamily__(tp, equivalenceNodes)
    return tp

#==============================================================================
# Conversion d'un arbre CGNS en un arbre de profile elsAxdt
# IN: t: arbre a convertir
# OUT: tp: arbre converti (copy)
#==============================================================================
def convert2elsAxdt(t, sameBase=0):
    """Convert a standard CGNS/python tree for use with elsA aerodynamic solver."""
    # Ajout noeud TurbulentDistanceIndex a -1, necessaire pour relire les
    # distances a la paroi.
    # Fonctionne pour du RANS sans transition.
    # Si transition ou si calcul LES, les valeurs des noeuds
    # TurbulentDistanceIndex sont utilisees
    tp = addTurbulentDistanceIndex__(t)

    # Creation de masques par fichier et destructions des noeuds OversetHoles
    # si on met un second argument a 'false'
    tp = buildMaskFiles__(tp)

    # Adaptation des raccords 'nearmatch' au profil elsAxdt
    tp = adaptNearmatch__(tp)

    # Adapt periodic match
    tp = adaptPeriodicMatch__(tp)

    # Construction de la liste "lZoneOverset" des zones de l'arbre contenant
    # un noeud 'Overset'
    tp = buildBCOverlap__(tp)

    # Copie des connectivites 'Overset' de Cassiopee dans des noeuds ZoneBC_t
    # de type BndOverlap
    tp = addNeighbours__(tp, sameBase)

    # BUG elsAxdt: ajout (ou eventuellement renommage) d'un noeud FamilyBC
    # aux familles qui n'en n'ont pas.
    tp = addFamilyBCNode__(tp)

    # Merge des familles pour alleger l arbre
    tp = cleanTree__(tp)

    # Interpolations : Prefixage nom de la zone donneuse par la base
    tp = addBaseToDonorZone__(tp)

    # Suppression des noeuds d'interpolation Cassiopee
    #tp = Internal.rmNodesByName(tp, 'ID_*')

    return tp

def _convert2elsAxdt(t, sameBase=0):
  # a faire!!
  return None

#=========================================================================================
def createElsaHybrid(t, method=0, axe2D=0):
    tp = Internal.copyRef(t)
    _createElsaHybrid(tp, method, axe2D)
    return tp

def _createElsaHybrid(t, method=0, axe2D=0):
    import converter
    zones = Internal.getZones(t)
    for z in zones:
         GEl = Internal.getElementNodes(z)
         NGON = 0; found = False
         for c in GEl:
             if c[1][0] == 22: found = True; break
             NGON += 1
         if found:
             node = GEl[NGON]
             CE = Internal.getNodeFromName1(node, 'ElementConnectivity')
             PE = Internal.getNodeFromName1(node, 'ParentElements')
             if PE is None:
                 Internal._adaptNFace2PE(z, remove=False)
                 PE = Internal.getNodeFromName1(node, 'ParentElements')
             er = Internal.getNodeFromName1(node, 'ElementRange')
             nfaces = er[1][1]-er[1][0]+1
             child = Internal.createUniqueChild(z, ':elsA#Hybrid', 'UserDefinedData_t')
             # to be removed (only used by elsA for nothing)
             sct = numpy.arange((nfaces), dtype=numpy.int32)
             Internal.newDataArray('SortedCrossTable', value=sct, parent=child)
             inct = numpy.empty((nfaces), dtype=numpy.int32)
             Internal.newDataArray('IndexNGONCrossTable', value=inct, parent=child)
             # OK
             ict = numpy.empty((nfaces), dtype=numpy.int32)
             bcct = numpy.empty((nfaces), dtype=numpy.int32)
             Internal.newDataArray('InversedCrossTable', value=ict, parent=child)
             Internal.newDataArray('BCCrossTable', value=bcct, parent=child)
             if axe2D > 0:
                 x = Internal.getNodeFromName2(z, 'CoordinateX')[1]
                 y = Internal.getNodeFromName2(z, 'CoordinateY')[1]
                 z = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
             else: x = None; y = None; z = None
             (iTRI, iQUADS, eTRI, eQUADS) = converter.createElsaHybrid(
                 CE[1], PE[1], ict, bcct, inct, method,
                 axe2D, x, y, z)
             if method == 0:
                 Internal.newDataArray('InternalTris', iTRI, parent=child)
                 Internal.newDataArray('InternalQuads', iQUADS, parent=child)
                 Internal.newDataArray('ExternalTris', eTRI, parent=child)
                 Internal.newDataArray('ExternalQuads', eQUADS, parent=child)
             else:
                 Internal.newDataArray('InternalElts', iTRI, parent=child)
                 Internal.newDataArray('ExternalElts', eTRI, parent=child)
         else:
             print 'Warning: createElsaHybrid: no NGON node found for zone %s. No :elsAHybrid node created.'%z[0]
    return None
