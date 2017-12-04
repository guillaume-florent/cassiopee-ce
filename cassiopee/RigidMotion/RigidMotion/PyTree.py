#
# Python Interface to compute/define rigid motion from PyTrees
#
import RigidMotion
__version__ = RigidMotion.__version__
 
try:
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Transform.PyTree as T
except:
    raise ImportError("RigidMotion: requires Converter, Transform modules.")

import numpy
import math

# Stocke les function deja definies
DEFINEDMOTIONS = {}

#=============================================================================
# Permet de definir un mouvement solide par des chaines (dependant de {t})
# definissant:
# la translation (tx,ty,tz), le centre de rotation (cx,cy,cz),
# et l'axe de rotation (ex,ey,ez)
#=============================================================================
def setPrescribedMotion1(t, name, tx="0", ty="0", tz="0",
                         cx="0", cy="0", cz="0",
                         ex="0", ey="0", ez="1", angle="0"):
    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp)
    for z in zones:
        # Recupere le conteneur TimeMotion
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is None:
            cont = ['TimeMotion', None, [], 'UseDefinedData_t']
            z[2].append(cont)

        # Le TimeRigidMotion name existe-t-il?
        motion = Internal.getNodeFromName1(cont, name)
        if motion is None:
            motion = [name, None, [], 'TimeRigidMotion_t']
            cont[2].append(motion)

        # Set it
        motion[2] = []
        motion[2].append(['MotionType', numpy.array([1], numpy.int32), [], 'DataArray_t'])
        motion[2].append(['tx', numpy.fromstring(tx, 'c'), [], 'DataArray_t'])
        motion[2].append(['ty', numpy.fromstring(ty, 'c'), [], 'DataArray_t'])
        motion[2].append(['tz', numpy.fromstring(tz, 'c'), [], 'DataArray_t'])
        motion[2].append(['cx', numpy.fromstring(cx, 'c'), [], 'DataArray_t'])
        motion[2].append(['cy', numpy.fromstring(cy, 'c'), [], 'DataArray_t'])
        motion[2].append(['cz', numpy.fromstring(cz, 'c'), [], 'DataArray_t'])
        motion[2].append(['ex', numpy.fromstring(ex, 'c'), [], 'DataArray_t'])
        motion[2].append(['ey', numpy.fromstring(ey, 'c'), [], 'DataArray_t'])
        motion[2].append(['ez', numpy.fromstring(ez, 'c'), [], 'DataArray_t'])
        motion[2].append(['angle', numpy.fromstring(angle, 'c'), [], 'DataArray_t'])
    return tp

#=============================================================================
# Permet de definir un mouvement de RotorMotion calcule par le CassiopeeSolver
#=============================================================================
def setPrescribedMotion2(
    t, name,
    transl_speed=(0.,0.,0.), psi0=0., psi0_b=0.,
    alp_pnt=(0.,0.,0.), alp_vct=(0.,1.,0.), alp0=0.,
    rot_pnt=(0.,0.,0.), rot_vct=(0.,0.,1.), rot_omg=0.,
    del_pnt=(0.,0.,0.), del_vct=(0.,0.,1.), del0=0.,
    delc=(0.,0.,0.), dels=(0.,0.,0.),
    bet_pnt=(0.,0.,0.), bet_vct=(0.,1.,0.), bet0=0.,
    betc=(0.,0.,0.), bets=(0.,0.,0.),
    tet_pnt=(0.,0.,0.), tet_vct=(1.,0.,0.), tet0=0.,
    tetc=(0.,), tets=(0.,),
    span_vct=(1.,0.,0.),
    pre_lag_pnt=(0.,0.,0.), pre_lag_vct=(1.,0.,0.), pre_lag_ang=0.,
    pre_con_pnt=(0.,0.,0.), pre_con_vct=(1.,0.,0.), pre_con_ang=0.):
    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp)
    for z in zones:
        # Recupere le conteneur TimeMotion
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is None:
            cont = ['TimeMotion', None, [], 'UseDefinedData_t']
            z[2].append(cont)

        # Le TimeRigidMotion name existe-t-il?
        motion = Internal.getNodeFromName1(cont, name)
        if motion is None:
            motion = [name, None, [], 'TimeRigidMotion_t']
            cont[2].append(motion)

        # Set it
        motion[2] = []
        motion[2].append(['MotionType', numpy.array([2], numpy.int32),
                          [], 'DataArray_t'])
        motion[2].append(['transl_speed', numpy.array([transl_speed[0], transl_speed[1], transl_speed[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['psi0', numpy.array([psi0], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['psi0_b', numpy.array([psi0_b], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['alp_pnt', numpy.array([alp_pnt[0], alp_pnt[1], alp_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['alp_vct', numpy.array([alp_vct[0], alp_vct[1], alp_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['alp0', numpy.array([alp0], numpy.float64),
                          [], 'DataArray_t'])
        
        motion[2].append(['rot_pnt', numpy.array([rot_pnt[0], rot_pnt[1], rot_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['rot_vct', numpy.array([rot_vct[0], rot_vct[1], rot_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['rot_omg', numpy.array([rot_omg], numpy.float64),
                          [], 'DataArray_t'])

        motion[2].append(['del_pnt', numpy.array([del_pnt[0], del_pnt[1], del_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['del_vct', numpy.array([del_vct[0], del_vct[1], del_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['del0', numpy.array([del0], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['delc', numpy.array([delc[0], delc[1], delc[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['dels', numpy.array([dels[0], dels[1], dels[2]], numpy.float64),
                          [], 'DataArray_t'])

        motion[2].append(['bet_pnt', numpy.array([bet_pnt[0], bet_pnt[1], bet_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['bet_vct', numpy.array([bet_vct[0], bet_vct[1], bet_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['bet0', numpy.array([bet0], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['betc', numpy.array([betc[0], betc[1], betc[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['bets', numpy.array([bets[0], bets[1], bets[2]], numpy.float64),
                          [], 'DataArray_t'])

        motion[2].append(['tet_pnt', numpy.array([tet_pnt[0], tet_pnt[1], tet_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['tet_vct', numpy.array([tet_vct[0], tet_vct[1], tet_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['tet0', numpy.array([tet0], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['tetc', numpy.array([tetc[0]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['tets', numpy.array([tets[0]], numpy.float64),
                          [], 'DataArray_t'])

        motion[2].append(['span_vct', numpy.array([span_vct[0], span_vct[1], span_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        
        motion[2].append(['pre_lag_pnt', numpy.array([pre_lag_pnt[0], pre_lag_pnt[1], pre_lag_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_lag_vct', numpy.array([pre_lag_vct[0], pre_lag_vct[1], pre_lag_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_lag_ang', numpy.array([pre_lag_ang], numpy.float64),
                          [], 'DataArray_t'])
        
        motion[2].append(['pre_con_pnt', numpy.array([pre_con_pnt[0], pre_con_pnt[1], pre_con_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_con_vct', numpy.array([pre_con_vct[0], pre_con_vct[1], pre_con_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_con_ang', numpy.array([pre_con_ang], numpy.float64),
                          [], 'DataArray_t'])
    return tp

#=============================================================================
# Permet de definir un mouvement de translation + rotation constant
#=============================================================================
def setPrescribedMotion3(
    t, name,
    transl_speed=(0.,0.,0.),
    axis_pnt=(0.,0.,0.), axis_vct=(0.,0.,0.),
    omega=0.):
    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp)
    for z in zones:
        # Recupere le conteneur TimeMotion
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is None:
            cont = ['TimeMotion', None, [], 'UseDefinedData_t']
            z[2].append(cont)

        # Le TimeRigidMotion name existe-t-il?
        motion = Internal.getNodeFromName1(cont, name)
        if motion is None:
            motion = [name, None, [], 'TimeRigidMotion_t']
            cont[2].append(motion)

        # Set it
        motion[2] = []
        motion[2].append(['MotionType', numpy.array([3], numpy.int32),
                          [], 'DataArray_t'])
        motion[2].append(['transl_speed', numpy.array([transl_speed[0], transl_speed[1], transl_speed[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['axis_pnt', numpy.array([axis_pnt[0], axis_pnt[1], axis_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['axis_vct', numpy.array([axis_vct[0], axis_vct[1], axis_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['omega', numpy.array([omega], numpy.float64),
                          [], 'DataArray_t'])
    return tp

#==============================================================================
# IN: m: motion node
# IN: string: chaine dependant de {t} a evaluer
# IN: time: instant
# OUT: valeur
#==============================================================================
def evalTimeString__(m, string, time):
    st = Internal.getNodeFromName1(m, string)
    st = st[1].tostring()
    st = st.replace('{t}', str(time))
    tx = eval(st)
    return tx

#==============================================================================
# IN: m: motion node
# IN: string: nom du noeud dont on retourne la valeur
#==============================================================================
def getNodeValue__(m, string):
    st = Internal.getNodeFromName1(m, string)
    return st[1]

#==============================================================================
def moveZone__(z, time):
    cont = Internal.getNodeFromName1(z, 'TimeMotion')
    if cont is not None:
        motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
        for m in motions:
            type = Internal.getNodeFromName1(m, 'MotionType')
            dtype = type[1][0]
            if (dtype == 1): # type 1: time string
                tx = evalTimeString__(m, 'tx', time)
                ty = evalTimeString__(m, 'ty', time)
                tz = evalTimeString__(m, 'tz', time)
                cx = evalTimeString__(m, 'cx', time)
                cy = evalTimeString__(m, 'cy', time)
                cz = evalTimeString__(m, 'cz', time)
                ex = evalTimeString__(m, 'ex', time)
                ey = evalTimeString__(m, 'ey', time)
                ez = evalTimeString__(m, 'ez', time)
                angle = evalTimeString__(m, 'angle', time)
                z = T.translate(z, (tx,ty,tz))
                if (angle != 0):
                    z = T.rotate(z, (cx,cy,cz), (ex-cx,ey-cy,ez-cz), angle)
            elif (dtype == 2): # type 2: rotation motion CassiopeeSolver
                try: import Cassiopee as K
                except: raise ImportError("evalPosition: motionRotor requires CassiopeeSolver.")
                import elsA_user as E # caveat
                import KBridge
                transl_speed = getNodeValue__(m, 'transl_speed')
                psi0 = getNodeValue__(m, 'psi0')
                psi0_b = getNodeValue__(m, 'psi0_b')
                alp_pnt = getNodeValue__(m, 'alp_pnt')
                alp_vct = getNodeValue__(m, 'alp_vct')
                alp0 = getNodeValue__(m, 'alp0')
                rot_pnt = getNodeValue__(m, 'rot_pnt')
                rot_vct = getNodeValue__(m, 'rot_vct')
                rot_omg = getNodeValue__(m, 'rot_omg')
                del_pnt = getNodeValue__(m, 'del_pnt')
                del_vct = getNodeValue__(m, 'del_vct')
                del0 = getNodeValue__(m, 'del0')
                delc = getNodeValue__(m, 'delc')
                dels = getNodeValue__(m, 'dels')
                bet_pnt = getNodeValue__(m, 'bet_pnt')
                bet_vct = getNodeValue__(m, 'bet_vct')
                bet0 = getNodeValue__(m, 'bet0')
                betc = getNodeValue__(m, 'betc')
                bets = getNodeValue__(m, 'bets')
                tet_pnt = getNodeValue__(m, 'tet_pnt')
                tet_vct = getNodeValue__(m, 'tet_vct')
                tet0 = getNodeValue__(m, 'tet0')
                tetc = getNodeValue__(m, 'tetc')
                tets = getNodeValue__(m, 'tets')
                span_vct = getNodeValue__(m, 'span_vct')
                pre_lag_ang = getNodeValue__(m, 'pre_lag_ang')
                pre_lag_pnt = getNodeValue__(m, 'pre_lag_pnt')
                pre_lag_vct = getNodeValue__(m, 'pre_lag_vct')
                pre_con_ang = getNodeValue__(m, 'pre_con_ang')
                pre_con_pnt = getNodeValue__(m, 'pre_con_pnt')
                pre_con_vct = getNodeValue__(m, 'pre_con_vct')

                if (DEFINEDMOTIONS.has_key(z[0]+'_'+m[0]) == False):
                    Func1 = E.function('rotor_motion', z[0]+'_'+m[0])
                    # angle initial du rotor par rapport au champ a l'infini
                    Func1.set("psi0", float(psi0[0]))
                    # Angle initial de la pale 
                    # (si z est l'axe de rotation cet angle vaut psi0_b+pi/2)
                    Func1.set("psi0_b", float(psi0_b[0]))
                    # parametres inclinaison du rotor
                    Func1.set("alp_pnt_x", float(alp_pnt[0]))
                    Func1.set("alp_pnt_y", float(alp_pnt[1]))
                    Func1.set("alp_pnt_z", float(alp_pnt[2]))
                    Func1.set("alp_vct_x", float(alp_vct[0]))
                    Func1.set("alp_vct_y", float(alp_vct[1]))
                    Func1.set("alp_vct_z", float(alp_vct[2]))
                    Func1.set("alp0", float(alp0[0]))
                    # parametres rotation uniforme du rotor
                    Func1.set("rot_pnt_x", float(rot_pnt[0]))
                    Func1.set("rot_pnt_y", float(rot_pnt[1]))
                    Func1.set("rot_pnt_z", float(rot_pnt[2]))
                    Func1.set("rot_vct_x", float(rot_vct[0]))
                    Func1.set("rot_vct_y", float(rot_vct[1]))
                    Func1.set("rot_vct_z", float(rot_vct[2]))
                    Func1.set("rot_omg", float(rot_omg[0]))
                    # parametres 1ere rotation: trainee
                    Func1.set("del_pnt_x", float(del_pnt[0]))
                    Func1.set("del_pnt_y", float(del_pnt[1]))
                    Func1.set("del_pnt_z", float(del_pnt[2]))
                    Func1.set("del_vct_x", float(del_vct[0]))
                    Func1.set("del_vct_y", float(del_vct[1]))
                    Func1.set("del_vct_z", float(del_vct[2]))
                    Func1.set("del0", float(del0[0]))
                    Func1.set("nhdel", 3)
                    Func1.set("del1c", float(delc[0]))
                    Func1.set("del1s", float(dels[0]))
                    Func1.set("del2c", float(delc[1]))
                    Func1.set("del2s", float(dels[1]))
                    Func1.set("del3c", float(delc[2]))
                    Func1.set("del3s", float(dels[2]))
                    # parametres 2eme rotation: battement
                    Func1.set("bet_pnt_x", float(bet_pnt[0]))
                    Func1.set("bet_pnt_y", float(bet_pnt[1]))
                    Func1.set("bet_pnt_z", float(bet_pnt[2]))
                    Func1.set("bet_vct_x", float(bet_vct[0]))
                    Func1.set("bet_vct_y", float(bet_vct[1]))
                    Func1.set("bet_vct_z", float(bet_vct[2]))
                    Func1.set("bet0", float(bet0[0]))
                    Func1.set("nhbet", 3)
                    Func1.set("bet1c", float(betc[0]))
                    Func1.set("bet1s", float(bets[0]))
                    Func1.set("bet2c", float(betc[1]))
                    Func1.set("bet2s", float(bets[1]))
                    Func1.set("bet3c", float(betc[2]))
                    Func1.set("bet3s", float(bets[2]))
                    # parametres 3eme rotation: pas cyclique 
                    Func1.set("tet_pnt_x", float(tet_pnt[0]))
                    Func1.set("tet_pnt_y", float(tet_pnt[1]))
                    Func1.set("tet_pnt_z", float(tet_pnt[2]))
                    Func1.set("tet_vct_x", float(tet_vct[0]))
                    Func1.set("tet_vct_y", float(tet_vct[1]))
                    Func1.set("tet_vct_z", float(tet_vct[2]))
                    Func1.set("tet0", float(tet0[0]))
                    Func1.set("nhtet", 1)
                    Func1.set("tet1c", float(tetc[0]))
                    Func1.set("tet1s", float(tets[0]))
                    Func1.set('span_vct_x', float(span_vct[0]))
                    Func1.set('span_vct_y', float(span_vct[1]))
                    Func1.set('span_vct_z', float(span_vct[2]))
                    Func1.set('pre_lag_ang', float(pre_lag_ang[0]))
                    Func1.set('pre_lag_pnt_x', float(pre_lag_pnt[0]))
                    Func1.set('pre_lag_pnt_y', float(pre_lag_pnt[1]))
                    Func1.set('pre_lag_pnt_z', float(pre_lag_pnt[2]))
                    Func1.set('pre_lag_vct_x', float(pre_lag_vct[0]))
                    Func1.set('pre_lag_vct_y', float(pre_lag_vct[1]))
                    Func1.set('pre_lag_vct_z', float(pre_lag_vct[2]))
                    Func1.set('pre_con_ang', float(pre_con_ang[0]))
                    Func1.set('pre_con_pnt_x', float(pre_con_pnt[0]))
                    Func1.set('pre_con_pnt_y', float(pre_con_pnt[1]))
                    Func1.set('pre_con_pnt_z', float(pre_con_pnt[2]))
                    Func1.set('pre_con_vct_x', float(pre_con_vct[0]))
                    Func1.set('pre_con_vct_y', float(pre_con_vct[1]))
                    Func1.set('pre_con_vct_z', float(pre_con_vct[2]))
                    DEFINEDMOTIONS[z[0]+'_'+m[0]] = Func1
                else: Func1 = DEFINEDMOTIONS[z[0]+'_'+m[0]]
                M = KBridge.evalKDesFunction(Func1, time)
                z = evalPosition___(z, None, M)
                z = T.translate(z, (transl_speed[0]*time,transl_speed[1]*time,transl_speed[2]*time))
            elif (dtype == 3): # type 3: constant transl + rotation
                transl_speed = getNodeValue__(m, 'transl_speed')
                axis_pnt = getNodeValue__(m, 'axis_pnt')
                axis_vct = getNodeValue__(m, 'axis_vct')
                omega = getNodeValue__(m, 'omega')
                tx = transl_speed[0]*time
                ty = transl_speed[1]*time
                tz = transl_speed[2]*time
                z = T.translate(z, (tx,ty,tz))
                cx = axis_pnt[0]+tx
                cy = axis_pnt[1]+ty
                cz = axis_pnt[2]+tz
                ex = axis_vct[0]
                ey = axis_vct[1]
                ez = axis_vct[2]
                angle = omega[0]*time*180./math.pi
                z = T.rotate(z, (cx,cy,cz), (ex-cx,ey-cy,ez-cz), angle)
    return z
    
#==============================================================================
# Evalue la position reelle de la zone a l'instant t
# Le mouvement est stockee dans chaque zone.
# Les coordonnees de la zone sont modifiees
#==============================================================================
def evalPosition__(a, time):
    tp = Internal.copyRef(a)
    type = Internal.typeOfNode(tp)
    if (type == 1): return moveZone__(tp, time)
    elif (type == 2):
        for i in xrange(len(tp)): tp[i] = moveZone__(tp[i], time)
    else:
        bases = Internal.getNodesFromType(tp, 'CGNSBase_t')
        for b in bases:
            c = 0
            for n in b[2]:
                if n[3] == 'Zone_t': b[2][c] = moveZone__(n, time)
                c += 1
    return tp

#==============================================================================
# Evalue la position reelle de la zone a l'instant t
# Le mouvement est definie dans la fonction F.
# Les coordonnees de la zone sont modifiees
#==============================================================================
def evalPosition___(a, time, F):
    """Move the mesh with defined motion to time t. Return an array with
    moved mesh coordinates.
    Usage: evalPosition2(a, time, F)"""
    return C.TZA(a, 'nodes', 'nodes', RigidMotion.evalPosition, None, time, F)

#==============================================================================
# Evalue la position reelle de la zone a l'instant t
#==============================================================================
def evalPosition(a, time, F=None):
    if F is None: return evalPosition__(a, time)
    else: return evalPosition___(a, time, F)
