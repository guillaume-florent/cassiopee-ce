"""Conversion module for Cassiopee package.
"""
__version__ = '2.5'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud, Benoit Rodriguez, Simon Verley, Bruno Maugars, Thomas Renaud"
#
# Python Interface for conversion between  array / file / CGNS
#
import numpy
import converter
import KCore

# -- Create an array --
# Les champs sont mis a zero, sauf si pour les champs cellN et cellNF
def array(vars, n1, n2, sub):
    """Create a structured or unstructured array.
    Usage: array(vars, ni, nj, nk)
    or array(vars, npoints, nelts, eltType)"""
    if isinstance(sub, str):
        return arrayNS(vars, n1, n2, sub)
    else:
        return arrayS(vars, n1, n2, sub)

# -- Create a structured array --
def arrayS(vars, ni, nj, nk):
    """Create a structured array.
    Usage: array(vars, ni, nj, nk)"""
    ni = int(ni); nj = int(nj); nk = int(nk)
    vars = vars.replace(' ','')
    l = len(vars)
    if (vars[l-1] == ',' or vars[0] == ','):
        print "Warning: array: your var string is suspicious."
    vl = vars.split(','); v = len(vl)
    a = numpy.zeros((v, ni*nj*nk), numpy.float64)
    for i in xrange(v):
        if (vl[i] == 'cellN' or vl[i] == 'cellNF'): a[i,:] = 1.
    return [vars, a, ni, nj, nk]

# -- Create an unstructured array --
def arrayNS(vars, npoints, nelts, eltType, sizeFN=0, sizeEF=0):
    """Create a unstructured array.
    Usage: array(vars, npoints, nelts, eltType)"""
    npoints = int(npoints); nelts = int(nelts)
    vars = vars.replace(' ','')
    l = len(vars)
    if (vars[l-1] == ',' or vars[0] == ','):
        print "Warning: array: your var string is suspicious."
    vl = vars.split(','); v = len(vl)
    if eltType == 'NODE' or eltType == 'NODE*': nt = 1
    elif eltType == 'BAR' or eltType == 'BAR*': nt = 2
    elif eltType == 'TRI' or eltType == 'TRI*': nt = 3
    elif eltType == 'QUAD' or eltType == 'QUAD*': nt = 4
    elif eltType == 'TETRA' or eltType == 'TETRA*': nt = 4
    elif eltType == 'PYRA' or eltType == 'PYRA*': nt = 5
    elif eltType == 'PENTA' or eltType == 'PENTA*': nt = 6
    elif eltType == 'HEXA' or eltType == 'HEXA*': nt = 8
    elif eltType == 'NGON' or eltType == 'NGON*':
        raise ValueError("array: this function doesnt work for NGONs.")
    else:
        raise ValueError("array: wrong element type: %s."%eltType)

    if eltType[len(eltType)-1] == '*':
        a = numpy.zeros((v, nelts), numpy.float64)
    else:
        a = numpy.zeros((v, npoints), numpy.float64)
    for i in xrange(v):
        if vl[i] == 'cellN' or vl[i] == 'cellNF': a[i,:] = 1.

    c = numpy.ones((nt, nelts), numpy.int32)
    return [vars, a, c, eltType]

# -- Get the number of points of an array --
def getNPts(a):
    """Return the total number of points."""
    if isinstance(a[0], list):
        b = 0
        for i in a: b += getNPts__(i)
        return b
    else: return getNPts__(a)

def getNPts__(a):
    npts = 0
    if len(a) == 5: # structure
        npts = a[2]*a[3]*a[4]
    elif len(a) == 4: # non structure
        if a[3][-1] != '*':
            npts = a[1].shape[1]
        else: npts = 0 # don't know
    return npts

# -- Get the number of cells of an array --
def getNCells(a):
    """Return the total number of cells."""
    if isinstance(a[0], list):
        b = 0
        for i in a: b += getNCells__(i)
        return b
    else: return getNCells__(a)

def getNCells__(a):
    ncells = 0
    if len(a) == 5: # structure
        nic = a[2]-1; njc = a[3]-1; nkc = a[4]-1
        nic = max(nic, 1); njc = max(njc, 1); nkc = max(nkc, 1)
        ncells = nic*njc*nkc
    elif len(a) == 4: # non structure
        if (a[3] == 'NGON' or a[3] == 'NGON*'):
            offset = a[2][0,1]
            ncells = a[2][0,offset+2]
        else: ncells = a[2].shape[1]
    return ncells

# -- Copy an array --
def copy(array):
    """Copy  an  array.
    Usage: a2 = copy(a1)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.copy(i))
        return b
    else:
        return converter.copy(array)

def setPartialFields(a, ad, listIndices):
    """Set values defined in a 1D array for points in listIndices
    Usage: b=setPartialFields(a,ad,listIndices)"""
    return converter.setPartialFields(a,ad,listIndices)

# -- rmVars --
def rmVars(a, var):
    """Remove variables."""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(rmVars__(i, var))
        return b
    else:
        return rmVars__(a, var)

def rmVars__(a, var):
    varstring = a[0]
    vs = varstring.split(',')
    eVars = []
    if isinstance(var, list):
        for v in vs:
            present = False
            for w in var:
                if v == w:
                    present = True; break
            if not present: eVars.append(v)
    else:
        for v in vs:
            if v != var: eVars.append(v)
    return extractVars(a, eVars)

# -- Extract a list of variables from an array --
def extractVars(array, vars):
    """Extract variables from a.
    Usage: a = extractVars(array, ['x','ro'])"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            varsi = []
            if isinstance(vars, list):
                for j in vars:
                    p = KCore.isNamePresent(i, j)
                    varsi.append(p+1)
            else:
                p = KCore.isNamePresent(i, vars)
                varsi.append(p+1)
            b.append(converter.extractVars(i, varsi))
        return b
    else:
        varsi = []
        for i in vars:
            p = KCore.isNamePresent(array, i)
            varsi.append(p+1)
        return converter.extractVars(array, varsi)

# -- Add variables --
def addVars(array, add=''):
    """Add variables to an array.
    Usage: a = addVars(array, 'varString')
    or: a = addVars([a1, a2, a3])"""
    if add != '':
        if isinstance(array[0], list):
            b = []
            for i in array: b.append(addVar(i, add))
            return b
        else:
            return addVar(array, add)
    else:
        if isinstance(array[0][0], list):
            l = 0
            for i in array:
                if l == 0: l = len(i)
                elif len(i) != l:
                    raise ValueError("addVars: can not add variables.")
            c = 0
            b = []
            for i in array[0]:
                l = []
                for i in array: l.append(i[c])
                c += 1
                b.append(addVars1(l))
            return b
        else:
            return addVars1(array)

# -- Add variables defined by a string --
def addVar(array, vars):
    """Add variables to an array.
    Usage: a = addVar(array, vars)"""
    if isinstance(vars, list):
        a = array
        for i in vars: a = converter.addVar(a, i)
        return a
    else:
        return converter.addVar(array, vars)

# -- Add variables defined by arrays --
def addVars1(arrays):
    """Put all variables in an array. arrays must be a list of arrays
    defined on the same grid.
    Usage: a = addVars(arrays)"""
    return converter.addVars(arrays)

# -- randomize a variable --
def randomizeVar(array, var, deltaMin, deltaMax):
    """Randomize a field defined by var within a range [a-deltaMin, a+deltaMax]
    Usage: a = randomizeVar(array, var, deltaMin, deltaMax)"""
    if isinstance(array[0], list):
        out = []
        for i in array: out.append(converter.randomizeVar(i, var, deltaMin, deltaMax))
        return out
    else:
        return converter.randomizeVar(array, var, deltaMin, deltaMax)

# -- Init variables --
def initVars(array, varNameString, val1=[], val2=[]):
    """Init variables defined by VarNameString.
    Usage: a = initVars(array, varNameString, val)
    or: a = initVars(array, varNameString, F, [strings])"""
    if val1 == []:
        if isinstance(array[0], list):
            b = []
            for i in array:
                b.append(initVarByEq(i, varNameString))
            return b
        else: return initVarByEq(array, varNameString)

    elif callable(val1):
        if isinstance(array[0], list):
            b = []
            for i in array:
                posvar = KCore.isNamePresent(i, varNameString)
                if posvar == -1: i = addVar(i, varNameString)
                b.append(initVar(i, varNameString, val1, val2))
            return b
        else:
            posvar = KCore.isNamePresent(array, varNameString)
            if posvar == -1:
                array = addVar(array, varNameString)
            return initVar(array, varNameString, val1, val2)
    else:
        if isinstance(array[0], list):
            b = []
            for i in array:
                posvar = KCore.isNamePresent(i, varNameString)
                if posvar == -1:
                    i = addVar(i, varNameString)
                b.append(converter.initVars(i, varNameString, val1))
            return b
        else:
            posvar = KCore.isNamePresent(array, varNameString)
            if posvar == -1:
                array = addVar(array, varNameString)
            return converter.initVars(array, varNameString, val1)

# -- Init variable by eq --
def initVarByEq(array, eq):
    a = copy(array)
    # Extrait les variables de a
    varstring = a[0]
    vars = varstring.split(',')

    eq = eq.replace('centers:', '')

    # Split suivant ; si plusieurs formules sont definies
    eq = eq.split(';')

    from numpy import *
    for eq0 in eq:
        # Extrait la variable a initialiser de eq
        s = eq0.split('=')
        if len(s) != 2:
            print 'Error: initVars: equation is incorrect.'; return a

        var = s[0]; var = var.replace('{', ''); var = var.replace('}', '')
        var = var.lstrip(); var = var.rstrip()
        varp = KCore.isNamePresent(a, var)
        if varp == -1:
            a = addVars(a, var); varp = KCore.isNamePresent(a, var)

        # Initialisation de la variable
        ap = a[1]
        loc = s[1]
        c = 0
        for v in vars:
            loc = loc.replace('{%s}'%v, 'ap[%d,:]'%c); c += 1

        # Evaluation
        ap[varp,:] = eval(loc)
    return a

# -- Init variable with a function --
def initVar(array, varNameString, F, strings):
    """Init variables to the value val in array containing
    their name in the array.
    Usage: a = initVar(array, varNameString, F, strings)"""
    b = copy(array)
    posvar = KCore.isNamePresent(b, varNameString)+1

    pos = []
    for i in strings:
        p = KCore.isNamePresent(b, i)+1
        if p == 0:
            raise TypeError("initVar: can't find %s in array."%i)
        else: pos.append(p)

    n = b[1]
    nsize = n.shape[1]
    l = len(strings)

    if l == 0:
        for i in xrange(nsize):
            n[posvar-1,i] = F()
    else:
        for i in xrange(nsize):
            x = [ n[pos[j]-1,i] for j in xrange(l)]
            n[posvar-1,i] = F(*x)
    return b

# Converti l'extension en nom de format
def convertExt2Format__(fileName):
    """Convertit un fichier en format suivant son extension."""
    import os.path
    ext = os.path.splitext(fileName)
    extension = ext[len(ext)-1]; extension = extension.lower()
    if extension == '.plt': format = 'bin_tp'
    elif extension == '.dat' or extension == '.tp': format = 'fmt_tp'
    elif extension == '.v3d': format = 'bin_v3d'
    elif extension == '.fv3d': format = 'fmt_v3d'
    elif extension == '.mesh': format = 'fmt_mesh'
    elif extension == '.msh': format = 'fmt_gmsh'
    elif extension == '.stl': format = 'fmt_stl'
    elif extension == '.fstl': format = 'fmt_stl'
    elif extension == '.bstl': format = 'bin_stl'
    elif extension == '.fig': format = 'fmt_xfig'
    elif extension == '.svg': format = 'fmt_svg'
    elif extension == '.pov': format = 'fmt_pov'
    elif extension == '.cgns': format = 'bin_cgns'
    elif extension == '.adf': format = 'bin_cgns'
    elif extension == '.hdf': format = 'bin_hdf'
    elif extension == '.hdf5': format = 'bin_hdf'
    elif extension == '.pickle': format = 'bin_pickle'
    elif extension == '.df3': format = 'bin_df3'
    elif extension == '.3ds': format = 'bin_3ds'
    elif extension == '.ply': format = 'bin_ply'
    elif extension == '.obj': format = 'fmt_obj'
    elif extension == '.gts': format = 'fmt_gts'
    elif extension == '.png': format = 'bin_png'
    elif extension == '.d': format = 'fmt_cedre'
    elif extension == '.su2': format = 'fmt_su2'
    elif extension == '.gbin': format = 'bin_plot3d'
    elif extension == '.gfmt': format = 'fmt_plot3d'
    elif extension == '.iges' or extension == '.igs': format = 'fmt_iges'
    elif extension[0:4] == '.ref': format = 'bin_pickle'
    else: format = 'unknown'
    return format

def convertFile2Arrays(fileName, format=None, nptsCurve=20, nptsLine=2,
                       density=-1., zoneNames=None, BCFaces=None):
    """Read file and return arrays containing file data.
    Usage: a = convertFile2Arrays(fileName, options)"""
    try: import locale; locale.setlocale(locale.LC_NUMERIC, 'C') # force .
    except: pass
    if format is None: format = convertExt2Format__(fileName); autoTry = True
    else: autoTry = False
    try: file = open(fileName, 'r')
    except: raise IOError("convertFile2Arrays: file %s not found."%fileName)
    file.close()
    if format == 'bin_pickle':
        import cPickle as pickle
        print 'Reading \''+fileName+'\'...',
        try:
            file = open(fileName, 'rb')
            a = pickle.load(file)
            file.close()
        except:
            raise TypeError("convertFile2Arrays: file %s can not be read."%fileName)
        else:
            print 'done.'
            return a
    elif format == 'fmt_iges':
        try: import OCC
        except: raise ImportError("convertFile2Arrays: IGES reader requires OCC module.")
        a = OCC.convertIGES2Arrays(fileName, h=0., chordal_err=0.)
        for c in xrange(len(a)): zoneNames.append('zone%d'%c)
        return a
    elif format == 'fmt_free':
        print 'Reading '+fileName+' (fmt_free)...',
        try:
            file = open(fileName, 'r')
            f = file.read()
            file.close()
            f = f.split('\n')
            np = len(f) # nbre d'enregistrements
            if f[np-1] == '' or f[np-1] == ' ': np = np-1
            nv = 1 # nbre de variables
            if np > 2: line = f[2]
            elif np > 1: line = f[1]
            elif np > 0: line = f[0]
            line = line.split(' ')
            nv = len(line)
            # varstring
            varString = ''
            for i in xrange(nv-1): varString += 'v%d,'%(i+1)
            if nv > 0: varString += 'v%d'%nv
            a = array(varString, np, 1, 1)
            pt = a[1]
            for i in xrange(np):
                line = f[i]; line = line.split(' ')
                for n in xrange(len(line)): pt[n,i] = float(line[n])
            print 'done.'
            zoneNames.append('zone0')
            return [a]
        except:
            raise TypeError("convertFile2Arrays: file %s can not be read."%fileName)
    else:
        try:
            return  converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine, density, zoneNames, BCFaces)
        except:
            if not autoTry: raise
            else: pass
            format = checkFileType(fileName)
            try: 
               return  converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine, density, zoneNames, BCFaces)
            except:   
                FORMATS = ['bin_ply', 'fmt_tp', 'fmt_v3d',
                'bin_tp', 'bin_v3d', 'fmt_mesh',
                'fmt_gmsh', 'bin_gmsh', 'fmt_stl', 'bin_stl',
                'fmt_xfig', 'fmt_svg', 'bin_3ds',
                'fmt_obj', 'fmt_gts' , 'fmt_pov']
                success = 0
                for fmt in FORMATS:
                    try:
                        a = converter.convertFile2Arrays(fileName, fmt, nptsCurve, nptsLine, density, zoneNames, BCFaces)
                        success = 1; break
                    except: pass
                    if success == 1: return a
                    else: return converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine, density, zoneNames, BCFaces)

def convertArrays2File(arrays, fileName, format=None, isize=4, rsize=8,
                       endian='big', colormap=0, dataFormat='%.9e ',
                       zoneNames=[], BCFaces=[]):
    """Write arrays to output file.
    Usage: convertArrays2File(arrays, fileName, format, options)"""
    try: import locale; locale.setlocale(locale.LC_NUMERIC, 'C') # force .
    except: pass
    if arrays != [] and not isinstance(arrays[0], list): arrays = [arrays]
    znames = []
    if len(zoneNames) == 0:
        for iz in xrange(len(arrays)): znames.append('Zone%d'%(iz+1))
    else:
        znames = zoneNames
        if len(zoneNames) != len(arrays):
            raise ValueError("convertArrays2File: zoneNames list, if not empty, must have the same length of arrays list.")
    if format is None:
        format = convertExt2Format__(fileName)
    if format == 'bin_pickle':
        import cPickle as pickle
        file = open(fileName, 'wb')
        print 'Writing \''+fileName+'\'...',
        pickle.dump(arrays, file, protocol=pickle.HIGHEST_PROTOCOL); file.close()
        print 'done.'
    elif format == 'fmt_free':
        file = open(fileName, 'w')
        print 'Writing %s (fmt_free)...'%fileName,
        for a in arrays:
            out = ''
            pt = a[1]
            nv = pt.shape[0]
            s = pt.shape[1]
            for i in xrange(s):
                for n in xrange(nv-1):
                    out += dataFormat%(pt[n,i])
                    out += ' '
                out += dataFormat%(pt[nv-1,i])
                out += '\n'
            file.write(out)
        file.close()
        print 'done.'
    else:
        converter.convertArrays2File(arrays, fileName, format, isize, rsize,
                                     endian, colormap, dataFormat,
                                     znames, BCFaces)

def diffArrays(arrays1, arrays2, arrays3=[]):
    """Diff arrays defining solutions.
    Usage: diffArrays(arrays1, arrays2, [arrays3])"""
    return converter.diffArrays(arrays1, arrays2, arrays3)

def getValue(array, ind):
    """Return the values of an array for a point of index ind or (i,j,k)...
    Usage: getValue(array, ind)"""
    if isinstance(array[0], list):
        raise TypeError("getValue: only for one array.")

    if isinstance(ind, tuple):
        if len(array) != 5: # structure
            raise TypeError("getValue: (i,j,k) indexing is only for structured array.")
        ni = array[2]; nj = array[3]
        if len(ind) == 3: index = (ind[0]-1)+(ind[1]-1)*ni+(ind[2]-1)*ni*nj
        elif len(ind) == 2: index = (ind[0]-1)+(ind[1]-1)*ni
        else: raise ValueError("getValue: too much values in index tuple.")
    else: index = ind

    n = array[1]; nfld = n.shape[0]
    v = []
    for nf in xrange(nfld): v.append(n[nf, index])
    return v

def setValue(array, ind, values):
    """Set the values in an array for a point of index ind or (i,j,k)...
    Usage: setValue(array, ind, values)"""
    if isinstance(array[0], list):
        raise TypeError("setValue: only for one array.")

    if isinstance(ind, tuple):
        if len(array) != 5: # structure
            raise TypeError("setValue: (i,j,k) indexing is only for structured array.")
        ni = array[2]; nj = array[3]
        if len(ind) == 3:
            index = (ind[0]-1)+(ind[1]-1)*ni+(ind[2]-1)*ni*nj
        elif len(ind) == 2:
            index = (ind[0]-1)+(ind[1]-1)*ni
        else:
            raise ValueError("setValue: too much values in index tuple.")
    else:
        index = ind
    ar = array[1]
    nf = ar.shape[0]
    nf2 = len(values)
    if nf2 != nf: raise ValueError("setValue: values is badly dimensioned.")
    for i in xrange(nf): ar[i, index] = values[i]
    return None

def getArgMin(array, varName):
    """Get array value where the variable defined by varName is minimum.
    Usage: getArgMin(array, varName)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = converter.getArgMin(array[0], varName)
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getArgMin(i, varName)
            if a[pos2] < b[pos1]:
                b = a
                pos1 = pos2
        return b
    else:
        return converter.getArgMin(array, varName)

def getMinValue(array, varName):
    """Get the minimum value of variable defined by varName in array.
    Usage: getMinValue(array, varName)"""
    if not isinstance(varName, list): varNames = [varName]
    else: varNames = varName
    out = []
    for v in varNames:
        if isinstance(array[0], list): pos = KCore.isNamePresent(array[0], v)
        else: pos = KCore.isNamePresent(array, v)
        out.append(getArgMin(array, varName)[pos])
    if len(out) == 1: return out[0]
    else: return out

def getArgMax(array, varName):
    """Get array value where the variable defined by varName is maximum.
    Usage: getArgMax(array, varName)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = converter.getArgMax(array[0], varName)
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getArgMax(i, varName)
            if a[pos2] > b[pos1]:
                b = a
                pos1 = pos2
        return b
    else:
        return converter.getArgMax(array, varName)

def getMaxValue(array, varName):
    """Get the maximum value of variable defined by varName in array.
    Usage: getMaxValue(array, varName)"""
    if not isinstance(varName, list): varNames = [varName]
    else: varNames = varName
    out = []
    for v in varNames:
        if isinstance(array[0], list): pos = KCore.isNamePresent(array[0], v)
        else: pos = KCore.isNamePresent(array, v)
        out.append(getArgMax(array, varName)[pos])
    if len(out) == 1: return out[0]
    else: return out

def getMeanValue(array, varName):
    """Get the mean value of the variable defined by varName in an array.
    Usage: getMeanValue(array, varName)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = 0.
        lb = 0.
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getMeanValue(i, varName)
            la = i[1].shape[1]*1.
            b = lb*b + la*a
            lb = lb + la
            b = b / lb
        return b
    else:
        return converter.getMeanValue(array, varName)

def getMeanRangeValue(array, varName, rmin, rmax):
    """Get the mean value of the variable defined by varName for a sorted range in an array.
    Usage: getMeanRangeValue(array, varName, rmin, rmax)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = 0.
        lb = 0.
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getMeanRangeValue(i, varName, rmin, rmax)
            la = i[1].shape[1]*1.*(rmax-rmin)
            b = lb*b + la*a
            lb = lb + la
            b = b / lb
        return b
    else:
        return converter.getMeanRangeValue(array, varName, rmin, rmax)

def normL0(array, varName):
    """Get the L0 norm of the field defined by varName in the array.
    If celln exists in the array, the norm for blanked points is not computed.
    Usage: normL0(array, varName)"""
    if isinstance(array[0], list):
        norm = 0
        for i in array:
            norm = max(converter.normL0(i, varName), norm)
        return norm
    else:
        return converter.normL0(array, varName)

def normL2(array, varName):
    """Get the L2 norm of the field defined by varName in the array.
    If celln exists in the array, the norm for blanked points is not computed.
    Usage: normL2(array, varName)"""
    if isinstance(array[0], list):
        norm = 0.; nsize = 0
        for i in array:
            nsize0 = i[1].shape[1]
            nsize += nsize0
            norm0 = converter.normL2(i, varName)
            norm += norm0*norm0*nsize0
        if nsize > 0: norm = numpy.sqrt(norm/nsize)
        return norm
    else:
        return converter.normL2(array, varName)

def normalize(array, vars):
    """Get the normalisation of the fields defined by vars in the array.
    Usage: normalize(array, vars)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.normalize(i, vars))
        return b
    else:
        return converter.normalize(array, vars)

def magnitude(array, vars):
    """Get the magnitude of the fields defined by vars in the array.
    Usage: magnitude(array, vars)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.magnitude(i, vars))
        return b
    else:
        return converter.magnitude(array, vars)

def convertHexa2Struct(array):
    """Convert an hexa array in a struct array.
    Usage: convertHexa2Struct(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.convertHexa2Struct(i))
        return b
    else:
        return converter.convertHexa2Struct(array)

# Convert a BAR array without ramifications, closed into an i-array
def convertBAR2Struct(array):
    """Convert a BAR array without ramifications, closed into an i-array.
    Usage: convertBAR2Struct(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            if len(i) == 4: b.append(converter.convertBAR2Struct(i))
            else: b.append(i)
        return b
    else:
        if len(array) == 4: return converter.convertBAR2Struct(array)
        else: return array

def convertTri2Quad(array, alpha=30.):
    """Convert a TRI array to a QUAD array.
    Usage: convertTri2Quad(array, alpha)"""
    if isinstance(array[0], list):
        b = []; c = []
        for i in array:
            res = converter.convertTri2Quad(i, alpha)
            b += res[0]; c += res[1]
        return b, c
    else:
        return converter.convertTri2Quad(array, alpha)

def conformizeNGon(array, tol=1.e-6):
    """Conformize topologically a NGON array.
    Usage: conformizeNGon(array, tol)"""
    if len(array) == 4 and array[3] == 'NGON':
        return converter.conformizeNGon(array, tol)
    else: return array


# -- interne --
def convertArray2Tetra1__(array, arrayC=[], split='simple'):
    try: sub = array[3]
    except: raise TypeError("convertArray2Tetra: arg must be an array.")

    if isinstance(sub, str): t = sub
    else: t = 'STRUCT'

    if split == 'simple': # no points added
        if t == 'STRUCT':
            return converter.convertStruct2Tetra(array)
        elif t == 'HEXA' or t == 'QUAD':
            return converter.convertHexa2Tetra(array)
        elif t == 'PENTA':
            return converter.convertPenta2Tetra(array)
        elif t == 'PYRA':
            return converter.convertPyra2Tetra(array)
        elif t == 'TETRA' or t == 'TRI' or t == 'NODE' or t == 'BAR':
            return array
        elif t == 'NGON':
            try: import Transform as T
            except: raise ImportError("convertArray2Tetra: requires Transform for NGONs.")
            tmp = T.breakElements(array)
            brd = []
            for i in tmp:
                if i[3] != 'NGON': brd.append(convertArray2Tetra1__(i))
                else:
                    brd.append(convertArray2Tetra1__(i, split="withBarycenters"))
            brd = T.join(brd)
            return brd
        elif t == 'HEXA*' or t == 'QUAD*':
            tmp = center2Node(array)
            tmp = converter.convertHexa2Tetra(tmp)
            return node2Center(tmp)
        elif t == 'PENTA*':
            tmp = center2Node(array)
            tmp = converter.convertPenta2Tetra(tmp)
            return node2Center(tmp)
        elif t == 'PYRA*':
            tmp = center2Node(array)
            tmp = converter.convertPyra2Tetra(tmp)
            return node2Center(tmp)
        elif t == 'TETRA*' or t == 'TRI*' or t == 'NODE*' or t == 'BAR*':
            return array
        else:
            raise TypeError("convertArray2Tetra: type of element is not taken into account: %s."%t)

    elif split == 'withBarycenters': # new points added at centers of elements and faces
        if t == 'STRUCT':
            if arrayC == []:
                return converter.convertStruct2TetraBary(array)
            else:
                return converter.convertStruct2TetraBaryBoth(array,arrayC)

        elif t == 'HEXA' or t == 'QUAD' or t == 'PENTA' or t == 'PYRA' or t == 'BAR':
            if arrayC == []:
                return converter.convertArray2TetraBary(array)
            else:
                return converter.convertArray2TetraBaryBoth(array, arrayC)

        elif t == 'TETRA' or t == 'TRI' or t == 'NODE':
            if arrayC == []: return array
            else: return [array,arrayC]

        elif t == 'NGON':
            if arrayC == []: return converter.convertNGon2TetraBary(array)
            else: 
                return converter.convertNGon2TetraBaryBoth(array, arrayC)
        elif t == 'HEXA*' or t == 'PENTA*' or t == 'PYRA*' or t == 'BAR*':
            tmp = center2Node(array)
            tmp =  converter.convertArray2TetraBary(tmp)
            return node2Center(tmp)
        elif t == 'TETRA*' or t == 'TRI*' or t == 'NODE*':
            return array
        else:
            raise TypeError("convertArray2Tetra: with split='withBarycenters', type of element is not taken into account: %s."%t)
    else:
        raise TypeError("convertArray2Tetra: type of split unknown: %s."%split)

# -- Convert array(s) to tetra
def convertArray2Tetra(array, split='simple'):
    """Convert a array in an unstructured tetra array.
    Unstructured array is triangular in 2D and tetrahedral in 3D.
    Usage: convertArray2Tetra(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array: 
            b.append(convertArray2Tetra1__(i, split=split))       
        return b
    else: return convertArray2Tetra1__(array, split=split)

# -- interne --
def convertArray2Hexa1__(array):
    try: sub = array[3]
    except: raise TypeError("convertArray2Hexa: arg must be an array.")
    if isinstance(sub, str): t = sub
    else: t = 'STRUCT'
    if t == 'STRUCT': return converter.convertStruct2Hexa(array)
    elif t == 'HEXA' or t == 'QUAD' or t == 'BAR' or t == 'NODE':
        return array
    elif t == 'HEXA*' or t == 'QUAD*' or t == 'BAR*' or t == 'NODE*':
        return array
    elif t == 'TRI' or t == 'TETRA' or t == 'PENTA':
        return converter.convertUnstruct2Hexa(array)
    elif t == 'NGON':
        try: import Transform as T
        except: raise ImportError("convertArray2Hexa: requires Transform for NGONs.")
        tmp = T.breakElements(array)
        brd = []
        for i in tmp:
            if (i[3] != 'NGON'): brd.append(convertArray2Hexa1__(i))
        brd = T.join(brd)
        return brd
    elif t == 'TRI*' or t == 'TETRA*' or t == 'PENTA*':
        return converter.convertUnstruct2Hexa(array)
    else: raise TypeError("convertArray2Hexa: type of element is not taken into account: %s."%t)

# -- convert arrays(s) to hexa
def convertArray2Hexa(array):
    """Convert a array in an unstructured hexa array.
    Unstructured array can be quad in 2D and hexa in 3D.
    Usage: convertArray2Hexa(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(convertArray2Hexa1__(i))
        return b
    else: return convertArray2Hexa1__(array)

def convertArray2NGon1__(array):
    try: sub = array[3]
    except: raise TypeError("convertArray2NGon: arg must be an array.")
    if isinstance(sub, str): t = sub
    else: t = 'STRUCT'
    if t == 'STRUCT': return converter.convertStruct2NGon(array)
    elif t == 'NGON': return array
    else: return converter.convertUnstruct2NGon(array)

def convertArray2NGon(array):
    """Convert a array in a NGON array.
    Usage: convertArray2NGon( array ) """
    if isinstance(array[0], list):
        b = []
        for i in array: b.append(convertArray2NGon1__(i))
        return b
    else: return convertArray2NGon1__(array)

def node2Center(array, accurate=0):
    """Convert array defined on nodes to array defined on centers.
    Usage: node2Center(array,accurate=0)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.node2Center(i,accurate))
        return b
    else:
        return converter.node2Center(array,accurate)

def center2Node(array, cellNType=0):
    """Convert array defined on centers to array defined on nodes.
    Usage: center2Node(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.center2Node(i, cellNType))
        for nob in xrange(len(b)):
            if len(b[nob]) == 4:
                eltType = b[nob][3]
                b[nob][3] = eltType.split('*')[0]
        return b
    else:
        b = converter.center2Node(array, cellNType)
        if len(b) == 4:
            eltType = b[3]; b[3] = eltType.split('*')[0]
        return b

def node2ExtCenter(array):
    """Convert array defined on nodes to an array defined on extended centers.
    Usage: node2ExtCenter(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.node2ExtCenter(i))
        return b
    else:
        return converter.node2ExtCenter(array)

def extCenter2Node(array):
    """Convert array defined on extended centers to an array defined on nodes.
    Usage: extCenter2Node(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.extCenter2Node(i))
        return b
    else:
        return converter.extCenter2Node(array)

def center2ExtCenter(array):
    """Convert array defined for centers to an array defined for extended centers.
    Usage: center2ExtCenter(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.center2ExtCenter(i))
        return b
    else:
        return converter.center2ExtCenter(array)

# -- convert arrays(s) to node array(s).
def convertArray2Node(array):
    """Convert an array in an unstructured node array.
    Usage: convertArray2Node(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.convertArray2Node(i))
        return b
    else:
        return converter.convertArray2Node(array)

#==============================================================================
def addGhostCellsNGon(arrayN, arrayC=[],depth=2):
    """Add depth layers of ghost cells in an NGON array.
    Usage: addGhostCellsNGon(arrayN,arrayC,depth)"""
    if isinstance(arrayN[0], list):
        b = []
        nzones = len(arrayN); nzonesC = len(arrayC)
        if nzones == nzonesC:
            for noz in xrange(nzones):
                if arrayN[noz] != [] and arrayC[noz] == []:
                    res = converter.addGhostCellsNGonNodes(arrayN[noz],depth)
                    b.append(res)
                elif arrayN[noz] == [] and arrayC[noz] != []:
                    res = converter.addGhostCellsNGonCenters(arrayC[noz],depth)
                    b.append(res)
                else:
                    res = converter.addGhostCellsNGonBoth(arrayN[noz], arrayC[noz],depth)
                    b.append(res)
        else:
            for noz in xrange(nzones):
                b.append(converter.addGhostCellsNGon(arrayN[noz],depth))
        return b
    else:
        if arrayN != [] and arrayC==[]: return converter.addGhostCellsNGonNodes(arrayN,depth)
        elif arrayN == [] and arrayC != []: return converter.addGhostCellsNGonCenters(arrayC,depth)
        else: return converter.addGhostCellsNGonBoth(arrayN, arrayC, depth)

def rmGhostCellsNGon(arrayN, arrayC = [], depth=2):
    """Delete depth layers of cells at exterior borders of a NGON mesh.
    Usage: rmGhostCellsNGon(arrayN, arrayC, depth)"""
    if isinstance(arrayN[0], list):
        b = []
        nzones = len(arrayN); nzonesC = len(arrayC)
        if nzones == nzonesC:
            for noz in xrange(nzones):
                if arrayN[noz] != [] and arrayC[noz] == []:
                    res = converter.rmGhostCellsNGonNodes(arrayN[noz],depth)
                    b.append(res)
                elif arrayN[noz] == [] and arrayC[noz] != []:
                    res = converter.rmGhostCellsNGonCenters(arrayC[noz],depth)
                    b.append(res)
                else:
                    res = converter.rmGhostCellsNGonBoth(arrayN[noz], arrayC[noz],depth)
                    b.append(res)
        else:
            for noz in xrange(nzones):
                b.append(converter.rmGhostCellsNGon(arrayN[noz],depth))
        return b
    else:
        if arrayN != [] and arrayC==[]: return converter.rmGhostCellsNGonNodes(arrayN,depth)
        elif arrayN == [] and arrayC != []: return converter.rmGhostCellsNGonCenters(arrayC,depth)
        else: return converter.rmGhostCellsNGonBoth(arrayN, arrayC, depth)

#==============================================================================
# Retourne une liste des noms de variables presents dans l'arrays
#==============================================================================
def getVarNames(a):
    """Get variable names.
    Usage: getVarNames(a)"""
    if isinstance(a[0], list):
        allvars = []
        for i in a:
            v = i[0].split(",")
            allvars.append(v)
    else: allvars = a[0].split(",")
    return allvars

#==============================================================================
# Fonctions de preconditionement (hook)
# IN: function: le nom de la fonction qui va utiliser le hook
#==============================================================================
def createHook(a, function='None'):
    """Create a hook for a given function.
    Usage: hook = createHook(a, function)"""
    if function == 'None': return None
    elif function == 'faceCenters': # 0
        # Retourne un KDT pour les centres des faces
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(converter.registerFaces(i))
            return b
        else:
            return converter.registerFaces(a)
    elif function == 'nodes': # 2
        # Retourne un KDT pour les noeuds
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(converter.registerNodes(i))
            return b
        else:
            return converter.registerNodes(a)
    elif function == 'elementCenters': # 3
        # Retourne un KDT pour les centres des elements
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(converter.registerElements(i))
                #b.append(converter.registerElements(convertArray2NGon(i)))
            return b
        else:
            return converter.registerElements(a)
            #return converter.registerElements(convertArray2NGon(a))
    elif function == 'extractMesh': # 1
        # Retourne un ADT pour les elements
        try: import Post as P
        except: inl = a
        else:
            inl, modified = P.growOfEps__(a, 1.e-6, nlayers=2, planarity=False)
        return converter.registerCells(inl)

    elif function == 'adt': # 1 ADT pour les interpolations
        return converter.registerCells(a)

    else: raise ValueError("function %s is invalid."%function)

#===============================================================================
# Fonctions de preconditionement (hook)
# IN: function: le nom de la fonction qui va utiliser le hook
# IN: a: liste des zones utilisees pour faire le preconditionnement global
# IN: extended: 0: pas d indirection sur les zones
#               1: sauvegarde de indirZones:
#                   indirection sur le no de la zone pour chq pt du hook
# OUT: hook(,indirZones)
#===============================================================================
def createGlobalHook(a, function='None', extended=0):
    """Create a hook for a given function.
    Usage: hook = createGlobalHook(a, function)"""
    if function == 'None': return None
    elif function == 'faceCenters': # 0
        # Retourne un KDT pour les centres des faces
        if not isinstance(a[0],list): return converter.registerAllFaces([a],extended)
        else: return converter.registerAllFaces(a,extended)
    elif function == 'nodes': # 2
        # Retourne un KDT pour les noeuds
        if not isinstance(a[0],list): return converter.registerAllNodes([a],extended)
        else: return converter.registerAllNodes(a,extended)
    elif function == 'elementCenters': # 3
        # Retourne un KDT pour les centres des elements
        if not isinstance(a[0],list): return converter.registerAllElements(convertArray2NGon([a]),extended)
        else: return converter.registerAllElements(convertArray2NGon(a),extended)
    elif function == 'extractMesh': # 1
        print 'function=extractMesh not implemented for global hook.'
    else: raise ValueError("function is invalid.")

#==============================================================================
def freeHook(hook):
    """Free hook.
    Usage: freeHook(hook)"""
    if isinstance(hook, list):
        for i in hook: converter.freeHook(i)
    else: converter.freeHook(hook)

#==============================================================================
# Fonctions d'identification geometrique
#==============================================================================
def identifyNodes(hook, a, tol=1.e-11):
    """Identify nodes of a in KDT. return identified indices.
    Usage: identifyNodes(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.identifyNodes(hook, i, tol))
        return b
    else:
        return converter.identifyNodes(hook, a, tol)

def identifyFaces(hook, a, tol=1.e-11):
    """Identify face centers of a in KDT. return identified indices.
    Usage: identifyFaces(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.identifyFaces(hook, convertArray2NGon(i), tol))
        return b
    else:
        return converter.identifyFaces(hook, convertArray2NGon(a), tol)

def identifyElements(hook, a, tol=1.e-11):
    """Identify element centers of a in KDT. return identified indices.
    Usage: identifyElement(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            if len(i) == 5: i = convertArray2Hexa(i)
            b.append(converter.identifyElements(hook, i, tol))
        return b
    else:
        if len(a) == 5: a = convertArray2Hexa(a)
        return converter.identifyElements(hook, a, tol)

#=============================================================================
def identifySolutions(coordsRcv, solDnr, hookDnr, vars = [], tol=1.e-13):
    """Identify points in a hook to mesh points and set the solution if donor
    and receptor points are distant from tol.
    Usage: identifySolutions(coordsRcv, solDnr, hookDnr, vars, tol)"""
    if vars != []: solDnr = extractVars(solDnr,vars)

    if isinstance(coordsRcv[0], list): # receptor is a list of zones
        if isinstance(solDnr[0], list):
            res = converter.identifySolutions(hookDnr, solDnr, coordsRcv, tol)
        elif not isinstance(solDnr[0], list): # une seule zone
            res = converter.identifySolutions(hookDnr, [solDnr], coordsRcv, tol)

    else: # receptor is a single zone
        if isinstance(solDnr[0], list):
            res = converter.identifySolutions(hookDnr, solDnr, [coordsRcv], tol)[0]
        elif not isinstance(solDnr[0], list): # une seule zone
            res = converter.identifySolutions(hookDnr, [solDnr], [coordsRcv], tol)[0]
    return res

#==============================================================================
# Fonctions d'identification du noeud/element/face le plus proche
#==============================================================================
def nearestNodes(hook, a):
    """Find in KDT nearest points to nodes of a. return identified indices.
    Usage: nearestNodes(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.nearestNodes(hook, i))
        return b
    else:
        return converter.nearestNodes(hook, a)

def nearestFaces(hook, a):
    """Find in KDT nearest points to face centers of a. return identified indices.
    Usage: nearestFaces(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.nearestFaces(hook, convertArray2NGon(i)))
        return b
    else:
        return converter.nearestFaces(hook, convertArray2NGon(a))

def nearestElements(hook, a):
    """Find in KDT nearest points to element centers of a. return identified indices.
    Usage: nearestElement(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            if len(i) == 5: i = convertArray2Hexa(i)
            b.append(converter.nearestElements(hook, i))
        return b
    else:
        if len(a)==5: a = convertArray2Hexa(a)
        return converter.nearestElements(hook, a)

# Retourne -1: la variable n'est presente dans aucun array
# Retourne 0: la variable est presente dans au moins un array
# Retourne 1: la variable est presente dans tous les arrays
def isNamePresent(a, varname):
    """Test if varName is present in a."""
    if isinstance(a[0], list):
        one = 0
        for i in a:
            p = KCore.isNamePresent(i, varname)
            if p != -1: one += 1
        if one == len(a): return 1
        elif one == 0: return -1
        else: return 0
    else: # un seul array
        p = KCore.isNamePresent(a, varname)
        if p == -1: return -1
        else: return 1

#==============================================================================
# Client/Server - send
#==============================================================================
def send(data, host='localhost', rank=0, port=15555):
    """Send data to socket."""
    import socket; import Compressor
    port = port+rank; sizeBuf = 1024

    #print 'connecting to port',port
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
    except:
        print 'Send: can not connect to %s [%d]. Nothing sent.'%(host,port)
        return

    data = Compressor.pack(data) # serialize
    size = len(data)

    # Blocks
    header = (size,sizeBuf)
    #print 'sending', header
    header = Compressor.pack(header)
    header = header.ljust(255)
    
    s.send(header)
    nbytes = 0
    while nbytes < size:
        if nbytes+sizeBuf > size:
            nbytes += s.send(data[nbytes:nbytes+sizeBuf])
        else: 
            nbytes += s.send(data[nbytes:])
        #print 'send',nbytes,size
    s.close()

#==============================================================================
# Client/server - createSockets
#==============================================================================
def createSockets(nprocs=1, port=15555):
    """Create sockets for communication."""
    import socket
    sockets = []
    for i in xrange(nprocs):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(('', port+i))
        s.listen(5)
        sockets.append(s)
    return sockets

#==============================================================================
# Client/server - listen
#==============================================================================
def listen(s):
    """Listen for sends."""
    import socket
    import Compressor
    while True:
        #s.listen(5)
        s.setblocking(1)
        client, address = s.accept()
        data = None
        nb = client.recv(255)
        if nb != "":
            nb = nb.rstrip()
            (size,sizeBuf) = Compressor.unpack(nb)
            #print 'Received ',size
            data = ''
            nbytes = 0
            while nbytes < size:
                if nbytes+sizeBuf < size:
                    received = client.recv(sizeBuf)
                else:
                    received = client.recv(size-nbytes)
                data += received; nbytes += len(received)
            data = Compressor.unpack(data)
            client.close()
            return data

#=============================================================================
# checkfile type
#=============================================================================
def checkFileType(fileName):
  """Find file type."""
  import binascii as b
  import os
  eol = "0a"
  inria = "MeshVersionUnformatted"
  try: file = open(fileName, 'r')
  except: raise IOError("checkFileType: file %s not found."%fileName)

  fileSize = os.path.getsize(fileName)
  header = file.read(512)  #lecture des 512 premiers octets
  file.seek(0)
  beader = b.b2a_hex(header)
  eolx1 = beader.find(eol, 0)
  eolx2 = beader.find(eol, eolx1 + 1)
  eol1 = (eolx1 +1)/2
  eol2 = (eolx2 +1)/2
  file.seek(0)
  dt = numpy.dtype('<i4')
  ieader = numpy.fromfile(fileName,dtype=dt,count=128,sep="")
  # Overflow sur ce clacul de sizet (uniquement util pour bin_stl)
  #ntri = numpy.frombuffer(header[80:],dtype='i4',count=1)  # format bin_stl, nombre de triangles
  #sizet = ntri[0]*50+84  #format bin_stl 80 octets d entete/nombre de triangles/50 octets par triangles
  sizet = 0

  if header[1:4] == 'HDF':
    fileType = 'bin_hdf'
  elif header[4:7] == 'ADF':
    fileType = 'bin_adf'
  elif header[0:5] == '#!TDV':
    fileType = 'bin_tp'
  elif header[0:5] == 'TITLE' or header[0:5] == 'title' or header[0:9] == "VARIABLES" or header[0:9] == "variables" or header[0:8] == "FILETYPE" or header[0:8] == "filetype":
    fileType = 'fmt_tp'
  elif ((beader[0:8] == "04000000" or beader[0:8] == "00000004") and (header[16:18] == 'va' or  header[16:17] == 'x' or header[16:17] == 'y' or header[16:17] == 'z' or header[16:18] == 'VA' or header[16:17] == 'X' or header[16:17] == 'Y' or header[16:17] == 'Z')): 
    fileType = 'bin_v3d' 
  elif ((beader[0:8] == "08000000" or beader[0:8] == "00000008") and (header[20:22] == 'va' or  header[20:21] == 'x' or header[20:21] == 'y' or header[20:21] == 'z' or header[20:22] == 'VA' or header[20:21] == 'X' or header[20:21] == 'Y' or header[20:21] == 'Z')): 
    fileType = 'bin_v3d' 
  elif (beader[10:12] == eol and (beader[50:52] == "78" or beader[50:52] == "79" or beader[50:52] == "80" or  beader[12:14] == "78" or   beader[12:14] == "79" or beader[12:14] == "80"  or beader[50:52] == "58" or beader[50:52] == "59" or beader[50:52] == "60" or  beader[12:14] == "58" or   beader[12:14] == "59" or beader[12:14] == "60" or beader[46:52] == "766172" or beader[46:52] == "564152" or beader[12:18] == "766172" or beader[12:18] == "564152")):
    fileType = 'fmt_v3d' 
  elif (header.find("MeshVersionUnformatted") != -1):
    fileType = 'bin_mesh'
  elif (header.find("MeshVersionFormatted") != -1):
    fileType = 'fmt_mesh'
  elif  (header.find("$MeshFormat") != -1):
    EndMesh=header.find("$EndMeshFormat")
    if EndMesh == 20:
      fileType = 'fmt_gmsh'
    elif EndMesh == 25:
      fileType = 'bin_gmsh'
  elif header.find("NDIME=") != -1:
    fileType = 'fmt_su2'
  elif header.find("DONNEES GENERALES") != -1 or header.find("--------") != -1:
    fileType = 'fmt_cedre'
  elif header[0:8] == 'CEDRE_IO':
    fileType = 'bin_cedre'
  elif fileSize == sizet:
    fileType = 'bin_stl'
  elif header.find("solid") == 0:
    fileType = 'fmt_stl'
  elif header.find("v ") != -1:
    fileType = 'fmt_obj'
  elif beader.find("4d4d") == 0:
    fileType = 'bin_3ds'
  elif header.find("ply") != -1:
    fileType = 'bin_ply'
  elif header.find("#FIG") == 0:
    fileType = 'fmt_xfig'
  elif header.find("<svg") != -1:
    fileType = 'fmt_svg'
  elif header.find("PNG") == 1:
    fileType = 'bin_png'
  elif ieader[0] == 4:
      ninjnk = ieader[1] * 3 * 4   # 3 pour ni,nj,nk, 4 pour 4 octets cas 3D multibloc
      ninj   = ieader[1] * 2 * 4 # 3 pour ni,nj, 4 pour 4 octets cas 2D multibloc
      if ieader[3] == ninjnk or ieader[3] == ninj: 
        fileType = 'bin_plot3d'
  elif ieader[0] == 12 or ieader[0] == 8:  # cas 2D ou 3d monobloc
        fileType = 'bin_plot3d'
  elif eol1 != 0:
      i = 0
      ligne0=[]
      ligne1=[]
      ninjnk_size = 0
      for line in file:
        if i == 0:
          ligne0= line.split()
          npi = int(ligne0[0])
        else:
          newline = line.split()
          sfloat=newline[0].count(".")
          if sfloat == 0: 
             ligne1.extend(newline)
          else: break
        i += 1
      ninjnk_size = ninjnk_size + len(ligne1)
      if (ninjnk_size == 2 * npi or ninjnk_size == 3 * npi):
        fileType = 'fmt_plot3d'
      else: 
        fileType = 'unknown'
  file.close()
  return fileType 
