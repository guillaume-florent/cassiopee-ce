# - naca (array) -
import Geom as D
import Converter as C

a = D.naca(12.)
C.convertArrays2File([a], 'out.plt', 'bin_tp')
