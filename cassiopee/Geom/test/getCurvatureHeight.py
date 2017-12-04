# - getCurvatureHeight (array) -
import Converter as C
import Geom as D
import Transform as T

a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a = T.join (a1, a2)
hmax = D.getCurvatureHeight( a )
a = C.addVars([a,hmax])
C.convertArrays2File([a], 'out.plt')
