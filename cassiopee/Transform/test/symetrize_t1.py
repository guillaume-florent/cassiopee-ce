# - symetrize (array) -
import Generator as G
import Transform as T
import Converter as C
import KCore.test as test

test.stdTestA(T.symetrize, (0.,0.,0.), (1,0,0), (0,0,1))
