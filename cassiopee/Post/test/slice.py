# - slice (array) -
import Generator as G
import Converter as C
import Post as P

a = G.cart((0,0,0), (1,1,1), (30,30,30))
a = C.initVars(a, '{f} = 3*{x}')
p = P.slice(a, type='cone', eq='2*{x}+{r}+1.=0')
print p
C.convertArrays2File(p, 'out.plt')
