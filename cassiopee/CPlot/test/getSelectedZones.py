# - getSelectedZones (array) -
import Generator as G
import CPlot
import time

a1 = G.cart( (0,0,0), (1,1,1), (5,5,5) )
a2 = G.cart( (7,0,0), (1,1,1), (3,3,3) )
CPlot.display([a1, a2])

ret = []
while (ret == []):
    ret = CPlot.getSelectedZones(); time.sleep(2.)
print 'Zones have been selected: ', ret
