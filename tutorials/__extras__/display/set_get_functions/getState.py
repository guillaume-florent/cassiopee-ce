# - getState (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (5,5,5) )
CPlot.display([a])

print 'dim=',CPlot.getState('dim')
print 'mode=',CPlot.getState('mode')
print 'displayBB=',CPlot.getState('displayBB')
print 'displayInfo=',CPlot.getState('displayInfo')
print 'meshStyle=',CPlot.getState('meshStyle')
print 'solidStyle=',CPlot.getState('solidStyle')
print 'isoStyle=',CPlot.getState('isoStyle')
print 'win=',CPlot.getState('win')
print 'posCam=',CPlot.getState('posCam')
print 'posEye=',CPlot.getState('posEye')
print 'dirCam=',CPlot.getState('dirCam')
