# - getDistantIndex (array) -
import Geom as D

a = D.line((0.,0.,0.), (1.,0.,0), 100)
print 'distant Index:', D.getDistantIndex(a, 25, 0.2)
print 'distant Index:', D.getDistantIndex(a, 25, -0.2)
