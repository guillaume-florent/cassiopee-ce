# - naca (array) -
import Geom as D
import KCore.test as test

a = D.naca(12.)
test.testA([a],1)

a = D.naca(12., N=501)
test.testA([a],2)

test.writeCoverage(100)
