# - sphere6 (array) -
import Geom as D
import Converter as C
import KCore.test as test

a = D.sphere6((0,0,0), 1., 20)
test.testA(a, 1)
test.writeCoverage(100)
