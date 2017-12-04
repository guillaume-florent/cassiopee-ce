# - convertArrays2File -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test
LOCAL = test.getLocal()

# Create test meshes
cart1 = G.cart((0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
cart2 = G.cartTetra((0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
cart3 = G.cartTetra((0,0,0), (0.1, 0.2, 1.), (11, 11, 1))
cart4 = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
l1 = D.line((0,0,50), (1,1,50))
l2 = D.line((0,0,1), (1,1,1))

# bin_tp
C.convertArrays2File([cart1, cart2], LOCAL+'/out.plt', 'bin_tp')
test.testF('out.plt', 1)

# fmt_tp
C.convertArrays2File([cart1], LOCAL+'/out.tp', 'fmt_tp', dataFormat='%.9e ')
test.testF('out.tp', 2)

# fmt_tp
C.convertArrays2File([cart1, cart2], LOCAL+'/out2.tp', 'fmt_tp', dataFormat='%.9e ',
                     zoneNames=['ZoneA', 'ZoneB'])
test.testF('out2.tp', 21)

# fmt_tp
C.convertArrays2File([cart4], LOCAL+'/out3.tp', 'fmt_tp',
                     dataFormat='%.9e ')
test.testF('out3.tp', 22)

# bin_v3d
C.convertArrays2File([cart1], LOCAL+'/outbe.v3d', 'bin_v3d')
test.testF('outbe.v3d', 3)

# bin_v3d little endian
C.convertArrays2File([cart1], LOCAL+'/outle.v3d', 'bin_v3d', endian='little')
test.testF('outle.v3d', 31)
    
# fmt_v3d
C.convertArrays2File([cart1], LOCAL+'/out.dat', 'fmt_v3d', dataFormat='%16.9e')
test.testF('out.dat', 4)
    
# bin_plot3d big endian
C.convertArrays2File([cart1], LOCAL+'/outbe.dat', 'bin_plot3d')
test.testF('outbe.dat.gbin', 5)

# bin_plot3d little endian
C.convertArrays2File([cart1], LOCAL+'/outle.dat', 'bin_plot3d', endian='little')
test.testF('outle.dat.gbin', 51)

# fmt_pov
C.convertArrays2File([cart3], LOCAL+'/out.pov', 'fmt_pov', dataFormat='%f ')
test.testF('out.pov', 6)
    
# fmt_mesh
C.convertArrays2File([cart3], LOCAL+'/out.mesh', 'fmt_mesh', dataFormat='%.9e ')
test.testF('out.mesh', 7)

# fmt_mesh
C.convertArrays2File([cart3], LOCAL+'/out.su2', 'fmt_su2', dataFormat='%.9e ')
test.testF('out.su2', 71)

# bin_stl
C.convertArrays2File([cart3], LOCAL+'/out.stl', 'bin_stl')
test.testF('out.stl', 75)

# fmt_obj
C.convertArrays2File([cart3], LOCAL+'/out.obj', 'fmt_obj', dataFormat='%f ')
test.testF('out.obj', 77)

# bin_pickle
C.convertArrays2File([cart3], LOCAL+'/out.pickle', 'bin_pickle')
test.testF('out.pickle', 8)

# fmt_xfig
C.convertArrays2File([l1], LOCAL+'/out.fig', 'fmt_xfig')
test.testF('out.fig', 9)

# fmt_svg
C.convertArrays2File([l2], LOCAL+'/out.svg', 'fmt_svg', dataFormat='%.9f ')
test.testF('out.svg', 10)

# fmt_cedre
C.convertArrays2File([cart4], LOCAL+"/out.d", "fmt_cedre")
test.testF('out.d', 11)
