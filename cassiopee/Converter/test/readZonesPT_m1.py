# - readZones (pyTree) -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

# Cree le fichier test HDF
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (10,10,10))
    b = G.cart((12,0,0), (1,1,1), (10,10,10))
    t = C.newPyTree(['Base',a,b])
    C.convertPyTree2File(t, 'test.hdf')
Cmpi.barrier()

# Relit des zones par procs
t = Cmpi.convertFile2SkeletonTree('test.hdf')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')
t = Cmpi.readZones(t, 'test.hdf', rank=Cmpi.rank)
if Cmpi.rank == 0: test.testT(t, 1)

# Cree le fichier test ADF
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (10,10,10))
    b = G.cart((12,0,0), (1,1,1), (10,10,10))
    t = C.newPyTree(['Base', a,b])
    C.convertPyTree2File(t, 'test.adf')
Cmpi.barrier()

# Relit des zones par procs
t = Cmpi.convertFile2SkeletonTree('test.adf')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')
t = Cmpi.readZones(t, 'test.adf', rank=Cmpi.rank)
if Cmpi.rank == 0: test.testT(t, 2)
