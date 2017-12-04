# - redispatch (pyTree) -
import Converter.PyTree as C
import Distributor2.PyTree as D2
import Distributor2.Mpi as D2mpi
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Connector.PyTree as X
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# Case
N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in xrange(N):
    a = G.cart((pos,0,0), (1,1,1), (10+i, 10, 10))
    pos += 10 + i - 1
    t[2][1][2].append(a)
t = X.connectMatch(t)
if Cmpi.rank == 0: C.convertPyTree2File(t, 'in.cgns')
Cmpi.barrier()

# lecture du squelette
a = Cmpi.convertFile2SkeletonTree('in.cgns')

# equilibrage 1
(a, dic) = D2.distribute(a, NProc=Cmpi.size, algorithm='fast', useCom=0)

# load des zones locales dans le squelette
a = Cmpi.readZones(a, 'in.cgns', rank=Cmpi.rank)

# equilibrage 2 (a partir d'un squelette charge)
(a, dic) = D2.distribute(a, NProc=Cmpi.size, algorithm='gradient1', 
                         useCom='match')

a = D2mpi.redispatch(a)
if Cmpi.rank == 0: test.testT(a, 1)

# force toutes les zones sur 0
zones = Internal.getNodesFromType(a, 'Zone_t')
for z in zones:
    nodes = Internal.getNodesFromName(z, 'proc')
    Internal.setValue(nodes[0], 0)

a = D2mpi.redispatch(a)
if Cmpi.rank == 0: test.testT(a, 2)
