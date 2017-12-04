# - addGhostCells with periodicity -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import math
import KCore.test as test
ni = 11; nj = 5; nk = 3
dk = 1./max(1,(nk-1))
a = G.cart((0.2,0,0),(1./(ni-1), 1./(nj-1),dk),(ni,nj,nk))
a = C.initVars(a,'VelocityX={CoordinateX}')
a = C.initVars(a,'VelocityY={CoordinateY}')
a = C.initVars(a,'VelocityZ=1.')
a = C.addBC2Zone(a,'per1','BCMatch',[1,1,1,nj,1,nk], zoneDonor=a,
                 rangeDonor=[ni,ni,1,nj,1,nk], trirac=[1,2,3],
                 rotationCenter=(0,0,0), rotationAngle=(0.,0.,0.),
                 translation=(+1.,0.,0.))
a = C.addBC2Zone(a,'per2','BCMatch',[ni,ni,1,nj,1,nk], zoneDonor=a,
                 rangeDonor=[1,1,1,nj,1,nk], trirac=[1,2,3],
                 rotationCenter=(0,0,0), rotationAngle=(0.,0.,0.),
                 translation=(-1.,0.,0.))

t = C.newPyTree(['Base',a])
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=1)
test.testT(t,1)

# cas periodicite par rotation
ni = 91; nj = 21; nk = 1
alpha = 90.
a = G.cylinder((0,0,0),0.5,1.,0.,alpha,1.,(ni,nj,nk))
a = C.initVars(a,'VelocityX=0.')
a = C.initVars(a,'VelocityY=0.')
a = C.initVars(a,'VelocityZ=1.')
vx = C.getField('VelocityX',a)[0]
vy = C.getField('VelocityY',a)[0]
nic = vx[2]; njc = vx[3]
for j in xrange(njc):
    for i in xrange(nic):
        i0 = i * math.pi/2./90.
        vx[1][0,i+j*nic] = math.cos(i0)
        vy[1][0,i+j*nic] = math.sin(i0)
C.setFields([vx],a,loc='nodes')
C.setFields([vy],a,loc='nodes')
a = C.addBC2Zone(a,'per1','BCMatch',[1,1,1,nj,1,nk], zoneDonor=a,
                 rangeDonor=[ni,ni,1,nj,1,nk], trirac=[1,2,3],
                 rotationCenter=(0,0,0), rotationAngle=(0.,0.,-alpha),
                 translation=(0.,0.,0.))
a = C.addBC2Zone(a,'per2','BCMatch',[ni,ni,1,nj,1,nk], zoneDonor=a,
                 rangeDonor=[1,1,1,nj,1,nk], trirac=[1,2,3],
                 rotationCenter=(0,0,0), rotationAngle=(0.,0.,alpha),
                 translation=(0.,0.,0.))
t = C.newPyTree(['Base',a])
t2 = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=1)
test.testT(t2,2)

a = G.cylinder((0,0,0),0.5,1.,0.,alpha,1.,(ni,nj,nk))
a = C.initVars(a,'VelocityX=0.')
a = C.initVars(a,'VelocityY=0.')
a = C.initVars(a,'VelocityZ=1.')
vx = C.getField('VelocityX',a)[0]
vy = C.getField('VelocityY',a)[0]
nic = vx[2]; njc = vx[3]
for j in xrange(njc):
    for i in xrange(nic):
        i0 = i * math.pi/2./90.
        vx[1][0,i+j*nic] = math.cos(i0)
        vy[1][0,i+j*nic] = math.sin(i0)
C.setFields([vx],a,loc='nodes')
C.setFields([vy],a,loc='nodes')
a = C.addBC2Zone(a,'per1','BCMatch',[1,1,1,nj,1,nk], zoneDonor=a,
                 rangeDonor=[ni,ni,1,nj,1,nk], trirac=[1,2,3],
                 rotationCenter=(0,0,0), rotationAngle=(0.,0.,-alpha),
                 translation=(0.,0.,0.))
a = C.addBC2Zone(a,'per2','BCMatch',[ni,ni,1,nj,1,nk], zoneDonor=a,
                 rangeDonor=[1,1,1,nj,1,nk], trirac=[1,2,3],
                 rotationCenter=(0,0,0), rotationAngle=(0.,0.,alpha),
                 translation=(0.,0.,0.))
t2 = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=0)
test.testT(t2,3)
