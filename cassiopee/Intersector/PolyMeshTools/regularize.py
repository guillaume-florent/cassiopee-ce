import Converter.PyTree as C
import Intersector.PyTree as XOR
import Converter.Internal as I
import time
import os, sys

if len(sys.argv) is not 4 :
    print "ARG ERROR : phmeshfile step# verbose#"
    sys.exit()

FILE=sys.argv[1]
STEP=int(sys.argv[2])
VERBOSE=int(sys.argv[3])

def stats(a):
  nzones = 0 ; npTot = 0 ; ncellsTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
      dim = I.getZoneDim(z)
      np = dim[1]
      ncells = dim[2]
      npTot += np
      ncellsTot += ncells
      nzones += 1
  print 'Nombre de zones : %d\n'%nzones
  print 'Nombre de points : %d\n'%npTot
  print 'Nombre de cellules : %d\n'%ncellsTot

def nb_cells(a):
  ncellsTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
      dim = I.getZoneDim(z)
      np = dim[1]
      ncells = dim[2]
      ncellsTot += ncells
  return ncellsTot

def nb_faces(a):
  ncfacesTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
  	GEl = I.getElementNodes(z)
  	NGON = 0; found = False
  	for c in GEl:
  		if c[1][0] == 22: found = True; break
  		NGON += 1
  	if found:
  		node = GEl[NGON]
  		er = I.getNodeFromName1(node, 'ElementRange')
        ncfacesTot += er[1][1]-er[1][0]+1

  return ncfacesTot

def aglomerate(t, vr, vm):
	nb_cells0 = nb_cells(t)
	carry_on=1
	i=0
	while (carry_on == 1):  
		print "iter %s"%i
		t=XOR.agglomerateSmallCells(t, vmin=vm, vratio=vr)
		nb_cells1 = nb_cells(t)
		if (nb_cells1 == nb_cells0): carry_on=0
		nb_cells0 = nb_cells1
		i=i+1

	return t

def collapse_bad_pgs(t):
	IMAX=20
	nb_faces0 = nb_faces(t)
	carry_on=1
	i=0
	while (carry_on == 1):  
		print "iter %s"%i
		
		q=XOR.collapseUncomputableFaces(t)

		# check for open cells
		err=XOR.checkCellsClosure(q)
		if (err==1):
                    print "ERROR : collasping has generated open cells"
                    return (1,t)
                # check for self intersections
                q=XOR.selfX(q)
                nbe = nb_cells(q)
                if (nbe > 0):
                    err=1
                    print "ERROR : collapsing has generated new Xs"
                    return (1,t)
                
                t=q
                
		nb_faces1 = nb_faces(t)
		if (nb_faces1 == nb_faces0): carry_on=0
		if (i == IMAX): carry_on=0
		if (err==1): carry_on=0
		nb_faces0 = nb_faces1
		i=i+1

	return (0,t)

t0 = time.time()
t = C.convertFile2PyTree(FILE)

err=0

if (VERBOSE >= 2):
	t2 = time.time()
	q=XOR.extractUncomputables(t)
	print ' - extractUncomputables CPU time : ',time.time()-t2,'s'
	C.convertPyTree2File(q, "uncomputable_0.cgns")

	t2 = time.time()
	q=XOR.extractNonStarCells(t)
	print ' - extractNonStarCells CPU time : ',time.time()-t2,'s'
	C.convertPyTree2File(q, "nonstar_0.cgns")

t1 = time.time()

if (STEP == 0):
	print "FIRST ROUND OF CLEANING : UNCOMPUTABLES"
	t2 = time.time()
	print "agglomerate cells with a ratio of 1000"
	t=aglomerate(t, vr=1000, vm=0.)
	print "agglomerate cells with a ratio of 100"
	t=aglomerate(t, vr=100, vm=0.)
	print "agglomerate cells with a ratio of 10"
	t=aglomerate(t, vr=10, vm=0.)
	print ' - agglomerating CPU time : ',time.time()-t2,'s'
	t2 = time.time()
	print " simplify Cells"
	t=XOR.simplifyCells(t)
	print ' - simplifyCells CPU time : ',time.time()-t2,'s'

	if (VERBOSE >= 1):
		C.convertPyTree2File(t, "cleaning1_result.cgns")

	print ' - FIRST ROUND CPU time : ',time.time()-t1,'s'

t1 = time.time()
if (STEP <= 1):
	print "SECOND ROUND OF CLEANING : COLLAPSING (remaining) UNCOMPUTABLES"

	(err,t)=collapse_bad_pgs(t)

	if (VERBOSE >= 1):
		C.convertPyTree2File(t, "cleaning2_result.cgns")

	print ' - SECOND ROUND CPU time : ',time.time()-t1,'s'

t1 = time.time()
if (STEP <= 2 and err == 0):
	print "THIRD ROUND OF CLEANING : SPLITTING NON-STAR CELLS"
	t2 = time.time()
	set = 1 # 0 for concave cells or 1 for non-centroid-star_shaped cells
	policy = 2 #0 : convexify concave pgs on PH set. 1 : starify concave pgs on PH set. 2 : starify any pgs at concave-chains ends
	print "prepareCellsSplit"
	t = XOR.prepareCellsSplit(t, PH_set = set, split_policy = policy, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8)
	print ' - prepareCellsSplit CPU time : ',time.time()-t2,'s'
	if (VERBOSE >= 3):
		C.convertPyTree2File(t, "prepareCellsSplit.cgns")
	t2 = time.time()
	print "splitNonStarCells"
	t = XOR.splitNonStarCells(t)
	print ' - splitNonStarCells CPU time : ',time.time()-t2,'s'
	if (VERBOSE >= 3):
		C.convertPyTree2File(t, "splitNonStarCells.cgns")
	t2 = time.time()
	print "simplifyCells"
	t = XOR.simplifyCells(t, 0)# do not treat externals
	print ' - simplifyCells CPU time : ',time.time()-t2,'s'
	if (VERBOSE >= 3):
		C.convertPyTree2File(t, "simplifyCells.cgns")

	if (VERBOSE >= 1):
		C.convertPyTree2File(t, "cleaning3_result.cgns")

	if (VERBOSE >= 2):
		t2 = time.time()
		print "extractNonStarCells"
		q=XOR.extractNonStarCells(t)
		C.convertPyTree2File(q, "nonstar_1.cgns")
		print ' - extractNonStarCells CPU time : ',time.time()-t2,'s'

t1 = time.time()
if (STEP <= 3 and err == 0):
	print "FOURTH ROUND OF CLEANING : COLLAPSING (generated) UNCOMPUTABLES"

	t=collapse_bad_pgs(t)

	if (VERBOSE >= 1):
		C.convertPyTree2File(t, "cleaning4_result.cgns")
	print ' - FOURTH ROUND CPU time : ',time.time()-t1,'s'

t1 = time.time()
if (STEP <= 4 and err == 0):
	print "FIFTH ROUND OF CLEANING : AGGLOMERATING SPLIT CELLS"

	t2 = time.time()
	print "agglomerate cells with a ratio of 1000"
	t=aglomerate(t, vr=1000, vm=0.)
	print "agglomerate cells with a ratio of 100"
	t=aglomerate(t, vr=100, vm=0.)
	print "agglomerate cells with a ratio of 10"
	t=aglomerate(t, vr=10, vm=0.)
	print ' - agglomerating CPU time : ',time.time()-t2,'s'

	t2 = time.time()
	print "simplifyCells"
	t=XOR.simplifyCells(t, 0)# do not treat externals
	print ' - simplifyCells CPU time : ',time.time()-t2,'s'

	if (VERBOSE >= 1):
		C.convertPyTree2File(t, "cleaning5_result.cgns")
	if (VERBOSE >= 2):
		t2 = time.time()
		print "extractNonStarCells"
		q=XOR.extractNonStarCells(t)
		print ' - extractNonStarCells CPU time : ',time.time()-t2,'s'
		C.convertPyTree2File(q, "nonstar_2.cgns")
	print ' - FIFTH ROUND CPU time : ',time.time()-t1,'s'

print ' TOTAL CPU TIME : ', time.time()-t0, 's'

if (err == 0): C.convertPyTree2File(t, "result.cgns")
else : print "ERROR due to collapsing"
