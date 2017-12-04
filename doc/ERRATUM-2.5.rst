KCore/Dist.py
-------------

(lignes 1008 et lignes 1091), ajouter une étoile à libgomp.so


if useOMP:
               l = checkLibFile__('libgomp.so*', additionalLibPaths)
               if l is None:
                    l = checkLibFile__('libgomp.a', additionalLibPaths)
               if l is not None:
                    libs += ['gomp']; paths += [l]
               else: ret = False


CPlot/apps/tkCassiopee.py
-------------------------

ligne 38 : Supprimer tkFastSolver.

SOLVERAPPS = ['tkInit', 'tkDistributor', 'tkDist2Walls', '---',
              'tkCassiopeeSolver', 'tkElsaSolver']
