import sys
import glob
import p_center_v1 as pc1
import p_center_v2 as pc2

myfiles = list(glob.glob('*.csv')) # list of *.csv files in the directory
if len(myfiles) == 0:
    print('No *.csv files found')
    sys.exit(1)

myfiles.sort()
TIEMPO = 3600  # Maximum execution time
CORTES = -1    # 0: not applicable, -1: default
HEUR = 0.05     # 0: not applicable, 0.05: default
PRE = -1       # 0: not applicable, -1: default
N = 2          # Problem dimension
P = [1, 2, 3, 4]  # Number of p-centers

if len(myfiles) == 0:
    print('No *.csv files found')
else:
    for p in P:
        for i in myfiles:
            fichdat = i[0:i.find('_')] # File name from the first position to the '_'
            fichinArcs = glob.glob(fichdat+'_*_arcs.txt')
            if len(fichinArcs) == 0:
                print('Cannot find any file '+fichdat+'*_arcs.txt')
                sys.exit(1)

            fichinArcs.sort()

            for k in fichinArcs:
                pc1.optimiza(i, k, N, p, TIEMPO, CORTES, HEUR, PRE)
                pc2.optimiza(i, k, N, p, TIEMPO, CORTES, HEUR, PRE)
