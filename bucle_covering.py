import sys
import glob
import numpy as np
import pandas as pd
import p_covering_v1 as pc1
import p_covering_v2 as pc2

myfiles = list(glob.glob('*.csv')) # list of *.csv files in the directory
if len(myfiles) == 0:
    print('No *.csv files found')
    sys.exit(1)

myradios = list(glob.glob('cov_*.txt')) # list of cov_*.txt files in the directory
if len(myradios) == 0:
    print('No cov_*.txt files found')
    sys.exit(1)

myfiles.sort()
TIEMPO = 3600  # Maximum execution time
CORTES = -1    # 0: not applicable, -1: default
HEUR = 0.05     # 0: not applicable, 0.05: default
PRE = -1       # 0: not applicable, -1: default
N = 2          # Problem dimension
P = [1, 2, 3, 4]  # Number of p-centers

for p in P:
    for i in myfiles:
        fichdat = i[0:i.find('_')] # File name from the first position to the '_'
        fichinArcs = glob.glob(fichdat+'_*_arcs.txt')
        if len(fichinArcs) == 0:
            print('Cannot find any file '+fichdat+'*_arcs.txt')
            sys.exit(1)

        fichinArcs.sort()

        for k in fichinArcs:
            da = k[k.find('_')+4:k.find('_',k.find('_')+1)]  # Arc density d=1 (0.07), d=2 (0.10), d=3 (0.15) and d=4 (0.20)
            for l in myradios:
                pc = l[l.find('_')+1:l.find('_',l.find('_')+1)]  # Percentage of demand we want to cover
                d = l[l.find('_',l.find('_')+1)+1:l.find('.',l.find('_'))]  # densities d=1 (0.07), d=2 (0.10), d=3 (0.15) and d=4 (0.20)

                if d == da:
                    # Read the file with the radii
                    radio = pd.read_csv('cov_'+pc+'_'+d+'.txt', header=None)

                    # Read the file with regions (Ejx) to see how many there are
                    aux = pd.read_csv(i)
                    data = np.array(aux)
                    r = len(data)  # number of neighborhoods must be read from the number of rows in the file
                    R = int(radio[p-1][r-2]) # Read the radius that corresponds according to the number of neighborhoods and p

                    pc1.optimiza(i, k, N, p, R, pc, TIEMPO, CORTES, HEUR, PRE)
                    pc2.optimiza(i, k, N, p, R, pc, TIEMPO, CORTES, HEUR, PRE)
