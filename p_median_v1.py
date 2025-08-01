def optimiza(datosrec, ficheroArcs, n, p, tiempo, cortes, heur, pre):

    # ========================================================================
    #
    #       P-MEDIAN
    #       DIMENSION 2
    #       EUCLIDEAN DISTANCE
    #       Neighborhoods with three areas (disks or rectangles).
    #
    #       Formulation: F2pMP
    #
    # ========================================================================

    from gurobipy import Model, GRB, quicksum, tuplelist
    import numpy as np
    import pandas as pd  # to read data from csv files
    import os

    modelo = 'P-median_v1'
    # abrev = 'PM1'
    print()
    print('============================================================================')
    print(modelo)

    # Name of the file with regions (Ejx)
    aux = pd.read_csv(datosrec)
    data = np.array(aux)
    print('Regions file: ', datosrec)

    # Name of the file with arcs (Ejy)
    print('Arcs file: ', ficheroArcs)

    Ej_dat = datosrec[datosrec.find('_')+1:datosrec.find('.',datosrec.find('_')+1)]
    Ej_arc = ficheroArcs[ficheroArcs.find('_')+1:ficheroArcs.find('_',ficheroArcs.find('_')+1)]

    RP = data[:, 0].astype(int)  # Integer vector with number of polygons in each region
    RD = data[:, 1].astype(int)  # Integer vector with number of disks in each region
    rn = (data[0, 0] + data[0, 1]).astype(int)  # number of figures in each neighborhood

    bigM = data[-1][0]  # The bigM comes in the last row of the data file

    data = np.delete(data, [0, 1], 1)

    data = np.delete(data, -1, 0)  # Delete the last row

    print('MAX_dist_aprox: ', bigM)

    N = range(n) # Dimension of the neighborhoods. Normally 2 (2D)
    r = len(data)  # number of neighborhoods must be read from the number of rows in the file
    print('Total regions: ', r)
    R = range(r)  # [0,1,..r-1]

    # Read the file with the arcs
    arcs = tuplelist((np.array(pd.read_csv(ficheroArcs, header=None))).tolist())
    arc = len(arcs)
    print('Total arcs: ', arc)
    # Calculate the arc density
    dens = round(arc/(r*(r-1)), 2)
    print('Density: ', dens)

    # Number of p-medians
    print('P-medians: ', str(p))
    print()

    # ------------------------------------------------------------------------
    #              CREATE A NEW MODEL AND VARIABLES
    # ------------------------------------------------------------------------
    mo = Model("P-MEDIAN-1")

    # ------------- variables that define each region or neighborhood --------------
    # u_i continuous variables represent points of the polygon.
    # R_i P_j C_k Environment i, polygon j coordinate k
    # w_i continuous variables represent points of the disk.
    # R_i RD_j C_k Environment i, disk j, coordinate k
    u = mo.addVars(R, range(rn), N, name="u")  # Ri RPj 0 and Ri RPj 1
    w = mo.addVars(R, range(rn), N, name="w")  # Ri RDj 0 and Ri RDj 1

    # t_ki=1 if the point is in polygon i of region k.
    # s_ki=1 if the point is in disk i of region k.
    t = mo.addVars(R, range(rn), vtype=GRB.BINARY, name="t")
    s = mo.addVars(R, range(rn), vtype=GRB.BINARY, name="s")

    # ------------------------  matching variables  -----------------------
    # x_ij=1 if neighborhood i is assigned to median j
    x = mo.addVars([(i, j) for i in R for j in R if i != j], vtype=GRB.BINARY, name="x")

    # y_j=1 if neighborhood j is a median
    y = mo.addVars(R, vtype=GRB.BINARY, name="y")

    # f_ij=1 if region i is connected to region j
    f = mo.addVars([(i, k, l) for k, l in arcs for i in R], vtype=GRB.BINARY, name="f")

    # dis_ij distances between zi and zj
    dis = mo.addVars([(i, k, l) for k, l in arcs for i in R], lb=0.0, name="dis")

    # z_i point in region i
    z = mo.addVars(R, N, name="z")
    zkl = mo.addVars([(i, k, l, d) for k, l in arcs for i in R for d in N], name="zij")

    # In the distance constraints between regions with lp norms
    # we will call rz=dij+M(1-fikl)
    rz = mo.addVars([(i, k, l) for k, l in arcs for i in R], name="rz")

    # -----------------defining objective------------------------------------

    obj = dis.sum()
    mo.setObjective(obj, GRB.MINIMIZE)

    # ---------------------------------------------------------
    #                  CONSTRAINTS
    # ----------------------------------------------------------

    # -----------------zi belongs to Ri------------------------
    # R is the num. of regions. RP[i] is the num. of polygons in region i
    mo.addConstrs(((z[i, d] == quicksum(u[i, j, d] for j in range(RP[i])) + quicksum(w[i, j, d] - data[i, RP[i]*(2*n)+j*(1+n)+d]*(1 - s[i, j]) for j in range(RD[i]))) for i in R for d in N), name='R1')

    # -------------constraints for RP and RD----------
    mo.addConstrs((-u[i, j, d] <= -data[i, (2*n)*j+d]*t[i, j] for i in R if RP[i] != 0 for j in range(RP[i]) for d in N), name='R2_ll')
    mo.addConstrs((u[i, j, d] <= data[i, (2*n)*j+n+d]*t[i, j] for i in R if RP[i] != 0 for j in range(RP[i]) for d in N), name='R2_ur')  # In 3D would be R2_tur (diagonal opposite to ll)

    for i in R:
        if RD[i] != 0:
            for j in range(RD[i]):
                mo.addQConstr((quicksum((w[i, j, d] - data[i, RP[i]*(2*n)+j*(1+n)+d])*(w[i, j, d] - data[i, RP[i]*(2*n)+j*(1+n)+d]) for d in N) <= data[i, RP[i]*(2*n)+j*(1+n)+n]*s[i, j]*data[i, RP[i]*(2*n)+j*(1+n)+n]*s[i, j]), name='R3')

    # --------------------- variables t and s ---------------------
    mo.addConstrs(((1 == quicksum(t[i, j] for j in range(RP[i])) + quicksum(s[i, j] for j in range(RD[i]))) for i in R), name='R4')

    # ----------------- P-median constraints ------------------
    mo.addConstr(p == y.sum(), name='R22')
    mo.addConstrs((x[i, j] <= y[j] for i in R for j in R if i != j), name='R23')
    mo.addConstrs(((1 == (y[i]+(quicksum(x[i, j] for j in R if j != i)))) for i in R), name='R24')
    mo.addConstrs((-x[i, hk] == ((quicksum(f[i, k, l] for k, l in arcs if k != i and k == hk) - quicksum(f[i, l, k] for l, k in arcs if k != i and k == hk))) for i in R for hk in R if i != hk), name='R35')
    mo.addConstrs(((1-y[i] == (quicksum(f[i, k, l] for k, l in arcs if k == i))) for i in R), name='R36')
    mo.addConstrs(((1-y[kh] >= (quicksum(f[i, k, l] for k, l in arcs if k == kh))) for i in R for kh in R if i != kh), name='R37') # MAY BE REDUNDANT?

    # --- Constraints for SOC of euclidean distances between regions ---
    mo.addConstrs(((rz[i, k, l] == dis[i, k, l] + bigM*(1-f[i, k, l])) for i in R for k, l in arcs), name='R40_0') # Studies were done with ==, but the paper indicates <=
    mo.addConstrs((zkl[i, k, l, d] >= z[k, d] - z[l, d] for i in R for k, l in arcs for d in N), name='R40_1')
    mo.addConstrs((zkl[i, k, l, d] >= -z[k, d] + z[l, d] for i in R for k, l in arcs for d in N), name='R40_2')

    for i in R:
        for k, l in arcs:
            mo.addQConstr((quicksum(zkl[i, k, l, d]*zkl[i, k, l, d] for d in N) <= rz[i, k, l]*rz[i, k, l]), name='R40_3')

    # ***********************************************
    #               Solving
    # ***********************************************

    # -------parameters----------
    mo.Params.timelimit = tiempo
    mo.Params.Cuts = cortes
    mo.Params.Heuristics = heur
    mo.Params.Presolve = pre

    mo.update()  # Updates the model with variables and constraints

    # -------Solving----------
    mo.optimize()
    btime = mo.Runtime

    # ***************************************************
    #                  outputs
    # ***************************************************

    if not os.path.exists('Salidas/Resultados'):
        os.makedirs('Salidas/Resultados')
    fichero = 'Salidas/Resultados/resultados.txt'

    outfile = open(fichero, 'a')

    # 3=Infeasible, 4=infeasible or unbounded and 5=unbounded
    if mo.status == 3:
        print('Model not solved because it is infeasible')
        print()
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:3d}  {:7.1f}    {:5.1f}  {:5.1f}  Infeasible'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre), file=outfile)
        return

    if mo.status == 4:
        print('Model not solved because it is infeasible or unbounded')
        print()
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:3d}  {:7.1f}    {:5.1f}  {:5.1f}  Infeas.-Unbounded'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre), file=outfile)
        return

    if mo.status == 5:
        print('Model not solved because it is unbounded')
        print()
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:7.1f}  {:5.1f}  {:5.1f}  {:5.1f}  Unbounded'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre), file=outfile)
        return

    if mo.status == GRB.Status.OPTIMAL:
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:7.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:8.2f}  {:5.1f}   {:5.5f}  {:10.0f}'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre, mo.ObjVal, btime, mo.MIPGap, mo.NodeCount), file=outfile)
        
    else:
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:7.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:8.2f}  {:5.1f}   {:5.5f}  {:10.0f}'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre, mo.ObjVal, btime, mo.MIPGap, mo.NodeCount), file=outfile)
        
    outfile.close()
