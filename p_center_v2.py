def optimiza(datosrec, ficheroArcs, n, p, tiempo, cortes, heur, pre):

    # ========================================================================
    #
    #       P-CENTER
    #       DIMENSION 2
    #       EUCLIDEAN DISTANCE
    #       Neighborhoods with three areas (disks or rectangles).
    #
    #       Formulation: PpCP
    #
    # ========================================================================

    from gurobipy import Model, GRB, quicksum, tuplelist
    import numpy as np
    import pandas as pd  # to read data from csv files
    import os

    modelo = 'P-center_v2'
    # abrev = 'PC2'
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

    data = np.delete(data, [0, 1], 1) # Delete the first row with the header
    data = np.delete(data, -1, 0)  # Delete the last row

    print('MAX_dist_aprox: ', bigM)

    N = range(n)  # dimension of the neighborhoods. Normally 2 (2D)
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

    # Number of p-centers
    print('P-centers: ', str(p))
    print()

    # ------------------------------------------------------------------------
    #              CREATE A NEW MODEL AND VARIABLES
    # ------------------------------------------------------------------------
    mo = Model("P-CENTER-2")

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
    # x_ij=1 if neighborhood i is assigned to center j
    x = mo.addVars([(i, j) for i in R for j in R if i != j], vtype=GRB.BINARY, name="x")

    # y_j=1 if neighborhood j is a center
    y = mo.addVars(R, vtype=GRB.BINARY, name="y")

    # h_ij=1 if neighborhood i is linked to neighborhood j
    h = mo.addVars([(k, l) for k, l in arcs], vtype=GRB.BINARY, name="h")

    # dd_i distance from i to its facility
    dd = mo.addVars([(i) for i in R], lb=0.0, name="dd")

    # dis_kl distances between zk and zl
    dis = mo.addVars([(k, l) for k, l in arcs], lb=0.0, name="dis")

    # z_i point in region i. zkl is an auxiliary variable to calculate distance between neighborhoods k and l.
    z = mo.addVars(R, N, name="z")
    zkl = mo.addVars([(k, l, d) for k, l in arcs for d in N], name="zij")

    # In the distance constraints between regions with lp norms we will call rz=dij+M(1-fikl)
    rz = mo.addVars([(k, l) for k, l in arcs], name="rz")

    # Zmin is to minimize the maximum distance
    Zmin = mo.addVar(name="Zmin")

    # -----------------defining objective------------------------------------

    mo.setObjective(Zmin, GRB.MINIMIZE)

    # ---------------------------------------------------------
    #                  CONSTRAINTS
    # ----------------------------------------------------------

    # ----------------- zi belongs to Ri ------------------------
    # R is the num. of regions. RP[i] is the num. of polygons in region i
    mo.addConstrs(((z[i, d] == quicksum(u[i, j, d] for j in range(RP[i])) + quicksum(w[i, j, d] - data[i, RP[i]*(2*n)+j*(1+n)+d]*(1 - s[i, j]) for j in range(RD[i]))) for i in R for d in N), name='R1')

    # ------------- constraints for RP and RD ----------
    mo.addConstrs((-u[i, j, d] <= -data[i, (2*n)*j+d]*t[i, j] for i in R if RP[i] != 0 for j in range(RP[i]) for d in N), name='R2_ll')
    mo.addConstrs((u[i, j, d] <= data[i, (2*n)*j+n+d]*t[i, j] for i in R if RP[i] != 0 for j in range(RP[i]) for d in N), name='R2_ur')  # En 3D seria la R2_tur (diagonal opuesta a ll)

    for i in R:
        if RD[i] != 0:
            for j in range(RD[i]):
                mo.addQConstr((quicksum((w[i, j, d] - data[i, RP[i]*(2*n)+j*(1+n)+d])*(w[i, j, d] - data[i, RP[i]*(2*n)+j*(1+n)+d]) for d in N) <= data[i, RP[i]*(2*n)+j*(1+n)+n]*s[i, j]*data[i, RP[i]*(2*n)+j*(1+n)+n]*s[i, j]), name='R3')

    # ---------------------- variables t and s ----------------------
    mo.addConstrs(((1 == quicksum(t[i, j] for j in range(RP[i])) + quicksum(s[i, j] for j in range(RD[i]))) for i in R), name='R4')

    # ----------------- P-center constraints ------------------
    mo.addConstr(p == y.sum(), name='R22')
    mo.addConstrs((x[i, j] <= y[j] for i in R for j in R if i != j), name='R23')
    mo.addConstrs(((Zmin >= dd[i]) for i in R), name='R58')
    mo.addConstrs(((dd[i] >= dd[j] + dis[i, j] - (r-p)*bigM*(1 - h[i, j])) for i, j in arcs), name='R59')
    mo.addConstrs((((quicksum(x[i, j] for j in R if j != i)) == (quicksum(h[k, l] for k, l in arcs if k == i))) for i in R), name='R60')
    mo.addConstrs(((1 == (y[i] + quicksum(h[i, l] for k, l in arcs if k == i))) for i in R), name='R61')
    mo.addConstrs(((h[k, i] <= ((quicksum(h[j, l] for j, l in arcs if l != k and j == i)) + y[i])) for k, i in arcs), name='R62')
    mo.addConstrs((x[i, k] >= h[i, j] + x[j, k] - 1 for i, j in arcs for k in R if i != k and j != k), name='R63')
    mo.addConstrs((x[i, j] >= h[i, j] + y[j] - 1 for i, j in arcs), name='R64')

    # --- Constraints for SOC of euclidean distances between regions ---
    mo.addConstrs(((rz[k, l] == dis[k, l] + bigM*(1-h[k, l])) for k, l in arcs), name='R59_0')
    mo.addConstrs((zkl[k, l, d] >= z[k, d] - z[l, d] for k, l in arcs for i in R for d in N), name='R59_1')
    mo.addConstrs((zkl[k, l, d] >= -z[k, d] + z[l, d] for k, l in arcs for i in R for d in N), name='R59_2')

    for k, l in arcs:
            mo.addQConstr((quicksum(zkl[k, l, d]*zkl[k, l, d] for d in N) <= rz[k, l]*rz[k, l]), name='R59_3')

    # ***********************************************
    #               Solving
    # ***********************************************

    # -------parameters----------
    mo.Params.timelimit = tiempo
    mo.Params.Cuts = cortes
    mo.Params.Heuristics = heur
    mo.Params.Presolve = pre

    mo.update()  # Updates the model with variables and constraints

    # -------- Solving --------
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
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:3d}  {:7.1f}    {:5.1f}  {:5.1f}  Infactible'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre), file=outfile)
        return

    if mo.status == 4:
        print('Model not solved because it is infeasible or unbounded')
        print()
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:3d}  {:7.1f}    {:5.1f}  {:5.1f}  Infact.-No_acotado'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre), file=outfile)
        return

    if mo.status == 5:
        print('Model not solved because it is unbounded')
        print()
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:7.1f}  {:5.1f}  {:5.1f}  {:5.1f}  No_acotado'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre), file=outfile)
        return

    if mo.status == GRB.Status.OPTIMAL:
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:7.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:8.2f}  {:5.1f}   {:5.5f}  {:10.0f}'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre, mo.ObjVal, btime, mo.MIPGap, mo.NodeCount), file=outfile)
        
    else:
        print('{:17} {:5} {:5} {:3d} {:2d}  {:5.2f}  {:7.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:8.2f}  {:5.1f}   {:5.5f}  {:10.0f}'.format(modelo, Ej_dat, Ej_arc, r, p, dens, tiempo, cortes, heur, pre, mo.ObjVal, btime, mo.MIPGap, mo.NodeCount), file=outfile)
        
    outfile.close()
