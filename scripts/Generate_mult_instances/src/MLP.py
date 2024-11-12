# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import itertools as itt
import time

import numpy as np
import scipy.optimize as sp


class struct:
    pass


"""
MultiLinear Program class
"""


class MLP(object):
    """
    Instantiate the class
        nvar - number of variables
        var_lb - array of lower bounds on variables
        var_ub - array of upper bounds on variables
        var_type - array of variable type - 'C' or 'B'
        obj_sense - string defining objective sense - Min or Max
        obj_coeff - dict storing as: key -- the index of variables in monomial and value -- the coefficient in the objective
        obj_offset - constant term in the objective
        ncon - number of constraints
        cons_coeff - array of dict storing as: key -- the index of variables in monomial and value -- the coefficient in each constraint
        cons_ub - upper bound on constraints
        name - instance name
    """

    def __init__(
        self,
        mlp_dir=None,
        mlp_file=None,
        nvar=None,
        var_lb=None,
        var_ub=None,
        var_type=None,
        obj_sense="Min",
        obj_coeff=None,
        obj_offset=0.0,
        ncon=None,
        cons_coeff=None,
        cons_ub=None,
        name=None,
    ):

        if mlp_file is None:
            self.nvar = nvar
            self.var_lb = var_lb
            self.var_ub = var_ub
            self.var_type = var_type
            self.obj_sense = obj_sense
            self.obj_coeff = obj_coeff
            self.obj_offset = obj_offset
            self.ncon = ncon
            self.cons_coeff = cons_coeff
            self.cons_ub = cons_ub
            self.name = name

            self.obj_maxd = max([len(mon) for mon in self.obj_coeff.keys()])
            self.cons_maxd = [max([len(mon) for mon in self.cons_coeff[ic].keys()]) for ic in range(self.ncon)]

            self.time_readmlp = 0.0
        else:
            st = time.perf_counter()
            self.readMLP(mlp_dir, mlp_file)
            self.time_readmlp = time.perf_counter() - st

        st = time.perf_counter()
        self.finalizeMLP()
        self.time_finalizemlp = time.perf_counter() - st

    """
    Read the polynomial with monomials defined in the string array
        lines - string array in which each entry defines a monomial
        nmon - number of monomials in the polynomial
    """

    def readPolynomial(self, lines, nmon):
        maxd = 0
        coeff = {}
        for im in range(nmon):
            words = lines[im].split("\n")[0]
            iw = 0
            while not (words[iw] == "]"):
                iw = iw + 1
            val = float(words[iw + 1 :])
            inds = []
            inds_st = words[1:iw].split(",")
            inds = []
            for ii in range(len(inds_st)):
                inds.append(int(inds_st[ii]))
            coeff.update({tuple(np.sort(inds)): val})
            maxd = max(maxd, len(inds))
        return coeff, maxd

    """
    Reads the MLP instance from text file
        dir - folder name
        fname - file name
    """

    def readMLP(self, dir, fname):

        words = fname.split(".")
        name = words[0]

        fp = open(dir + fname, "r")

        lines = fp.readlines()

        cnt = 0
        # num variables
        words = lines[cnt].split("\n")[0].split()
        nvar = int(words[1])
        cnt = cnt + 1

        # num constraints
        words = lines[cnt].split("\n")[0].split()
        ncon = int(words[1])
        cnt = cnt + 1

        # objective sense
        words = lines[cnt].split("\n")[0].split()
        obj_sense = words[1]
        cnt = cnt + 1

        # variables info
        cnt = cnt + 1
        vlb = np.zeros((nvar, 1))
        vub = np.zeros((nvar, 1))
        vtype = []
        for iv in range(nvar):
            words = lines[cnt].split("\n")[0].split()
            vlb[iv] = float(words[0])
            vub[iv] = float(words[1])
            if words[2] == "Cont":
                vtype.append("C")
            if words[2] == "Bin":
                vtype.append("B")
            cnt = cnt + 1

        # read objective
        words = lines[cnt].split("\n")[0].split()
        nmon = int(words[1])
        cnt = cnt + 1
        words = lines[cnt].split("\n")[0].split()
        obj_offset = float(words[1])
        cnt = cnt + 1
        obj_coeff, obj_maxd = self.readPolynomial(lines[cnt : cnt + nmon], nmon)
        cnt = cnt + nmon

        # read constraints
        cons_coeff = []
        cons_ub = []
        cons_maxd = []
        for ic in range(ncon):
            words = lines[cnt].split("\n")[0].split()
            nmon = int(words[1])
            cnt = cnt + 1
            words = lines[cnt].split("\n")[0].split()
            ub = float(words[1])
            cnt = cnt + 1
            con_coeff, con_maxd = self.readPolynomial(lines[cnt : cnt + nmon], nmon)
            cons_coeff.append(con_coeff)
            cons_ub.append(ub)
            cons_maxd.append(con_maxd)
            cnt = cnt + nmon

        fp.close()

        self.nvar = nvar
        self.var_lb = vlb
        self.var_ub = vub
        self.var_type = vtype
        self.ncon = ncon
        self.obj_sense = obj_sense
        self.obj_coeff = obj_coeff
        self.obj_offset = obj_offset
        self.obj_maxd = obj_maxd
        self.cons_coeff = cons_coeff
        self.cons_ub = cons_ub
        self.cons_maxd = cons_maxd
        self.name = name

    """
    Computes other quantities needed for minimum and best bound linearization
    """

    def finalizeMLP(self):

        # find the unique monmials
        self.uniqueMonomials()

        # construct the unified linearization diagram
        self.all_nodes_uld = self.generateNodes()
        self.all_arcs_uld, self.all_triples_uld = self.generateArcsTriples(self.all_nodes_uld)

        if self.ncon > 0:
            self.obj_nodes_uld = self.generateNodes(self.obj_coeff.keys())
            self.obj_arcs_uld, self.obj_triples_uld = self.generateArcsTriples(self.obj_nodes_uld)

        self.all_obj_coeff = {mon: 0.0 for mon in self.all_nodes_uld}
        for mon in self.obj_coeff.keys():
            if self.obj_sense == "Min":
                self.all_obj_coeff[mon] = self.obj_coeff[mon]
            else:
                self.all_obj_coeff[mon] = -self.obj_coeff[mon]

        # construct the linearization diagram for each monomial
        self.monomial_nodes_ld = {}
        self.monomial_arcs_ld = {}
        self.monomial_triples_ld = {}
        for d in range(2, self.maxd + 1):
            for mon in self.unique_monomials[d]:
                self.monomial_nodes_ld.update({mon: self.generateNodes([mon])})
                arcs, triples = self.generateArcsTriples(self.monomial_nodes_ld[mon])
                self.monomial_arcs_ld.update({mon: arcs})
                self.monomial_triples_ld.update({mon: triples})

        # set lower and upper bounds for all variables
        self.findLowerUpperBounds()

        self.lam1bnd, self.lam2bnd, self.lam3bnd = self.computeUpperBound()

    """
    Identify the set of unique monomials that occur in objective or constraints
    """

    def uniqueMonomials(self):

        if self.ncon == 0:
            maxd = self.obj_maxd
        else:
            maxd = max(self.obj_maxd, max(self.cons_maxd))

        unique_monomials = {d: [] for d in range(1, maxd + 1)}

        unique_monomials[1] = [(iv,) for iv in range(1, self.nvar + 1)]

        for ic in range(self.ncon + 1):
            if ic < self.ncon:
                mons = self.cons_coeff[ic].keys()
            else:
                mons = self.obj_coeff.keys()
            for mon in mons:
                dmon = len(mon)
                if not (mon in unique_monomials[dmon]):
                    unique_monomials[dmon].append(mon)

        self.maxd = maxd
        self.unique_monomials = unique_monomials

    """
    Find all lower and upper bounds on the variables
    """

    def findLowerUpperBounds(self):

        self.all_var_lb = {}
        self.all_var_ub = {}
        set_all_lb = False
        set_all_ub = False
        if all(self.var_lb == np.zeros((self.nvar, 1))):
            self.all_var_lb = {mon: 0.0 for mon in self.all_nodes_uld}
            set_all_lb = True
        if all(self.var_ub == np.ones((self.nvar, 1))):
            self.all_var_ub = {mon: 1.0 for mon in self.all_nodes_uld}
            set_all_ub = True

        # if all bounds are 0 and 1 return
        if set_all_lb and set_all_ub:
            self.zerone_bdns = True
            return

        self.zerone_bdns = False

        # compute the lower and upper bounds by progressing up the unified linearization diagram
        for d in range(self.maxd):
            mons_d = [mon for mon in self.all_nodes_uld if len(mon) == d + 1]
            # degree-1 monomials
            if d == 0:
                self.all_var_lb = {mon: self.var_lb[mon[0] - 1] for mon in mons_d}
                self.all_var_ub = {mon: self.var_ub[mon[0] - 1] for mon in mons_d}
            # monomials of degree 2 and higher
            else:
                for mon in mons_d:
                    triples = [t for t in self.all_triples_uld if t[2] == mon]
                    self.all_var_lb[mon] = np.Infinity
                    self.all_var_ub[mon] = -np.Infinity
                    mon1 = (mon[0],)
                    if len(mon[1:]) == 1:
                        mon2 = (mon[1],)
                    else:
                        mon2 = tuple(mon[1:])
                    lb1 = self.all_var_lb[mon1]
                    ub1 = self.all_var_ub[mon1]
                    lb2 = self.all_var_lb[mon2]
                    ub2 = self.all_var_ub[mon2]
                    all_prods = [lb1 * lb2, lb1 * ub2, ub1 * lb2, ub1 * ub2]
                    self.all_var_lb[mon] = min(self.all_var_lb[mon], np.min(all_prods))
                    self.all_var_ub[mon] = max(self.all_var_ub[mon], np.max(all_prods))

    """
    Returns the bit-vector encoding of the monomial mon
    example: n = 5, mon = x1x2x5, deg = [1,1,0,0,1]
    """

    def getDegree(self, mon):

        deg = np.zeros((self.nvar, 1), dtype=int)
        for ik in mon:
            deg[ik - 1] = 1
        return deg

    def deg1LTDeg2(self, deg1, deg2):

        lt = True
        for i in range(self.nvar):
            if deg1[i] > deg2[i]:
                break
            if deg2[i] > deg1[i]:
                lt = False
                break
        return lt

    """
    Returns the monomial corresponding to product between var1 and var 2 with degree representation deg1 and deg2
    deg_mons are the degree representations of the monomials
    example: n = 5, deg1 = [1,1,0,0], deg2 = [0,0,0,1] deg_mons = [[1,1,1,0,0], [1,1,0,1,1],[0,1,1,1,0]]. returns (1,2,5) since x1x2x5 is an intermdediate of x1x2x4x5
    example: n = 5, deg1 = [1,1,0,0], deg2 = [0,0,0,1] deg_mons = [[1,1,1,0,0], [1,1,0,1,0], [0,1,1,1,0]]. returns None since x1x2x5 is not an intermediate monomials
    """

    def getProduct(self, mon1, mon2):

        deg1 = self.getDegree(mon1)
        deg2 = self.getDegree(mon2)
        deg = deg1 + deg2
        sumdeg = sum(deg)[0]
        if sumdeg <= self.maxd and max(deg)[0] == 1:
            flag = False
            for d in range(sumdeg, self.maxd + 1):
                for mon in self.unique_monomials[d]:
                    if all(deg <= self.getDegree(mon)):
                        flag = True
                        break
            if flag:
                inds = np.where(deg == 1)[0] + 1
                if self.deg1LTDeg2(deg1, deg2):
                    return mon1, mon2, tuple(inds)
                else:
                    return mon2, mon1, tuple(inds)

        return None, None, None

    """
    Generate all possible monomials i.e. all nodes in the unified linearization diagram
    """

    def generateNodes(self, monomials=None):

        # populate with unary monomials
        if monomials is None:
            # unified linearization diagram
            nodes = [(i,) for i in range(1, self.nvar + 1)]
            # populate with monomials of higher degrees
            for dmon in range(2, self.maxd + 1):
                for imon in self.unique_monomials[dmon]:
                    for d in range(2, dmon + 1):
                        inds = itt.combinations(imon, d)
                        for ind in inds:
                            if not (ind in nodes):
                                nodes.append(ind)
        else:
            nodes = []
            for monomial in monomials:
                # linearization diagram for monomial
                for i in monomial:
                    if not ((i,) in nodes):
                        nodes.append((i,))
                dmon = len(monomial)
                for d in range(2, dmon + 1):
                    inds = itt.combinations(monomial, d)
                    for ind in inds:
                        if not (ind in nodes):
                            nodes.append(ind)

        return nodes

    """
    Generate all the pairwise products i.e. all the arcs in the unified linearization diagram
    """

    def generateArcsTriples2(self, nodes_list):

        # if 1 > 0:
        #     return null, null

        n_nodes = len(nodes_list)

        arcs = []
        triples = []
        for i in range(n_nodes):
            if len(nodes_list[i]) <= 1:
                continue

            elif len(nodes_list[i]) == 2:
                y1 = (nodes_list[i][0],)
                y2 = (nodes_list[i][1],)
                # check if product y1*y2 exists and add to triples
                ny1, ny2, y3 = self.getProduct(y1, y2)

                # add to arcs
                if not ((y1, y3) in arcs):
                    # print("add",y1,y3)
                    arcs.append((y1, y3))
                if not ((y2, y3) in arcs):
                    # print("add",y2,y3)
                    arcs.append((y2, y3))
                # add to triples
                triples.append((ny1, ny2, y3))

            elif len(nodes_list[i]) == 3:
                Y = [
                    (
                        (
                            nodes_list[i][0],
                            nodes_list[i][2],
                        ),
                        (nodes_list[i][1],),
                    ),
                    (
                        (
                            nodes_list[i][0],
                            nodes_list[i][1],
                        ),
                        (nodes_list[i][2],),
                    ),
                    (
                        (
                            nodes_list[i][1],
                            nodes_list[i][2],
                        ),
                        (nodes_list[i][0],),
                    ),
                ]
                for y1, y2 in Y:
                    # check if product y1*y2 exists and add to triples
                    ny1, ny2, y3 = self.getProduct(y1, y2)

                    # add to arcs
                    if not ((y1, y3) in arcs):
                        # print("add",y1,y3)
                        arcs.append((y1, y3))
                    if not ((y2, y3) in arcs):
                        # print("add",y2,y3)
                        arcs.append((y2, y3))
                    # add to triples
                    triples.append((ny1, ny2, y3))

            elif len(nodes_list[i]) == 4:
                Y = [
                    (
                        (
                            nodes_list[i][0],
                            nodes_list[i][1],
                        ),
                        (
                            nodes_list[i][2],
                            nodes_list[i][3],
                        ),
                    ),
                    (
                        (
                            nodes_list[i][0],
                            nodes_list[i][2],
                        ),
                        (
                            nodes_list[i][1],
                            nodes_list[i][3],
                        ),
                    ),
                    (
                        (
                            nodes_list[i][0],
                            nodes_list[i][3],
                        ),
                        (
                            nodes_list[i][1],
                            nodes_list[i][2],
                        ),
                    ),
                    (
                        (nodes_list[i][0],),
                        (
                            nodes_list[i][1],
                            nodes_list[i][2],
                            nodes_list[i][3],
                        ),
                    ),
                    (
                        (nodes_list[i][1],),
                        (
                            nodes_list[i][0],
                            nodes_list[i][2],
                            nodes_list[i][3],
                        ),
                    ),
                    (
                        (nodes_list[i][2],),
                        (
                            nodes_list[i][0],
                            nodes_list[i][1],
                            nodes_list[i][3],
                        ),
                    ),
                    (
                        (nodes_list[i][3],),
                        (
                            nodes_list[i][0],
                            nodes_list[i][1],
                            nodes_list[i][2],
                        ),
                    ),
                ]
                for y1, y2 in Y:
                    # check if product y1*y2 exists and add to triples
                    ny1, ny2, y3 = self.getProduct(y1, y2)

                    # add to arcs
                    if not ((y1, y3) in arcs):
                        # print("add",y1,y3)
                        arcs.append((y1, y3))
                    if not ((y2, y3) in arcs):
                        # print("add",y2,y3)
                        arcs.append((y2, y3))
                    # add to triples
                    triples.append((ny1, ny2, y3))

        return arcs, triples

    """
    Generate all the pairwise products i.e. all the arcs in the unified linearization diagram
    """

    def generateArcsTriples(self, nodes_list):

        arcs, triples = self.generateArcsTriples2(nodes_list)

        return arcs, triples

        n_nodes = len(nodes_list)

        arcs = []
        triples = []
        for i1 in range(n_nodes):
            y1 = nodes_list[i1]
            for i2 in range(i1 + 1, n_nodes):
                y2 = nodes_list[i2]

                # check if product y1*y2 exists and add to triples
                ny1, ny2, y3 = self.getProduct(y1, y2)

                if y3 is None:
                    continue

                # add to arcs
                if not ((y1, y3) in arcs):
                    arcs.append((y1, y3))
                if not ((y2, y3) in arcs):
                    arcs.append((y2, y3))
                # add to triples
                triples.append((ny1, ny2, y3))

        return arcs, triples

    """
    Compute the upper bounds for the multipliers of the network flow constraints
    """

    def computeUpperBound(self):

        ylist = self.all_nodes_uld
        triples = self.all_triples_uld

        lam1bnd = {t: 0.0 for t in triples}
        lam2bnd = {t: 0.0 for t in triples}

        # upper bound on lam3
        lam3bnd = sum([max(0, -self.obj_coeff[k]) for k in self.obj_coeff.keys()])

        degs = [len(mon) for mon in ylist]
        sorted_index = np.argsort(degs)

        # one inequality for each monomial
        for ind in sorted_index:

            mon = ylist[ind]

            tail1_set = [t for t in triples if mon == t[0]]
            tail2_set = [t for t in triples if mon == t[1]]
            head_set = [t for t in triples if mon == t[2]]

            if tail1_set == [] and tail2_set == []:
                continue

            bnd = lam3bnd + self.all_obj_coeff[mon]

            for t in head_set:
                bnd += lam1bnd[t] + lam2bnd[t]

            for t in tail1_set:
                lam1bnd[t] = bnd
            for t in tail2_set:
                lam2bnd[t] = bnd

        return lam1bnd, lam2bnd, lam3bnd

    """
    Adds McCormick inequalities for var3 = var1*var2 to the lp
    """

    def addMcCormick(self, var1, var2, var3, y, lp):

        if self.zerone_bdns:
            lp.addConstr(y[var3] <= y[var1], name="McC1" + str(var1) + "_" + str(var2) + "_" + str(var3))
            lp.addConstr(y[var3] <= y[var2], name="McC2" + str(var1) + "_" + str(var2) + "_" + str(var3))
            lp.addConstr(y[var3] >= y[var1] + y[var2] - 1, name="McC3" + str(var1) + "_" + str(var2) + "_" + str(var3))
        else:
            lb1 = self.all_var_lb[var1]
            ub1 = self.all_var_ub[var1]
            lb2 = self.all_var_lb[var2]
            ub2 = self.all_var_ub[var2]

            lp.addConstr(
                y[var3] - lb1 * y[var2] - lb2 * y[var1] + lb1 * lb2 >= 0,
                name="McC1" + str(var1) + "_" + str(var2) + "_" + str(var3),
            )
            lp.addConstr(
                y[var3] - ub1 * y[var2] - ub2 * y[var1] + ub1 * ub2 >= 0,
                name="McC2" + str(var1) + "_" + str(var2) + "_" + str(var3),
            )
            lp.addConstr(
                y[var3] - lb1 * y[var2] - ub2 * y[var1] + lb1 * ub2 <= 0,
                name="McC3" + str(var1) + "_" + str(var2) + "_" + str(var3),
            )
            lp.addConstr(
                y[var3] - ub1 * y[var2] - lb2 * y[var1] + ub1 * lb2 <= 0,
                name="McC4" + str(var1) + "_" + str(var2) + "_" + str(var3),
            )

    """
    Forms and solves the root node LP relaxation for a given linearization using scipy
        lin - choice of linearization, if None then use all arcs
        print_flag - controls Gurobi output
    """

    def solveLPscipy(self, lin=None, print_flag=0):

        # n = self.nvar        # number of variables

        # m = len(self.coeff)  # number of monomials

        if self.all_nodes_uld is None:
            self.generate_all_monomials()
        ylist = self.all_nodes_uld

        # generate all inequalities
        if lin is None:
            if self.all_triples_uld is None:
                self.generate_all_products()
            lin = self.all_triples_uld
            print()

        n_ylist = len(ylist)
        n_arcs = len(lin)

        c = np.zeros((n_ylist, 1))
        for mon in ylist:
            c[ylist.index(mon)] = self.all_obj_coeff[mon]
        Ain = np.zeros((3 * n_arcs, n_ylist))
        bin = np.zeros((3 * n_arcs, 1))
        for ic in range(n_arcs):
            triple = lin[ic]
            ind_y1 = ylist.index(triple[0])
            ind_y2 = ylist.index(triple[1])
            ind_y3 = ylist.index(triple[2])
            Ain[3 * ic, ind_y1] = -1.0
            Ain[3 * ic, ind_y3] = 1.0
            Ain[3 * ic + 1, ind_y2] = -1.0
            Ain[3 * ic + 1, ind_y3] = 1.0
            Ain[3 * ic + 2, ind_y1] = 1.0
            Ain[3 * ic + 2, ind_y2] = 1.0
            Ain[3 * ic + 2, ind_y3] = -1.0
            bin[3 * ic + 2, 0] = 1.0
        bounds = (0, 1)

        # res = sp.optimize.linprog(c,A_ub=Ain,b_ub=bin,bounds=bounds)#,method='interior-point')
        res = sp.linprog(c, A_ub=Ain, b_ub=bin, bounds=bounds)  # ,method='interior-point')

        lp_results = struct()
        lp_results.nvar = n_ylist
        lp_results.ncon = 3 * n_arcs
        lp_results.obj = res.fun
        lp_results.cput = 0.0
        lp_results.solved = res.success

        if self.obj_sense == "Max":
            lp_results.obj = -1.0 * lp_results.obj
        lp_results.obj = lp_results.obj + self.obj_offset

        return lp_results

    """
    Return linearization with all triples
    """

    def fullLinearization(self):
        lin = None
        if self.all_triples_uld is None:
            self.all_arcs_uld, self.all_triples_uld = self.generateArcsTriples(self.all_nodes_uld)
        lin = self.all_triples_uld
        return lin

    """
    Heuristic linearization
        mlp - MLP object
        lin_type - "sequential"/"greedy"/"pairwise"
    """

    def heuristicLinearization(self, lin_type):

        lin = []
        nvar = self.nvar
        vars_indices = {(i + 1,): i for i in range(nvar)}
        maxd = self.maxd
        monomials = []
        for d in range(2, maxd + 1):
            for mon in self.unique_monomials[d]:
                monst = struct()
                monst.mon = mon
                monst.degree = list(np.reshape(self.getDegree(mon), (nvar,)))
                monomials.append(monst)

        # freq_pairs = { ((i1,),(i2,)) : 0 for i1 in range(nvar-1) for i2 in range(i1,nvar) }
        pairs_freq = {}

        for monst in monomials:
            mon = monst.mon
            nmon = len(mon)
            for i1 in range(nmon - 1):
                for i2 in range(i1 + 1, nmon):
                    y1 = mon[i1]
                    y2 = mon[i2]
                    if ((y1,), (y2,)) in pairs_freq.keys():
                        pairs_freq[((y1,), (y2,))] = pairs_freq[((y1,), (y2,))] + 1
                    else:
                        pairs_freq.update({((y1,), (y2,)): 1})

        while not (len(monomials) == 0):

            if lin_type == "greedy":
                maxtup = None
                maxval = 0
                for pair in pairs_freq.keys():
                    val = pairs_freq[pair]
                    if val > maxval:
                        maxtup = pair
                        maxval = val
                var1 = maxtup[0]
                var2 = maxtup[1]
            else:
                # identify the pair to linearize
                var1 = monomials[0].mon[0]
                var2 = monomials[0].mon[1]
            inds2remove = []

            if not (type(var1) == tuple):
                var1 = (var1,)
            if not (type(var2) == tuple):
                var2 = (var2,)

            y1, y2, y3 = self.getProduct(var1, var2)
            if y3 is None:
                raise Exception("This should not happen!")
            # new triple in the linearization
            lin.append((y1, y2, y3))
            # update the variables
            vars_indices.update({y3: nvar})
            nvar = nvar + 1

            ind1 = vars_indices[y1]
            ind2 = vars_indices[y2]

            # scan the monomials and replace the occurences of y1,y2 with y3
            # modify the .mon field by removing y1,y2 and replacing with y3
            # modify the .degree field by appending another entry for y3, removing y1,y2
            # if mon is of length 1, mark the monomial for deletion
            cnt = 0
            for monst in monomials:

                if monst.degree[ind1] == 1 and monst.degree[ind2] == 1:

                    # update pairs_freq dict
                    if lin_type == "greedy":
                        nmon = len(monst.mon)
                        for i1 in range(nmon - 1):
                            for i2 in range(i1 + 1, nmon):
                                yi1 = monst.mon[i1]
                                yi2 = monst.mon[i2]
                                if not (type(yi1) == tuple):
                                    yi1 = (yi1,)
                                if not (type(yi2) == tuple):
                                    yi2 = (yi2,)
                                if yi1 == y1 or yi1 == y2 or yi2 == y1 or yi2 == y2:
                                    nyi1, nyi2, nyi1i2 = self.getProduct(yi1, yi2)
                                    pairs_freq[(nyi1, nyi2)] = pairs_freq[(nyi1, nyi2)] - 1

                        for i1 in range(nmon):
                            yi1 = monst.mon[i1]
                            if not (type(yi1) == tuple):
                                yi1 = (yi1,)
                            if not (yi1 == y1 or yi1 == y2):
                                nyi1, ny3, nyi13 = self.getProduct(yi1, y3)
                                if (nyi1, ny3) in pairs_freq:
                                    pairs_freq[(nyi1, ny3)] = pairs_freq[(nyi1, ny3)] + 1
                                else:
                                    pairs_freq.update({(nyi1, ny3): 1})

                    monst.degree.append(1)
                    monst.degree[ind1] = 0
                    monst.degree[ind2] = 0
                    if lin_type == "sequential":
                        mon = [y3]
                    else:
                        mon = []
                    for it in monst.mon:
                        tupit = it
                        if not (type(it) == tuple):
                            tupit = (it,)
                        if not (tupit == var1 or tupit == var2):
                            mon.append(it)
                    if not (lin_type == "sequential"):
                        mon.append(y3)
                    monst.mon = mon
                    if len(mon) == 1:
                        inds2remove.append(cnt)
                else:
                    monst.degree.append(0)

                cnt = cnt + 1

            # pop the monomials to remove in reverse order
            for i in range(len(inds2remove) - 1, -1, -1):
                monomials.pop(inds2remove[i])

        return lin

    """
    Read linearization from file
        dir - folder name
        fname - file name
    """

    def readLinearization(self, dir, fname):

        fp = open(dir + fname, "r")

        lines = fp.readlines()

        ntriples = int(lines[0])

        lin = []

        for i in range(1, ntriples + 1):

            triplei = []
            start = 0
            while start < len(lines[i]) - 1:
                left_inds = lines[i][start:].index("(") + start
                right_inds = lines[i][start:].index(")") + start
                tup = eval(lines[i][left_inds : right_inds + 1])
                triplei.append(tup)
                start = right_inds + 1
            lin.append((triplei[0], triplei[1], triplei[2]))

        fp.close()

        return lin
