# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import itertools
from collections import Counter

import gurobi_logtools as glt
import gurobipy as gb
import numpy as np

import src.utils as ut


def mycallback(model, where):

    if where != gb.GRB.Callback.MIPNODE:
        return

    nodecount = model.cbGet(gb.GRB.Callback.MIPNODE_NODCNT)
    if model._firstcall:
        if model._objtype == "minlin":
            v = model.cbGetNodeRel(model._v)
        if model._objtype == "bestbnd":
            lam3 = model.cbGetNodeRel(model._lam3)
            mu = {mon: model.cbGetNodeRel(model._mu[mon]) for mon in model._ylist}

        objval = 0.0
        if model._objtype == "minlin":
            for t in model._mlp.all_triples_uld:
                objval += v[t]
        if model._objtype == "bestbnd":
            for t in model._triples:
                objval += -lam3[t]
            for mon in model._ylist:
                objval += -mu[mon]

        model._rootobj = objval
        model._firstcall = False


"""
Find the minimum linearization for the MLP
    objtype - "minlin" or "bestbnd"
    max_triple - maximum number of terms in the linearization
    time_limit - time limit
    num_solns -  number of optimal solutions
    grb_opt_emphasis - True: optimalilty emphasis, False: feasibility emphasis
    init_lin - initial linearization
    print_flag - Gurobi output flag
    num_threads - Gurobi flag for number of threads to use
"""


def optimalLinearization(
    mlp,
    objtype,
    max_triples=0,
    time_limit=60,
    num_solns=100,
    grb_opt_emphasis=0,
    init_lin=None,
    print_flag=0,
    num_threads=1,
):

    mip = gb.Model()
    mip.Params.LogToConsole = 0  # Disable logging to console
    mip.setParam("TimeLimit", time_limit)

    if objtype == "bestbnd" and max_triples == 0:
        raise Exception("Max number of terms must be grater than 0")

    # variables
    u = {}
    for d in range(2, mlp.maxd + 1):
        for mon in mlp.unique_monomials[d]:
            u.update({mon: {}})
            for t in mlp.monomial_triples_ld[mon]:
                u[mon].update({t: mip.addVar(obj=0.0, lb=0.0, ub=1.0, name="u" + str(mon) + str(t))})
    v = mip.addVars(mlp.all_triples_uld, obj=0.0, vtype=gb.GRB.BINARY, name="v")

    filtered_triples = []
    if (objtype == "minlin") and mlp.maxd == 4:

        # Set of isolated triples for monomials of size 4
        isolated_triples = list()
        for mon in mlp.unique_monomials[4]:
            isolated_triples += list(itertools.combinations(mon, 3))

        # Count the occurrences of each element using Counter
        element_counts = Counter(isolated_triples)

        # Triples that appear only once in the set of monomials of size 4
        filtered_triples = [triple for triple in isolated_triples if element_counts[triple] == 1]

        # Intermediate triples (i.e., triples that are not monomials)
        intermediate_triples = list(set(filtered_triples))
        intermediate_triples = [t for t in intermediate_triples if not (t in mlp.unique_monomials[3])]

        # If intermediate triple is selected, it must be chosen at least twice
        head_set = [t for t in mlp.all_triples_uld if t[2] in intermediate_triples]

        # One inequality per triple
        for t in head_set:
            expr = gb.LinExpr(0.0)
            for mon in mlp.unique_monomials[4]:
                if t[2] in mlp.monomial_nodes_ld[mon]:
                    expr += u[mon][t]
            mip.addConstr(2 * v[t] >= expr, name="triple_usage_" + str(t))

        # Remove linearizations of isolated triples
        for mon_lin in mlp.unique_monomials[4]:
            for mon in mlp.monomial_nodes_ld[mon_lin]:

                if len(mon) == 3:
                    head_set = [t for t in mlp.monomial_triples_ld[mon_lin] if mon == t[2]]
                    # If head is an isolated triple that appears only once, filter it out
                    if mon in filtered_triples and mon in intermediate_triples:
                        mip.addConstr(
                            sum([v[t] for t in head_set]) == 0, name="filtered_v_" + str(mon_lin) + "_" + str(mon)
                        )
                        mip.addConstr(
                            sum([u[mon_lin][t] for t in head_set]) == 0,
                            name="filtered_u_" + str(mon_lin) + "_" + str(mon),
                        )

    # add the network flow model for each monomial
    for d in range(2, mlp.maxd + 1):
        for mon in mlp.unique_monomials[d]:
            addMonomialLinearizationModel(mip, mlp, u, v, mon, filtered_triples, objtype)

    if objtype == "bestbnd":
        if mlp.ncon > 0:
            # consider only the objective in the case of constrained MLP
            ylist = mlp.obj_nodes_uld
            triples = mlp.obj_triples_uld
        else:
            ylist = mlp.all_nodes_uld
            triples = mlp.all_triples_uld

        # multipliers for y <= 1 and linearization equalities
        mu = {}
        for mon in ylist:
            mu.update({mon: mip.addVar(lb=0.0, obj=0.0, name="mu" + str(mon))})
        lam1 = mip.addVars(triples, lb=0.0, obj=0.0, name="lam1")
        lam2 = mip.addVars(triples, lb=0.0, obj=0.0, name="lam2")
        lam3 = mip.addVars(triples, lb=0.0, obj=0.0, name="lam3")

        # 1 inequality for each monomial
        for mon in ylist:

            tail1_set = [t for t in triples if mon == t[0]]
            tail2_set = [t for t in triples if mon == t[1]]
            head_set = [t for t in triples if mon == t[2]]

            expr = gb.LinExpr(0.0)
            for t in tail1_set:
                expr += -lam1[t] + lam3[t]
            for t in tail2_set:
                expr += -lam2[t] + lam3[t]
            for t in head_set:
                expr += lam1[t] + lam2[t] - lam3[t]
            expr += mu[mon] + mlp.all_obj_coeff[mon]

            mip.addConstr(expr >= 0)

        # big-M constraints
        for t in triples:
            mip.addConstr(lam1[t] <= mlp.lam1bnd[t] * v[t])
            mip.addConstr(lam2[t] <= mlp.lam2bnd[t] * v[t])
            mip.addConstr(lam3[t] <= mlp.lam3bnd * v[t])

        # limit number of terms that are chosen
        mip.addConstr(sum([v[t] for t in triples]) <= max_triples)

    # objective
    obj = gb.LinExpr(0.0)
    if objtype == "minlin":
        for t in mlp.all_triples_uld:
            obj += v[t]
        mip.setObjective(obj, gb.GRB.MINIMIZE)
    elif objtype == "bestbnd":
        for mon in ylist:
            obj += -mu[mon]
        for t in triples:
            obj += -lam3[t]
        mip.setObjective(obj, gb.GRB.MAXIMIZE)

    # set initial linearization
    for t in mlp.all_triples_uld:
        v[t].start = 0.0
    if not (init_lin is None):
        for t in init_lin:
            v[t].start = 1.0

    mip._firstcall = True
    mip._rootobj = 0.0
    mip._mlp = mlp
    mip._objtype = objtype
    mip._v = v
    if objtype == "bestbnd":
        mip._ylist = ylist
        mip._triples = triples
        mip._lam3 = lam3
        mip._mu = mu

    mip.update()
    mip.Params.LogFile = "gurobi.log"

    log_file = open(mip.Params.LogFile, "w")

    mip.optimize(mycallback)
    results = glt.parse(["gurobi.log"])
    nodelog_progress = results.progress("nodelog")
    results = ut.gurobiResults(mip, 1)
    nsolns = mip.getAttr("SolCount")

    lins = []
    for i in range(nsolns):
        mip.setParam("SolutionNumber", i)
        lin = [t for t in mlp.all_triples_uld if v[t].getAttr("Xn") >= 0.5]
        lins.append(lin)

    results = ut.gurobiResults(mip, 1)
    results.obj = mip.objVal
    results.best_lins = lins
    results.objbnd = mip.obj_bound
    results.nsolns = nsolns
    results.cpu = mip.Runtime

    results.primal_list = nodelog_progress["Incumbent"].to_list()
    results.dual_list = nodelog_progress["BestBd"].to_list()
    results.gap_list = nodelog_progress["Gap"].to_list()
    results.time_list = nodelog_progress["Time"].to_list()

    return results


# Form network flow formulation of DD for a monomial
def addMonomialLinearizationModel(mip, mlp, u, v, mon_lin, objtype="minlin", filtered_triples=[]):

    triples = mlp.monomial_triples_ld[mon_lin]

    # Store all triples associated with this monomial
    mon_triples = list()
    non_isolated_triple = False

    for mon in mlp.monomial_nodes_ld[mon_lin]:

        if len(mon) == 1:
            continue

        tail1_set = [t for t in triples if mon == t[0]]
        tail2_set = [t for t in triples if mon == t[1]]
        head_set = [t for t in triples if mon == t[2]]

        mon_triples += tail1_set
        mon_triples += tail2_set
        mon_triples += head_set

        expr = gb.LinExpr(0.0)
        expr += sum([u[mon_lin][t] for t in tail1_set])
        expr += sum([u[mon_lin][t] for t in tail2_set])
        expr += -sum([u[mon_lin][t] for t in head_set])

        if mon == mon_lin:
            expr += 1
            if objtype == "minlin":
                mip.addConstr(expr == 0, name="ld_" + str(mon_lin) + "_" + str(mon))
            else:
                mip.addConstr(expr <= 0, name="ld_" + str(mon_lin) + "_" + str(mon))
        else:
            mip.addConstr(expr == 0, name="ld_" + str(mon_lin) + "_" + str(mon))

    for t in triples:
        mip.addConstr(u[mon_lin][t] <= v[t])

    mip.update()
