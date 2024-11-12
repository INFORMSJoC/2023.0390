# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import gurobipy as gb
import mock
import numpy as np
import pyomo.environ as pyo
from pyomo.core import TransformationFactory
from pyomo.opt import SolverFactory

import src.utils as ut

"""
Solve the QCP using Pyomo and Gurobi
"""


def solveQCP_pyomo_gurobi(mlp, time_limit, lin):

    qcp = pyo.ConcreteModel()

    qcp.b = pyo.VarList(domain=pyo.Binary)
    qcp.c = pyo.VarList(domain=pyo.NonNegativeReals)
    var_dict = dict()

    list_monomials = mlp.all_nodes_uld

    # Create variables
    is_discrete = False
    for mon in mlp.all_nodes_uld:
        # Binary variable
        if all([mlp.var_type[i - 1] == "B" for i in mon]):
            is_discrete = True
            var_dict[mon] = qcp.b.add()
        # Continuous variable
        else:
            var_dict[mon] = qcp.c.add()
        var_dict[mon].setlb(mlp.all_var_lb[mon])
        var_dict[mon].setub(mlp.all_var_ub[mon])

    qcp.constraints = pyo.ConstraintList()

    # Add quadratic constraints
    for t in lin:
        qcp.constraints.add(var_dict[t[0]] * var_dict[t[1]] == var_dict[t[2]])

    # Add original constraints in MLP
    for ic in range(mlp.ncon):
        qcp.constraints.add(
            sum([mlp.cons_coeff[ic][mon] * var_dict[mon] for mon in mlp.cons_coeff[ic].keys()]) <= mlp.cons_ub[ic]
        )

    # Add objective
    obj = mlp.obj_offset
    for mon in mlp.obj_coeff.keys():
        obj = obj + mlp.obj_coeff[mon] * var_dict[mon]

    if mlp.obj_sense == "Min":
        qcp.objective = pyo.Objective(expr=obj, sense=pyo.minimize)
    else:
        qcp.objectives = pyo.Objective(expr=obj, sense=pyo.maximize)

    qcp.z = pyo.Var(domain=pyo.NonNegativeReals)
    for mon in mlp.obj_coeff.keys():
        obj = obj + mlp.obj_coeff[mon] * var_dict[mon]

    solver = SolverFactory("gurobi")

    try:
        results = solver.solve(qcp, tee=False, options={"TimeLimit": time_limit, "NonConvex": 2})  # GUROBI

    except Exception as e:
        print("Exception", e)

    output = mock.Mock()
    output.runtime = results["Solver"][0]["Time"]
    output.objective = 0
    if mlp.obj_sense == "Min":
        output.objective = results.Problem[0]["Upper bound"]
    else:
        output.objective = results.Problem[0]["Lower bound"]
    output.lb = results.Problem[0]["Lower bound"]
    output.ub = results.Problem[0]["Upper bound"]

    return output
