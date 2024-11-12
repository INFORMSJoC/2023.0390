# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import time

import gurobipy as gb
import mock
import numpy as np
import pyomo.environ as pyo
from pyomo.core import TransformationFactory
from pyomo.opt import SolverFactory, SolverManagerFactory

import src.utils as ut

"""
Solve original problem WITHOUT LINEARIZATIONS using Pyomo, BARON, and CPLEX
"""


def solvePO_nolin(mlp, time_limit):

    qcp = pyo.ConcreteModel()

    qcp.b = pyo.VarList(domain=pyo.Binary)
    qcp.c = pyo.VarList(domain=pyo.NonNegativeReals)
    var_dict = dict()

    # Add original variables
    is_discrete = False
    for var in range(mlp.nvar):
        if mlp.var_type[var] == "B":
            var_dict[(var + 1,)] = qcp.b.add()
            is_discrete = True
        else:
            var_dict[(var + 1,)] = qcp.c.add()
        var_dict[(var + 1,)].setlb(mlp.var_lb[var][0])
        var_dict[(var + 1,)].setub(mlp.var_ub[var][0])

    # Add objective
    obj = mlp.obj_offset
    for mon in mlp.obj_coeff.keys():
        expr = 1
        for i in mon:
            expr = expr * var_dict[(i,)]
        obj = obj + mlp.obj_coeff[mon] * expr

    # Objective sense
    if mlp.obj_sense == "Min":
        qcp.objective = pyo.Objective(expr=obj, sense=pyo.minimize)
    else:
        qcp.objectives = pyo.Objective(expr=obj, sense=pyo.maximize)

    # Add original constraints in MLP
    qcp.constraints = pyo.ConstraintList()
    for ic in range(mlp.ncon):
        qcp.constraints.add(
            sum([mlp.cons_coeff[ic][mon] * var_dict[mon] for mon in mlp.cons_coeff[ic].keys()]) <= mlp.cons_ub[ic]
        )

    output = mock.Mock()

    # Solve linear relaxation first
    solver = SolverFactory("baron")

    # Solve original model
    results = solver.solve(
        qcp,
        tee=False,
        options={
            "MaxTime": time_limit,
            "LPSol": 3,
            "CplexLibName": "/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx/libcplex2211.dylib",
        },
    )

    # Extract results
    output.runtime = results["Solver"][0]["Time"]
    output.objective = 0
    if mlp.obj_sense == "Min":
        output.objective = results.Problem[0].upper_bound
    else:
        output.objective = results.Problem[0].lower_bound
    output.lb = results.Problem[0].lower_bound
    output.ub = results.Problem[0].upper_bound

    return output
