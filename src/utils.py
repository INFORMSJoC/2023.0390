# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import numpy as np


class struct:
    pass


"""
Creates and populates a structure with Gurobi results
    model - Gurobi model object
    ismip - Boolean flag indicating if the model is a MIP
    ismip - Boolean flag indicating if the model is a QCP
"""


def gurobiResults(model, ismip, isqcp=0):

    results = struct()

    results.nvar = model.getAttr("NumVars")
    results.ncon = model.getAttr("NumConstrs")
    if isqcp == 1:
        results.ncon = results.ncon + model.getAttr("NumQConstrs")
    results.obj = model.getAttr("ObjVal")
    if ismip == 1:
        results.nbinvar = model.getAttr("NumBinVars")
        results.objbnd = model.getAttr("ObjBound")
        results.mipgap = model.getAttr("MIPGap")
        results.nnodes = model.getAttr("NodeCount")
    results.cput = model.getAttr("Runtime")

    return results
