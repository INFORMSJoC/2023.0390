# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# ----------------------------------------------------------------------
# Run instance
# ----------------------------------------------------------------------

import os
import sys

path = os.path.abspath("../")
sys.path.append(path)

import sys

import src.MLP as MLP
import src.optimalLinearization as optl
import src.solvePO as PO
import src.solveQCP as QCP

# Parameters

# --------------------------------------------------------------------------------------------------------------


def main(argv):
    if len(sys.argv) < 2:
        print("\tsimple_tests.py [filename]")
        sys.exit()

    full_time = 300
    accumulated_time = 0
    filename = sys.argv[1]

    print("Reading file...")
    mlp = MLP.MLP(mlp_dir="", mlp_file=filename)

    # Sequential linearization
    print("Sequential Linearization")
    seq_lin = mlp.heuristicLinearization("sequential")
    seq_results = mlp.solveLP(lin=seq_lin, time_limit=30)
    seq_qcp_results = QCP.solveQCP_pyomo_gurobi(mlp, full_time, seq_lin)

    print(
        str(len(seq_lin))
        + ", "  # size of linearization
        + str(seq_results.cput)
        + ", "  # time to solve LP relaxation
        + str(seq_results.obj)
        + ", "  # LP bound
        + str(seq_qcp_results.runtime)
        + ", "  # time to solve QCP
        + str(seq_qcp_results.objective)
        + ", "  # solution of QCP
        + str(seq_qcp_results.lb)
        + ", "  # LB of QCP
        + str(seq_qcp_results.ub)
        + ", "
    )  # UB of QCP

    # Greedy linearization
    print("Greedy Linearization")
    lin = mlp.heuristicLinearization("greedy")
    greedy_results = mlp.solveLP(lin=lin, time_limit=30)
    greedy_qcp_results = QCP.solveQCP_pyomo_gurobi(mlp, full_time, lin)

    print(
        str(len(lin))
        + ", "  # size of linearization
        + str(greedy_results.cput)
        + ", "  # time to solve LP relaxation
        + str(greedy_results.obj)
        + ", "  # LP bound
        + str(greedy_qcp_results.runtime)
        + ", "  # time to solve QCP
        + str(greedy_qcp_results.objective)
        + ", "  # solution of QCP
        + str(greedy_qcp_results.lb)
        + ", "  # LB of QCP
        + str(greedy_qcp_results.ub)
        + ", "
    )  # UB of QCP

    # Full linearization
    print("Full Linearization")
    full = mlp.Linearization("full")
    lp_results = mlp.solveLP(lin=full, time_limit=30)
    print(
        lp_results.cpu,
        min(lp_results.obj, lp_results.objbnd),
        max(lp_results.obj, lp_results.objbnd),
        abs((lp_results.obj - lp_results.objbnd) / abs(lp_results.obj)),
        lp_results.obj,
    )

    # Min-size linearization
    print("Minimum Size MLP")
    minlin = optl.optimalLinearization(mlp, "minlin", time_limit=30, init_lin=lin)
    accumulated_time = minlin.cput
    if minlin.nsolns > 0:
        lin = minlin.best_lins[0]
    mlp_results = mlp.solveLP(lin=lin, time_limit=30)
    mlp_qcp_results = QCP.solveQCP_pyomo_gurobi(mlp, full_time - accumulated_time, lin)
    print(
        str(minlin.cput),  # time to solve MLP
        str(minlin.objbnd),  # lb
        str(minlin.obj),  # ub
        str(mlp_results.cput),  # time to solve LP relaxation
        str(mlp_results.obj),  # LP bound
        str(mlp_qcp_results.runtime),  # time to solve QCP
        str(mlp_qcp_results.objective),  # solution of QCP
        str(mlp_qcp_results.lb),  # LB of QCP
        str(mlp_qcp_results.ub),
    )  # UB of QCP

    # Best bound MLP
    print("Best Bound MLP")
    min_triples = int(round(minlin.obj))
    bestbnd = optl.optimalLinearization(mlp, "bestbnd", max_triples=min_triples, time_limit=30, init_lin=lin)
    if bestbnd.nsolns > 0:
        lin = bestbnd.best_lins[0]
    bestbnd_results = mlp.solveLP(lin=lin, time_limit=30)
    accumulated_time += bestbnd_results.cput
    bestbnd_qcp_results = QCP.solveQCP_pyomo_gurobi(mlp, full_time - accumulated_time, lin)
    print(
        str(bestbnd_results.cput),  # size of linearization
        str(bestbnd_results.obj),  # time to solve MLP
        str(bestbnd_qcp_results.runtime),  # time to solve QCP
        str(bestbnd_qcp_results.objective),  # solution of QCP
        str(bestbnd_qcp_results.lb),  # LB  of QCP
        str(bestbnd_qcp_results.ub),
    )  # UB  of QCP

    print("BARON with CPLEX")
    po_results = PO.solvePO_nolin(mlp, time_limit=full_time)
    print(
        "Time="
        + str(po_results.runtime)
        + ", Obj="
        + str(po_results.objective)
        + ", LB="  # solution of QCP
        + str(po_results.lb)
        + ", UB="  # LB of QCP
        + str(po_results.ub)
    )  # UB of QCP

    return


if __name__ == "__main__":
    main(sys.argv[1:])
