# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# ----------------------------------------------------------------------
# Run instance
# ----------------------------------------------------------------------

import json
import os
import pickle
import resource
import sqlite3
import subprocess
import sys
import time
import timeit
from math import *
from os import listdir
from os.path import isfile, join

import pandas as pd
import psutil

path = os.path.abspath("../")
sys.path.append(path)

import sys

import src.MLP as MLP
import src.optimalLinearization as optl
import src.solvePO as PO
import src.solveQCP as QCP

# Create wrapper functions using lambda functions
wrapper_baron = lambda param1, param2, param3: PO.solvePO_nolin(param1, param2)
wrapper_gurobi10 = lambda param1, param2, param3: QCP.solveQCP_pyomo_gurobi(param1, param2, param3)


def get_instance(cursor, conn, mlp_dir, mlp_file, type):

    mlp = None

    cursor.execute(
        """
        SELECT Instance.*
            FROM Instance
            WHERE Instance.name = ?
        """,
        (mlp_file,),
    )
    # Extract the ID from the result tuple
    # Fetch the result (in this case, fetch one row)
    result = cursor.fetchone()
    if result:
        instance_id = result[0]  # Assuming ID is the first column in the SELECT query result
    else:
        # Read MLP
        mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file, mlp_only=True)

        # Entry doesn't exist, add a new entry
        print("Instance", mlp_file)
        cursor.execute(
            "INSERT INTO Instance (name, type, n, m) VALUES (?, ?, ?, ?)",
            (mlp_file, type, mlp.nvar, len(mlp.obj_coeff)),
        )
        # Retrieve the ID of the newly added entry
        cursor.execute("SELECT last_insert_rowid()")
        instance_id = cursor.fetchone()[0]

        # Commit the transaction
        conn.commit()

        # Incomplete MLP description, so set it to None
        mlp = None

    return mlp, instance_id


# Sequential and greedy - no solver needed
def linearization(cursor, conn, mlp, instance_id, lin_type, mlp_dir, mlp_file):

    lin = None
    experiment_id = -1
    accumulated_time = 0
    cursor.execute(
        """
        SELECT BaseLinearization.*
            FROM BaseLinearization
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND BaseLinearization.linearization = ?
        """,
        (mlp_file, lin_type),
    )
    result = cursor.fetchone()
    if result is None:
        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        # Obtain sequential linearization
        print("\t\tLinearization", mlp_file, lin_type)
        lin = mlp.Linearization(lin_type)
        blob = pickle.dumps(lin)

        cursor.execute(
            "INSERT INTO BaseLinearization   (instance_id, linearization, ub,triples ) VALUES (?, ?, ?, ?)",
            (instance_id, lin_type, len(lin), blob),
        )
        cursor.execute("SELECT last_insert_rowid()")
        linearization_id = cursor.fetchone()[0]

        # Commit the transaction
        conn.commit()

    else:
        linearization_id = result[0]
        # Retrieve linearization from the database
        lin = pickle.loads(result[4])

    return mlp, lin, linearization_id, 0


def minimum_linearization(cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver, warm_start):

    cursor.execute(
        """
        SELECT BaseLinearization.*, MinLinearization.*
            FROM MinLinearization
            JOIN BaseLinearization ON MinLinearization.id = BaseLinearization.id
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND MinLinearization.solver = ? AND MinLinearization.warm_start = ? AND BaseLinearization.linearization = 'minlin'
        """,
        (mlp_file, solver, warm_start),
    )
    result = cursor.fetchone()
    if result is None:

        # Retrieve greedy linearization
        base_lin_id = None
        lin = None
        if warm_start:
            # print("warmstart = ",warm_start)
            mlp, lin, base_lin_id, _ = linearization(cursor, conn, mlp, instance_id, "greedy", mlp_dir, mlp_file)

        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        # Obtain sequential linearization
        print("\t\tMinimum Linearization", mlp_file, solver, warm_start)
        minlin = optl.optimalLinearization(mlp, "minlin", time_limit=30, print_flag=0, init_lin=lin)

        # Replace greedy linearization for minlin
        if minlin.nsolns > 0:
            lin = minlin.best_lins[0]

        # Solution gap
        gap = abs((minlin.obj - ceil(minlin.objbnd))) / abs(minlin.obj)

        # Update accumulated time
        accumulated_time = minlin.cpu
        blob = pickle.dumps(lin)

        cursor.execute(
            "INSERT INTO BaseLinearization   (instance_id, linearization, ub,triples ) VALUES (?, ?, ?, ?)",
            (instance_id, "minlin", minlin.obj, blob),
        )
        cursor.execute("SELECT last_insert_rowid()")
        lin_id = cursor.fetchone()[0]

        cursor.execute(
            "INSERT INTO MinLinearization   (id,     solver, warm_start, base_lin_id,  runtime, lb, gap ) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (lin_id, solver, warm_start, base_lin_id, accumulated_time, minlin.objbnd, gap),
        )
        cursor.execute("SELECT last_insert_rowid()")
        lin_id = cursor.fetchone()[0]

        # Commit the transaction
        conn.commit()

    else:
        lin_id = result[0]
        accumulated_time = result[4 + 1 + 4]
        lin = pickle.loads(result[4])

    return mlp, lin, lin_id, accumulated_time


def ev_minimum_linearization(cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver):

    warm_start = 0

    cursor.execute(
        """
        SELECT BaseLinearization.*, MinLinearization.*
            FROM MinLinearization
            JOIN BaseLinearization ON MinLinearization.id = BaseLinearization.id
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND MinLinearization.solver = ? AND MinLinearization.warm_start = ? AND BaseLinearization.linearization = 'minlin0'
        """,
        (mlp_file, solver, warm_start),
    )
    result = cursor.fetchone()
    if result is None:

        # Retrieve greedy linearization
        base_lin_id = None
        lin = None

        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        # Obtain sequential linearization
        print("\t\tMinimum Linearization without warm start", mlp_file, solver, warm_start)
        minlin = optl.optimalLinearization(mlp, "minlin", time_limit=60, print_flag=0, init_lin=lin)

        # Replace greedy linearization for minlin
        if minlin.nsolns > 0:
            lin = minlin.best_lins[0]

        # Solution gap
        gap = abs((minlin.obj - ceil(minlin.objbnd))) / abs(minlin.obj)

        # Update accumulated time
        accumulated_time = minlin.cpu
        blob = pickle.dumps(lin)

        cursor.execute(
            "INSERT INTO BaseLinearization   (instance_id, linearization, ub,triples ) VALUES (?, ?, ?, ?)",
            (instance_id, "minlin0", minlin.obj, blob),
        )
        cursor.execute("SELECT last_insert_rowid()")
        lin_id = cursor.fetchone()[0]

        cursor.execute(
            "INSERT INTO MinLinearization   (id,     solver, warm_start, runtime, lb, gap ) VALUES (?, ?, ?, ?, ?, ?)",
            (lin_id, solver, warm_start, accumulated_time, minlin.objbnd, gap),
        )
        cursor.execute("SELECT last_insert_rowid()")
        lin_id = cursor.fetchone()[0]

        # Commit the transaction
        conn.commit()

    else:
        lin_id = result[0]
        accumulated_time = result[4 + 1 + 4]
        lin = pickle.loads(result[4])

    return mlp, lin, lin_id, accumulated_time


def bestbound_linearization(cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver, warm_start):

    cursor.execute(
        """
        SELECT BaseLinearization.*, BBLinearization.*
            FROM BBLinearization
            JOIN BaseLinearization ON BBLinearization.id = BaseLinearization.id
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND BBLinearization.solver = ? AND BBLinearization.warm_start = ? AND BaseLinearization.linearization = 'bestbnd'
        """,
        (mlp_file, solver, warm_start),
    )
    result = cursor.fetchone()
    if result is None:

        min_triples = 0
        accumulated_time = 0

        # Just pick linearization size if warm start is not used
        if 0 == warm_start:
            mlp, lin, base_lin_id, accumulated_time = minimum_linearization(
                cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver, 0
            )
            # Size of the smallest linearization
            min_triples = len(lin)
            lin = None
            accumulated_time = 0

        # Pick linearization size and actual linearization
        else:
            mlp, lin, base_lin_id, accumulated_time = minimum_linearization(
                cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver, 1
            )
            # Size of the smallest linearization
            min_triples = len(lin)

        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        # Obtain BB linearization
        print("\t\tBest-Bound Linearization", mlp_file, solver, warm_start)
        start = timeit.timeit()
        bb_lin = optl.optimalLinearization(
            mlp, "bestbnd", max_triples=min_triples, time_limit=30, print_flag=0, init_lin=lin
        )
        end = timeit.timeit()

        # Replace linearization for bblin
        if bb_lin.nsolns > 0:
            lin = bb_lin.best_lins[0]

        # print("BB runtime",bb_lin.cpu)
        accumulated_time = bb_lin.cpu
        blob = pickle.dumps(lin)

        cursor.execute(
            "INSERT INTO BaseLinearization   (instance_id, linearization, ub,triples ) VALUES (?, ?, ?, ?)",
            (instance_id, "bestbnd", min_triples, blob),
        )
        cursor.execute("SELECT last_insert_rowid()")
        bb_id = cursor.fetchone()[0]

        gap = abs(bb_lin.obj - ceil(bb_lin.objbnd)) / abs(bb_lin.obj)

        cursor.execute(
            "INSERT INTO BBLinearization (id, solver, warm_start, base_lin_id, lin_size, runtime, primal, dual, gap) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (bb_id, solver, warm_start, base_lin_id, min_triples, accumulated_time, bb_lin.obj, bb_lin.objbnd, gap),
        )

        # Commit the transaction
        conn.commit()

    else:
        bb_id = result[0]
        lin = pickle.loads(result[4])
        accumulated_time = result[4 + 1 + 5]

    return mlp, lin, bb_id, accumulated_time


def ev_bestbound_linearization(cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver):

    warm_start = 0

    cursor.execute(
        """
        SELECT BaseLinearization.*, BBLinearization.*
            FROM BBLinearization
            JOIN BaseLinearization ON BBLinearization.id = BaseLinearization.id
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND BBLinearization.solver = ? AND BBLinearization.warm_start = ? AND BaseLinearization.linearization = 'bestbnd0'
        """,
        (mlp_file, solver, warm_start),
    )
    result = cursor.fetchone()
    if result is None:

        min_triples = 0
        accumulated_time = 0

        # Just pick linearization size if warm start is not used
        mlp, lin, base_lin_id, accumulated_time = ev_minimum_linearization(
            cursor, conn, mlp, instance_id, mlp_dir, mlp_file, solver
        )
        # Size of the smallest linearization
        min_triples = len(lin)
        lin = None
        accumulated_time = 0

        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        # Obtain BB linearization
        print("\t\tBest-Bound Linearization without warm start", mlp_file, solver, warm_start)
        start = timeit.timeit()
        bb_lin = optl.optimalLinearization(
            mlp, "bestbnd", max_triples=min_triples, time_limit=600, print_flag=0, init_lin=lin
        )
        end = timeit.timeit()

        # Replace linearization for bblin
        if bb_lin.nsolns > 0:
            lin = bb_lin.best_lins[0]

        # print("BB runtime",bb_lin.cpu)
        accumulated_time = bb_lin.cpu
        blob = pickle.dumps(lin)

        cursor.execute(
            "INSERT INTO BaseLinearization   (instance_id, linearization, ub,triples ) VALUES (?, ?, ?, ?)",
            (instance_id, "bestbnd0", min_triples, blob),
        )
        cursor.execute("SELECT last_insert_rowid()")
        bb_id = cursor.fetchone()[0]

        gap = abs(bb_lin.obj - ceil(bb_lin.objbnd)) / abs(bb_lin.obj)

        cursor.execute(
            "INSERT INTO BBLinearization (id, solver, warm_start, base_lin_id, lin_size, runtime, primal, dual, gap) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (bb_id, solver, warm_start, base_lin_id, min_triples, accumulated_time, bb_lin.obj, bb_lin.objbnd, gap),
        )

        # Commit the transaction
        conn.commit()

    else:
        bb_id = result[0]
        lin = pickle.loads(result[4])
        accumulated_time = result[4 + 1 + 5]

    return mlp, lin, bb_id, accumulated_time


def solve_lp(cursor, conn, lin, mlp, lin_id, mlp_dir, mlp_file, solver):

    # Check first if linearization has been computed
    cursor.execute("SELECT * FROM BaseLinearization WHERE id = ?", (lin_id,))
    entry = cursor.fetchone()
    if entry is None:
        print("Linearization unavailable. Exiting.")
        return mlp

    # LP relaxation
    cursor.execute(
        """
        SELECT LP.*
            FROM LP
            JOIN BaseLinearization ON LP.id = BaseLinearization.id
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND BaseLinearization.id = ?
        """,
        (mlp_file, lin_id),
    )
    result = cursor.fetchone()
    if result is None:

        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        print("\t\tLP Relaxation", mlp_file, entry[2], solver)

        lp_results = mlp.solveLP(lin=lin, time_limit=60)

        gap = abs((lp_results.obj - lp_results.objbnd) / abs(lp_results.obj))
        cursor.execute(
            "INSERT INTO LP  (id, solver, runtime, lb, ub, gap, obj) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (
                lin_id,
                solver,
                lp_results.cpu,
                min(lp_results.obj, lp_results.objbnd),
                max(lp_results.obj, lp_results.objbnd),
                gap,
                lp_results.obj,
            ),
        )
        conn.commit()

    return mlp


def solve_qcp_g11(
    cursor, conn, lin, mlp, lin_id, full_time, accumulated_time, mlp_dir, mlp_file, solver, function_call
):

    # Check first if linearization has been computed
    cursor.execute("SELECT * FROM BaseLinearization WHERE id = ?", (lin_id,))
    entry = cursor.fetchone()
    if entry is None:
        print("Linearization unavailable. Exiting.")
        return mlp

    # Solve QCP
    cursor.execute(
        """
        SELECT QCPG11.*, BaseLinearization.linearization
            FROM QCPG11
            JOIN BaseLinearization ON QCPG11.id = BaseLinearization.id
            JOIN Instance ON BaseLinearization.instance_id = Instance.id
            WHERE Instance.name = ? AND BaseLinearization.id = ?  AND QCPG11.solver = ?
        """,
        (mlp_file, lin_id, solver),
    )
    result = cursor.fetchone()
    if result is None:
        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file)

        print("\t\tQCPG11", mlp_file, entry[2], solver)
        results = function_call(mlp, full_time - accumulated_time, lin)

        gap = abs((results.ub - results.lb))
        objective = 0
        if mlp.obj_sense == "Min":
            gap = gap / abs(results.ub)
            objective = results.ub
        else:
            gap = gap / abs(results.lb)
            objective = results.lb
        cursor.execute(
            "INSERT OR REPLACE INTO QCPG11  (id, solver, runtime, lb, ub, gap, objective) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (lin_id, solver, results.runtime + accumulated_time, results.lb, results.ub, gap, objective),
        )
        conn.commit()

    return


def solve_po(cursor, conn, mlp, instance_id, full_time, mlp_dir, mlp_file, solver, function):
    # Solve PO without linearizations
    cursor.execute(
        """
        SELECT PO.*
            FROM PO
            JOIN Instance ON PO.id = Instance.id
            WHERE Instance.name = ? AND PO.solver = ?
        """,
        (mlp_file, solver),
    )
    result = cursor.fetchone()
    if result is None:
        if mlp is None:
            mlp = MLP.MLP(mlp_dir=mlp_dir, mlp_file=mlp_file, mlp_only=True)

        print("\t\tPO", mlp_file, solver)

        results = function(mlp, full_time, None)

        gap = abs((results.ub - results.lb))
        objective = 0
        if mlp.obj_sense == "Min":
            gap = gap / abs(results.ub)
            objective = results.ub
        else:
            gap = gap / abs(results.lb)
            objective = results.lb
        cursor.execute(
            "INSERT INTO PO  (id, solver, runtime, lb, ub, gap, objective) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (instance_id, solver, results.runtime, results.lb, results.ub, gap, objective),
        )
        conn.commit()

    return


def main(argv):

    full_time = 600
    mypath = "../data/"

    # Connect to the SQLite database
    conn = sqlite3.connect("../results/experiments.db")
    cursor = conn.cursor()

    solved_instances = dict()

    file1 = ""
    types = ["mult_d_3", "mult_d_4", "vision", "autocorr"]

    for type in types:

        onlyfiles = sorted([f for f in listdir(mypath + type) if isfile(join(mypath + type, f))])
        num_files = len(onlyfiles)
        my_dir = mypath + type + "/"

        count_files = 0
        for file in onlyfiles:
            count_files += 1
            print(mypath + type, file, count_files, "out of", num_files)

            # Get ID from instance
            mlp, instance_id = get_instance(cursor, conn, mypath + type + "/", file, type)

            # LINEARIZATIONS

            # Sequential linearization
            mlp, seq_lin, seq_id, _ = linearization(cursor, conn, mlp, instance_id, "sequential", my_dir, file)
            solve_lp(cursor, conn, seq_lin, mlp, seq_id, my_dir, file, "GUROBI11")
            solve_qcp_g11(cursor, conn, seq_lin, mlp, seq_id, full_time, 0, my_dir, file, "GUROBI11", wrapper_gurobi10)

            # Greedy linearization
            mlp, greedy_lin, greedy_id, _ = linearization(cursor, conn, mlp, instance_id, "greedy", my_dir, file)
            solve_lp(cursor, conn, greedy_lin, mlp, greedy_id, my_dir, file, "GUROBI11")
            solve_qcp_g11(
                cursor, conn, greedy_lin, mlp, greedy_id, full_time, 0, my_dir, file, "GUROBI11", wrapper_gurobi10
            )

            # Full linearization
            mlp, full_lin, full_id, _ = linearization(cursor, conn, mlp, instance_id, "full", my_dir, file)
            solve_lp(cursor, conn, full_lin, mlp, full_id, my_dir, file, "GUROBI11")
            solve_qcp_g11(
                cursor, conn, full_lin, mlp, full_id, full_time, 0, my_dir, file, "GUROBI11", wrapper_gurobi10
            )

            # Minimum linearization (with warm start)
            mlp, min_lin, min_lin_id, time_minlin = minimum_linearization(
                cursor, conn, mlp, instance_id, my_dir, file, "GUROBI11", 1
            )
            solve_lp(cursor, conn, min_lin, mlp, min_lin_id, my_dir, file, "GUROBI11")
            solve_qcp_g11(
                cursor,
                conn,
                min_lin,
                mlp,
                min_lin_id,
                full_time,
                time_minlin,
                my_dir,
                file,
                "GUROBI11",
                wrapper_gurobi10,
            )

            # Best-Bound linearization (with warm start)
            mlp, bb_lin, bb_id, time_bb = bestbound_linearization(
                cursor, conn, mlp, instance_id, my_dir, file, "GUROBI11", 1
            )
            solve_lp(cursor, conn, bb_lin, mlp, bb_id, my_dir, file, "GUROBI11")
            solve_qcp_g11(
                cursor,
                conn,
                bb_lin,
                mlp,
                bb_id,
                full_time,
                time_minlin + time_bb,
                my_dir,
                file,
                "GUROBI11",
                wrapper_gurobi10,
            )

            # EV COMPARISON
            # Minimum linearization (no warm start)
            mlp, _, _, _ = ev_minimum_linearization(cursor, conn, mlp, instance_id, my_dir, file, "GUROBI11")
            # Best-Bound linearization (without warm start)
            mlp, _, _, _ = ev_bestbound_linearization(cursor, conn, mlp, instance_id, my_dir, file, "GUROBI11")

            # PO without linearization (BARON+CPLEX)
            solve_po(cursor, conn, mlp, instance_id, full_time, my_dir, file, "BARON", wrapper_baron)

    # Close the cursor and the database connection
    cursor.close()
    conn.close()

    return


if __name__ == "__main__":
    main(sys.argv[1:])
