# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sqlite3

# Connect to a new SQLite database (creates a new file if it doesn't exist)
conn = sqlite3.connect("../results/experiments.db")

# Create a cursor object
cursor = conn.cursor()


# Create entry for each instance
cursor.execute(
    """CREATE TABLE IF NOT EXISTS Instance (
    id INTEGER PRIMARY KEY,
    name TEXT,
    type TEXT,
    n INTEGER,
    m INTEGER
)
"""
)


# Create entry for each instance
cursor.execute(
    """CREATE TABLE IF NOT EXISTS BaseLinearization (
    id INTEGER PRIMARY KEY,
    instance_id INTEGER,
    linearization STRING,
    ub INTEGER,
    triples BLOB,
    FOREIGN KEY (instance_id) REFERENCES Instance(id)
)
"""
)


# MIN-LINEARIZATION TABLE
cursor.execute(
    """
    CREATE TABLE IF NOT EXISTS MinLinearization (
        id INTEGER PRIMARY KEY,
        solver TEXT,
        warm_start INTEGER CHECK (warm_start IN (0,1)),
        base_lin_id INTEGER KEY,
        runtime REAL,
        lb INTEGER,
        gap REAL,
        FOREIGN KEY (id) REFERENCES BaseLinearization(id),
        FOREIGN KEY (base_lin_id) REFERENCES BaseLinearization(id)
    )
    """
)

# BB-LINEARIZATION TABLE
cursor.execute(
    """
    CREATE TABLE IF NOT EXISTS BBLinearization (
        id INTEGER PRIMARY KEY,
        solver TEXT,
        warm_start INTEGER CHECK (warm_start IN (0,1) ),
        base_lin_id INTEGER KEY,
        lin_size INTEGER,
        runtime REAL,
        primal REAL,
        dual REAL,
        gap REAL,
        FOREIGN KEY (id) REFERENCES BaseLinearization(id),
        FOREIGN KEY (base_lin_id) REFERENCES BaseLinearization(id)
    )
    """
)


# LP BOUNDS TABLES
cursor.execute(
    """
CREATE TABLE IF NOT EXISTS LP (
    id INTEGER PRIMARY KEY,
    solver TEXT,
    runtime REAL,
    lb INTEGER,
    ub INTEGER,
    gap REAL,
    obj REAL,
    FOREIGN KEY (id) REFERENCES BaseLinearization(id)
)
"""
)


# QCP TABLES
cursor.execute(
    """
CREATE TABLE IF NOT EXISTS QCPG11 (
    id INTEGER PRIMARY KEY,
    solver TEXT,
    runtime REAL,
    lb REAL,
    ub REAL,
    gap REAL,
    objective REAL,
    FOREIGN KEY (id) REFERENCES BaseLinearization(id)
)
"""
)


# PO TABLES
cursor.execute(
    """
CREATE TABLE IF NOT EXISTS PO (
    id INTEGER PRIMARY KEY,
    solver TEXT,
    runtime REAL,
    lb REAL,
    ub REAL,
    gap REAL,
    objective REAL,
    FOREIGN KEY (id) REFERENCES Instance(id)
)
"""
)


# Commit the transaction
conn.commit()


# Close the cursor and the database connection
cursor.close()
conn.close()
