# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import random
import sys
from os import listdir
from os.path import isfile, join

path = os.path.relpath("./src")
sys.path.append(path)

import src.generateMLP as gmlp

n_variables = [40]
n_monomials = [100]
degrees = [3, 4]
n_instances = 3
random.seed(42)


for n_var in n_variables:
    for n_mon in n_monomials:
        for d in degrees:
            gen = gmlp.GenerateMLP(n_var, n_mon, d, n_instances)
            gen.generate_all("./data")
