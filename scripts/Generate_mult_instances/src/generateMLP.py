# Copyright (c) 2022-2024 Mitsubishi Electric Research Laboratories (MERL).
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import random as rnd

import numpy as np


class GenerateMLP(object):
    """Constructor for generating MLP
    n - number of variables
    m - number of monomials
    d - degree of the monomials
    ni - number of instances
    """

    def __init__(self, n, m, d, ni):

        self.n = n
        self.m = m
        self.d = d
        self.ni = ni

        rnd.seed(rnd.sample(range(1000), 1)[0])

    def generate_all(self, fpath):

        for i in range(self.ni):
            fname = (
                fpath
                + "/mult_n_"
                + str(self.n)
                + "_m_"
                + str(self.m)
                + "_d_"
                + str(self.d)
                + "_s_"
                + str(i + 1)
                + ".dat"
            )
            monomials, coeff = self.generate()
            self.write_mlp(fname, monomials, coeff)

    def generate_all_const(self, fpath):

        for i in range(self.ni):

            for k in [0, self.n // 5, self.n // 2, self.n]:
                fname = (
                    fpath
                    + "/mult_const_"
                    + str(self.n)
                    + "_m_"
                    + str(self.m)
                    + "_d_"
                    + str(self.d)
                    + "_s_"
                    + str(i + 1)
                    + "_"
                    + str(k)
                    + ".dat"
                )

                monomials, coeff = self.generate()

                list_constraints = list()
                for _ in range(k):
                    mon = list(zip([x for x in range(1, self.n + 1) if rnd.choice((True, False))]))
                    c = [1 for x in mon]
                    b = 0.5 + 5 * len(mon) // 6

                    list_constraints.append((mon, c, b))
                self.write_mlp_const(fname, monomials, coeff, list_constraints)

    def generate(self):

        monomials = []
        m = 0
        n = 0
        var_list = list(range(1, self.n + 1))

        while m < self.m:
            # we have to ensure all n variables are in at least one monomial
            if len(var_list) > 0:
                if len(var_list) >= self.d:
                    tmp_list = var_list
                else:
                    tmp_list = var_list
                    for i in range(1, self.n + 1):
                        if not (i in var_list):
                            tmp_list.append(i)
                        if len(tmp_list) == self.d:
                            break
            else:
                tmp_list = range(1, self.n + 1)

            tmp_mon = rnd.sample(tmp_list, self.d)
            tmp_mon.sort()
            tmp_mon_list = tuple(tmp_mon)
            if not (tmp_mon_list in monomials):
                monomials.append(tmp_mon_list)
                m = m + 1

            if len(var_list) > 0:
                if len(var_list) <= self.d:
                    var_list = []
                else:
                    for i in tmp_mon_list:
                        var_list.remove(i)

        coeff = [rnd.uniform(-1.0, 1.0) for i in range(self.m)]

        return monomials, coeff

    def write_mlp(self, fname, monomials, coeff):

        fp = open(fname, "w")
        fp.write("#Variables " + str(self.n) + "\n")
        fp.write("#Constraints 0\n")
        fp.write("Objsense Min\n")
        fp.write("VariablesInfo\n")
        for i in range(self.n):
            fp.write("0.0 1.0 Cont\n")
        fp.write("Objective " + str(self.m) + "\n")
        fp.write("Offset 0.0\n")
        for i in range(self.m):
            mon = monomials[i]
            fp.write("[")
            for j in range(self.d):
                fp.write(str(mon[j]))
                if j < self.d - 1:
                    fp.write(", ")
                else:
                    fp.write("] " + str(round(coeff[i], 2)) + "\n")
        fp.close()

    def write_mlp_const(self, fname, monomials, coeff, list_constraints):

        fp = open(fname, "w")
        fp.write("#Variables " + str(self.n) + "\n")
        fp.write("#Constraints " + str(len(list_constraints)) + "\n")
        fp.write("Objsense Min\n")
        fp.write("VariablesInfo\n")
        for i in range(self.n):
            fp.write("0.0 1.0 Bin\n")
        fp.write("Objective " + str(self.m) + "\n")
        fp.write("Offset 0.0\n")
        for i in range(self.m):
            mon = monomials[i]
            fp.write("[")
            for j in range(self.d):
                fp.write(str(mon[j]))
                if j < self.d - 1:
                    fp.write(", ")
                else:
                    fp.write("] " + str(round(coeff[i], 2)) + "\n")
        for c in range(len(list_constraints)):
            mon = list_constraints[c][0]
            coeff = list_constraints[c][1]
            b = list_constraints[c][2]
            fp.write("Constraint" + str(c + 1) + " " + str(len(mon)) + "\n")
            fp.write("UB " + str(b) + "\n")
            # print("Constraint",c,mon,coeff,len(mon))
            for i in range(len(mon)):
                # mon = monomials[i]
                fp.write("[")
                for j in range(len(mon[i])):
                    fp.write(str(mon[i][j]))
                    if j < len(mon[i]) - 1:
                        fp.write(", ")
                    else:
                        fp.write("] " + str(coeff[i]) + "\n")
        fp.close()
