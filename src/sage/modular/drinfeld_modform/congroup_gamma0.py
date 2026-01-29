r"""
Congruence subgroup Gamma0(N)
"""

# ****************************************************************************
#       Copyright (C) 2026 Cécile Armana <armana@math.cnrs.fr>
#                          Xavier Caruso <xavier@caruso.ovh>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.lazy_attribute import lazy_attribute
from sage.groups.group import Group
from sage.structure.element import MultiplicativeGroupElement
from sage.categories.finite_fields import FiniteFields
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_element import Polynomial


class Gamma0Element(MultiplicativeGroupElement):
    def __init__(self, parent, mat):
        MultiplicativeGroupElement.__init__(self, parent)
        MS = parent.matrix_space()
        mat = self._mat = MS(mat)
        level = parent.level()
        if not mat.is_invertible():
            raise ValueError("not in Gamma0(%s)" % level)
        if mat[1,0] % level != 0:
            raise ValueError("not in Gamma0(%s)" % level)

    def __repr__(self):
        return self._mat.__repr__()

    def level(self):
        return self.parent().level()

    def _mul_(self, other):
        return self.parent()(self._mat * other._mat)

    def __invert__(self):
        return self.parent()(self._mat.inverse())


class Gamma0_class(Group):
    Element = Gamma0Element

    def __init__(self, level):
        Group.__init__(self)
        self._base = A = level.parent()
        if not (isinstance(level, Polynomial) and A.base_ring() in FiniteFields()):
            raise TypeError("base ring must be a polynomial ring over a finite field")
        self._level = level
        self._matrix_space = MatrixSpace(A, 2)
        self._q = A.base_ring().cardinality()

    def _repr_(self):
        return "Congruence Subgroup Gamma0(%s)" % self._level

    def base_ring(self):
        return self._base

    def level(self):
        return self._level

    @lazy_attribute
    def _level_factorized(self):
        return self._level.factor()

    def _an_element_(self):
        N = self.level()
        return self([1-N, -N, N, 1+N])

    def matrix_space(self):
        return self._matrix_space

    def genus(self):
        q = self._q
        L = self._level_factorized
        s = len(L)
        epsilon = kappa = r = 1
        for P, mult in L:
            d = P.degree()
            epsilon *= q**(d*(mult-1)) * (1 + q**d)
            kappa *= q**(d * (mult // 2)) + q**(d * ((mult-1) // 2))
            if d % 2 == 1:
                r = 0
        g = 1 + (epsilon - (q+1)*kappa - 2**(s-1) * (r*q*(q-1) + (q+1)*(q-2))) / (q**2 - 1)
        return g

    def ncusps(self):
        q = self._q
        L = self._level_factorized
        s = len(L)
        kappa = 1
        for P, r in L:
            d = P.degree()
            kappa *= q**(d * (r // 2)) + q**(d * ((r-1) // 2))
        h = 2**s + (kappa - 2**s) / (q - 1)
        return h
