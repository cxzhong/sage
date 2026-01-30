r"""
Congruence subgroup '\Gamma0(N)' for polynomial rings over finite fields
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
    r"""
    An element in a congruence subgroup.
    """
    def __init__(self, parent, mat):
        r"""
        Initialize this element of the congruence subgroup.

        TESTS:

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4+2*T+3)
            sage: g = G.an_element()
            sage: TestSuite(g).run()
        """
        MultiplicativeGroupElement.__init__(self, parent)
        MS = parent.matrix_space()
        mat = self._mat = MS(mat)
        level = parent.level()
        if not mat.is_invertible():
            raise ValueError("not in Gamma0(%s)" % level)
        if mat[1,0] % level != 0:
            raise ValueError("not in Gamma0(%s)" % level)

    def __repr__(self):
        r"""
        Return the string representation of the element ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4+2*T+3)
            sage: g = G.an_element()
            sage: g.__repr__()
            '[4*T^4 + 3*T + 3 4*T^4 + 3*T + 2]\n[  T^4 + 2*T + 3   T^4 + 2*T + 4]'
        """
        return self._mat.__repr__()

    def level(self):
        r"""
        Return the level of the congruence subgroup in which the element 
        ``self`` lives.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4+2*T+3)
            sage: g = G.an_element()
            sage: g.level()
            T^4 + 2*T + 3
        """
        return self.parent().level()

    def _mul_(self, other):
        r"""
        Return the product of the element ``self`` with ``other``.

        EXAMPLES::

            sage: A.<T> = GF(5)[], N = T^4+2*T+3
            sage: G = Gamma0(N)
            sage: g = G.an_element()
            sage: h = G(matrix([[1,0],[N,1]])
            sage: g._mul_(h)
        
        """
        return self.parent()(self._mat * other._mat)

    def __invert__(self):
        r"""
        Return the inverse of the element ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4+2*T+3)
            sage: g = G.an_element()
            sage: g.__invert__()
        """
        return self.parent()(self._mat.inverse())


class Gamma0_class(Group):
    r"""
    A congruence subgroup.
    """
    Element = Gamma0Element

    def __init__(self, level):
        r"""
        The congruence subgroup `\Gamma_0(N)`.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4+2*T+3)
            sage: G
            Congruence Subgroup Gamma0(T^4 + 2*T + 3)
            sage: TestSuite(G).run()
        """
        Group.__init__(self)
        self._base = A = level.parent()
        if not (isinstance(level, Polynomial) and A.base_ring() in FiniteFields()):
            raise TypeError("Base ring must be a polynomial ring over a finite field")
        self._level = level
        self._matrix_space = MatrixSpace(A, 2)
        self._q = A.base_ring().cardinality()

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3)._repr_()
            'Congruence Subgroup Gamma0(T^4 + 2*T + 3)'        
        """
        return "Congruence Subgroup Gamma0(%s)" % self._level
        
    def base_ring(self):
        r"""
        Return the base ring of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3).base_ring()
            Univariate Polynomial Ring in T over Finite Field of size 5
        """
        return self._base

    def level(self):
        r"""
        Return the level of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3).level()
            T^4 + 2*T + 3
        """
        return self._level

    @lazy_attribute
    def _level_factorized(self):
        r"""
        Return the factorization of the level of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3)._level_factorized()
            (T + 2)^2 * (T^2 + T + 2)
        """
        return self._level.factor()
        
    @lazy_attribute
    def _level_factorizedradical(self):
        r"""
        Return the factorized radical of the level of ``self`` i.e. the product 
        of the distinct irreducible factors of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3)._level_factorizedradical()
            (T + 2) * (T^2 + T + 2)
        """
        return self._level.radical().factor()

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3).an_element()
            [4*T^4 + 3*T + 3 4*T^4 + 3*T + 2]
            [  T^4 + 2*T + 3   T^4 + 2*T + 4]
        """
        N = self.level()
        return self([1-N, -N, N, 1+N])

    def matrix_space(self):
        r"""
        Return the matrix space of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4+2*T+3).matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial 
            Ring in T over Finite Field of size 5 
        """
        return self._matrix_space

    def genus(self):
        r""""
        Return the genus of ``self``, i.e. the genus of the Drinfeld modular 
        curve attached to ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T-1).genus()
            0
            sage: Gamma0(T^3+1).genus()
            5
            sage: Gamma0(T^5-3*T^4+1).genus()
            155
        
        If N is irreducible, the formula simplifies as follows:

            sage: def genus_irr(N):
            ....:     q = N.base_ring().cardinality()
            ....:     d = N.degree()
            ....:     if is_even(d):
            ....:          g = (q^d-q^2)/(q^2-1)
            ....:     else:
            ....:          g = (q^d-q)/(q^2-1)
            ....:     return g
            sage: N = T^6 + T^4 + 4*T^3 + T^2 + 2; N.is_irreducible()
            True
            sage: Gamma0(N).genus() == genus_irr(N)
            True
            sage: N = T^7+3*T+3; N.is_irreducible()
            True
            sage: Gamma0(N).genus() == genus_irr(N)
            True 

        REFERENCE:

        [Gek2001]_
        """
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
        r"""
        Return the number of cusps of ``self`` i.e. the number of cusps of the 
        Drinfeld modular curve attached to ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T-1).ncusps()
            2
            sage: Gamma0(T^3+1).ncusps()
            4
            sage: Gamma0(T^7-3*T^4+1).ncusps()
            8

        If N is irreducible, the number of cusps is 2. We can check it:

            sage: N = T^6 + T^4 + 4*T^3 + T^2 + 2; N.is_irreducible()
            True
            sage: Gamma0(N).ncusps() == 2
            True
            sage: N = T^7+3*T+3; N.is_irreducible()
            True
            sage: Gamma0(N).ncusps() == 2
            True 

        REFERENCE:

        [Gek2001]_
        """
        q = self._q
        L = self._level_factorized
        s = len(L)
        kappa = 1
        for P, r in L:
            d = P.degree()
            kappa *= q**(d * (r // 2)) + q**(d * ((r-1) // 2))
        h = 2**s + (kappa - 2**s) / (q - 1)
        return h
        
    def index(self):
        r"""
        Return the index of self as a subgroup of `\mathrm{GL}_{2}(A)` where `A`
        is the base ring of ``self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4-3*T^4+1).index()
            626

        If the level `N` of ``self``, is irreducible, the index is `1+q^\deg(N)`.
        We can check it:

            sage: A.<T> = GF(5)[]
            sage: N = T^4 + 4*T^2 + 4*T + 2; N.is_irreducible()
            True
            sage: q = A.base_ring().cardinality()
            sage: Gamma0(N).index() == 1+q**N.degree()
            True
        """
        q = self._q
        level = self._level
        L = list(self._level_factorizedradical)
        ind = q**level.degree()
        for P,r in L:
            ind *= 1+1/(q**P.degree())
        return ind     
        
