# distutils: language = c++
# distutils: libraries = CBLAS_LIBRARIES
# distutils: library_dirs = CBLAS_LIBDIR
# distutils: include_dirs = CBLAS_INCDIR
r"""
Dense matrices over `\ZZ/n\ZZ` for `n < 2^{8}` using LinBox's ``Modular<float>``

AUTHORS:
- Burcin Erocal
- Martin Albrecht
"""
# #############################################################################
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.finite_rings.stdint cimport *

from sage.libs.linbox.givaro cimport \
    Modular_float as ModField, \
    Poly1Dom, Dense

from sage.libs.linbox.linbox cimport \
    DenseMatrix_Modular_float as DenseMatrix, \
    reducedRowEchelonize

from sage.libs.linbox.fflas cimport \
    fgemm, pfgemm, fgemv, Det, pDet, Rank, pRank, ReducedRowEchelonForm, pReducedRowEchelonForm, applyP, \
    MinPoly, CharPoly, \
    ModFloatDensePolynomial as ModDensePoly

ctypedef Poly1Dom[ModField, Dense] ModDensePolyRing

# LinBox supports up to 2^11 using float but that's double dog slow,
# so we pick a smaller value for crossover
MAX_MODULUS = 2**8

include "matrix_modn_dense_template.pxi"

cdef class Matrix_modn_dense_float(Matrix_modn_dense_template):
    r"""
    Dense matrices over `\ZZ/n\ZZ` for `n < 2^{8}` using LinBox's ``Modular<float>``.

    These are matrices with integer entries mod ``n`` represented as
    floating-point numbers in a 32-bit word for use with LinBox routines.
    This could allow for ``n`` up to `2^{11}`, but for performance reasons
    this is limited to ``n`` up to `2^{8}`, and for larger moduli the
    ``Matrix_modn_dense_double`` class is used.

    Routines here are for the most basic access, see the
    ``matrix_modn_dense_template.pxi`` file for higher-level routines.
    """
    def __cinit__(self):
        """
        The Cython constructor.

        TESTS::

            sage: A = random_matrix(GF(7), 4, 4, implementation='linbox')
            sage: type(A[0,0])
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
        """
        self._get_template = self._base_ring.zero()

    cdef void set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value) noexcept:
        r"""
        Set the (i,j) entry of ``self`` to the int value.

        EXAMPLES::

            sage: A = random_matrix(GF(7), 4, 4, implementation='linbox'); A
            [3 1 6 6]
            [4 4 2 2]
            [3 5 4 5]
            [6 2 2 1]
            sage: A[0,0] = 12; A
            [5 1 6 6]
            [4 4 2 2]
            [3 5 4 5]
            [6 2 2 1]

            sage: B = random_matrix(Integers(100), 4, 4, implementation='linbox'); B
            [13 95  1 16]
            [18 33  7 31]
            [92 19 18 93]
            [82 42 15 38]
            sage: B[0,0] = 422; B
            [22 95  1 16]
            [18 33  7 31]
            [92 19 18 93]
            [82 42 15 38]
        """
        self._matrix[i][j] = <float>value

    cdef unsigned long get_unsafe_ui(self, Py_ssize_t i, Py_ssize_t j):
        cdef float result = (<Matrix_modn_dense_template>self)._matrix[i][j]
        return <int_fast64_t>result

    cdef void set_unsafe_ui(self, Py_ssize_t i, Py_ssize_t j, unsigned long value):
        self._matrix[i][j] = <float>value

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        r"""
        Set the (i,j) entry with no bounds-checking, or any other checks.

        Assumes that `x` is in the base ring.

        EXAMPLES::

            sage: A = random_matrix(GF(13), 4, 4, implementation='linbox'); A
            [ 0  0  2  9]
            [10  6 11  8]
            [10 12  8  8]
            [ 3  6  8  0]
            sage: K = A.base_ring()
            sage: x = K(27)
            sage: A[0,0] = x
            sage: A[0,0] == x
            True
            sage: l[1:] == A.list()[1:]
            True

            sage: B = random_matrix(Integers(200), 4, 4, implementation='linbox'); B
            [ 13  95 101 116]
            [118 133   7 131]
            [192  19 118 193]
            [ 82 142 115  38]
            sage: R = B.base_ring()
            sage: x = R(311)
            sage: B[0,0] = x
            sage: B.list()[0] == x
            True
            sage: l[1:] == B.list()[1:]
            True
        """
        self._matrix[i][j] = <float>(<IntegerMod_int>x).ivalue

    cdef IntegerMod_int get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        r"""
        Return the (i,j) entry with no bounds-checking.

        OUTPUT:

        A :class:`sage.rings.finite_rings.integer_mod.IntegerMod_int`
        object.

        EXAMPLES::

            sage: A = random_matrix(Integers(100), 4, 4, implementation='linbox'); A
            [ 4 95 83 47]
            [44 57 91 53]
            [75 53 15 39]
            [26 25 10 74]
            sage: a = A[0,0]; a
            4
            sage: a in A.base_ring()
            True

            sage: B = random_matrix(Integers(100), 4, 4, implementation='linbox'); B
            [13 95  1 16]
            [18 33  7 31]
            [92 19 18 93]
            [82 42 15 38]
            sage: b = B[0,0]; b
            13
            sage: b in B.base_ring()
            True
        """
        cdef float result = (<Matrix_modn_dense_template>self)._matrix[i][j]
        return (<Matrix_modn_dense_float>self)._get_template._new_c(<int_fast32_t>result)

    cdef copy_from_unsafe(self, Py_ssize_t iDst, Py_ssize_t jDst, src, Py_ssize_t iSrc, Py_ssize_t jSrc):
        r"""
        Copy the ``(iSrc, jSrc)`` entry of ``src`` into the ``(iDst, jDst)``
        entry of ``self``.

        INPUT:

        - ``iDst`` - the row to be copied to in ``self``.
        - ``jDst`` - the column to be copied to in ``self``.
        - ``src`` - the matrix to copy from. Should be a Matrix_modn_dense_float
                    with the same base ring as ``self``.
        - ``iSrc``  - the row to be copied from in ``src``.
        - ``jSrc`` - the column to be copied from in ``src``.

        TESTS::

            sage: m = matrix(GF(131),3,4,range(12))
            sage: m
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            sage: m.transpose()
            [ 0  4  8]
            [ 1  5  9]
            [ 2  6 10]
            [ 3  7 11]
            sage: m.matrix_from_rows([0,2])
            [ 0  1  2  3]
            [ 8  9 10 11]
            sage: m.matrix_from_columns([1,3])
            [ 1  3]
            [ 5  7]
            [ 9 11]
            sage: m.matrix_from_rows_and_columns([1,2],[0,3])
            [ 4  7]
            [ 8 11]
        """
        cdef Matrix_modn_dense_float _src = <Matrix_modn_dense_float>src
        self._matrix[iDst][jDst] = _src._matrix[iSrc][jSrc]
