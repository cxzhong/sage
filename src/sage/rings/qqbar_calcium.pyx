r"""
FLINT/Calcium backend for algebraic numbers

This module provides an optional engine for :mod:`sage.rings.qqbar` built on
FLINT's Calcium ``ca_t`` type. It is activated with
:func:`sage.rings.qqbar.set_algebraic_backend`; see there for details.
Nothing in this module is a stable public interface.

EXAMPLES::

    sage: from sage.rings.qqbar_calcium import ca_from_rational
    sage: ca_from_rational(QQ(1/3))
    Ca(0.3333333333333334?)

AUTHORS:

- Chenxin Zhong (2026-07-10): initial version (:issue:`42261`)
"""

# ****************************************************************************
#       Copyright (C) 2026 Chenxin Zhong <chenxin.zhong@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator

from cysignals.signals cimport sig_on, sig_off

from sage.libs.flint.types cimport (
    ca_t, ca_ctx_t, qqbar_t, qqbar_ptr, acb_t, truth_t,
    fmpz_t, fmpq_t, fmpz_poly_struct, slong,
    T_TRUE, T_FALSE, T_UNKNOWN)
from sage.libs.flint.ca cimport (
    ca_ctx_init, ca_ctx_clear, ca_init, ca_clear,
    ca_set_fmpq, ca_get_fmpq, ca_get_acb)
from sage.libs.flint.fmpq cimport fmpq_init, fmpq_clear, fmpq_set_mpq
from sage.libs.flint.acb cimport acb_init, acb_clear
from sage.rings.complex_arb cimport acb_to_ComplexIntervalFieldElement
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


cdef class CalciumContext:
    """
    Owner of the (unique, module-global) Calcium context object.

    TESTS::

        sage: from sage.rings.qqbar_calcium import CalciumContext
        sage: CalciumContext()
        <sage.rings.qqbar_calcium.CalciumContext object at ...>
    """
    cdef ca_ctx_t ctx

    def __cinit__(self):
        ca_ctx_init(self.ctx)

    def __dealloc__(self):
        ca_ctx_clear(self.ctx)


cdef CalciumContext _ctx = CalciumContext()


cdef inline object _truth(truth_t t):
    if t == T_TRUE:
        return True
    if t == T_FALSE:
        return False
    return None


cdef Ca _new_ca():
    return Ca.__new__(Ca)


cdef class Ca:
    """
    A wrapper around one Calcium ``ca_t`` value.

    Do not instantiate directly; use the ``ca_from_*`` module functions.

    TESTS::

        sage: from sage.rings.qqbar_calcium import ca_from_rational
        sage: a = ca_from_rational(2); a
        Ca(2)
    """
    cdef ca_t x
    cdef CalciumContext _ctxref

    def __cinit__(self):
        self._ctxref = _ctx
        ca_init(self.x, _ctx.ctx)

    def __dealloc__(self):
        ca_clear(self.x, self._ctxref.ctx)

    def __repr__(self):
        """
        TESTS::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: repr(ca_from_rational(QQ(-7/2)))
            'Ca(-3.5000000000000000?)'
        """
        return 'Ca({})'.format(self.enclosure(53))

    def enclosure(self, long prec=64):
        """
        Return a complex interval containing this value.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(QQ(1/3)).enclosure(53)
            0.3333333333333334?
            sage: ca_from_rational(QQ(1/3)).enclosure(200).diameter() < 2^-190
            True
        """
        cdef acb_t z
        cdef ComplexIntervalFieldElement res
        acb_init(z)
        sig_on()
        ca_get_acb(z, self.x, prec, _ctx.ctx)
        sig_off()
        res = ComplexIntervalField(prec)(0)
        acb_to_ComplexIntervalFieldElement(res, z)
        acb_clear(z)
        return res


def ca_from_rational(x):
    """
    Return a :class:`Ca` representing the rational number ``x``.

    INPUT:

    - ``x`` -- Integer, Rational, or Python int

    EXAMPLES::

        sage: from sage.rings.qqbar_calcium import ca_from_rational
        sage: ca_from_rational(QQ(2/3))
        Ca(0.6666666666666667?)
        sage: ca_from_rational(5)
        Ca(5)
    """
    cdef Rational q = <Rational> QQ(x)
    cdef Ca z = _new_ca()
    cdef fmpq_t t
    fmpq_init(t)
    fmpq_set_mpq(t, q.value)
    ca_set_fmpq(z.x, t, _ctx.ctx)
    fmpq_clear(t)
    return z
