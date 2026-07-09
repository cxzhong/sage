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
    ca_set_fmpq, ca_get_fmpq, ca_get_acb,
    ca_add, ca_sub, ca_mul, ca_div, ca_neg, ca_inv, ca_sqrt, ca_abs, ca_re,
    ca_im, ca_conj, ca_pow_fmpq, ca_pow_si, ca_check_is_zero, ca_check_equal,
    ca_check_is_real, ca_check_lt)
from sage.libs.flint.fmpq cimport (
    fmpq_init, fmpq_clear, fmpq_set_mpq, fmpq_get_mpq)
from sage.libs.flint.acb cimport acb_init, acb_clear
from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si
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

    cdef inline int _guard_zero(self, str what) except -1:
        # raise if self is (possibly) zero; used before division-like ops
        z = _truth(ca_check_is_zero(self.x, _ctx.ctx))
        if z is True:
            raise ZeroDivisionError(what)
        if z is None:
            raise ValueError('Calcium cannot decide whether the operand is zero')
        return 0

    def add(self, Ca o):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(2).add(ca_from_rational(1/2))
            Ca(2.5000000000000000?)
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_add(z.x, self.x, o.x, _ctx.ctx)
        sig_off()
        return z

    def sub(self, Ca o):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(2).sub(ca_from_rational(3))
            Ca(-1)
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_sub(z.x, self.x, o.x, _ctx.ctx)
        sig_off()
        return z

    def mul(self, Ca o):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(3).mul(ca_from_rational(1/3))
            Ca(1)
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_mul(z.x, self.x, o.x, _ctx.ctx)
        sig_off()
        return z

    def div(self, Ca o):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(1).div(ca_from_rational(3))
            Ca(0.3333333333333334?)

        Division by zero raises rather than producing a Calcium
        special value::

            sage: ca_from_rational(1).div(ca_from_rational(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero
        """
        o._guard_zero('division by zero')
        cdef Ca z = _new_ca()
        sig_on()
        ca_div(z.x, self.x, o.x, _ctx.ctx)
        sig_off()
        return z

    def neg(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(2).neg()
            Ca(-2)
        """
        cdef Ca z = _new_ca()
        ca_neg(z.x, self.x, _ctx.ctx)
        return z

    def invert(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(4).invert()
            Ca(0.25000000000000000?)
            sage: ca_from_rational(0).invert()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inversion of zero
        """
        self._guard_zero('inversion of zero')
        cdef Ca z = _new_ca()
        sig_on()
        ca_inv(z.x, self.x, _ctx.ctx)
        sig_off()
        return z

    def sqrt(self):
        """
        Principal square root.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: s = ca_from_rational(2).sqrt(); s
            Ca(1.414213562373095?)
            sage: s.mul(s).equal(ca_from_rational(2))
            True
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_sqrt(z.x, self.x, _ctx.ctx)
        sig_off()
        return z

    def conjugate(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: i = ca_from_rational(-1).sqrt()
            sage: i.conjugate().add(i).is_zero()
            True
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_conj(z.x, self.x, _ctx.ctx)
        sig_off()
        return z

    def abs(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(-3).abs()
            Ca(3)
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_abs(z.x, self.x, _ctx.ctx)
        sig_off()
        return z

    def real(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: i = ca_from_rational(-1).sqrt()
            sage: i.real().is_zero()
            True
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_re(z.x, self.x, _ctx.ctx)
        sig_off()
        return z

    def imag(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: i = ca_from_rational(-1).sqrt()
            sage: i.imag().equal(ca_from_rational(1))
            True
        """
        cdef Ca z = _new_ca()
        sig_on()
        ca_im(z.x, self.x, _ctx.ctx)
        sig_off()
        return z

    def pow_rational(self, e):
        """
        ``self ** e`` for rational ``e``.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(2).pow_rational(QQ(1/2)).mul(
            ....:     ca_from_rational(2).pow_rational(QQ(1/2))).equal(ca_from_rational(2))
            True
            sage: ca_from_rational(0).pow_rational(QQ(-1))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: negative power of zero
        """
        cdef Rational q = <Rational> QQ(e)
        cdef Integer n = q.numerator()
        if q < 0:
            self._guard_zero('negative power of zero')
        cdef Ca z = _new_ca()
        cdef fmpq_t t
        if q.denominator() == 1 and mpz_fits_slong_p(n.value):
            sig_on()
            ca_pow_si(z.x, self.x, <slong> mpz_get_si(n.value), _ctx.ctx)
            sig_off()
            return z
        fmpq_init(t)
        fmpq_set_mpq(t, q.value)
        sig_on()
        ca_pow_fmpq(z.x, self.x, t, _ctx.ctx)
        sig_off()
        fmpq_clear(t)
        return z

    def _binop(self, Ca o, op):
        """
        Dispatch on an ``operator`` module function.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: import operator
            sage: ca_from_rational(2)._binop(ca_from_rational(3), operator.mul)
            Ca(6)
        """
        if op is operator.add:
            return self.add(o)
        if op is operator.sub:
            return self.sub(o)
        if op is operator.mul:
            return self.mul(o)
        if op is operator.truediv:
            return self.div(o)
        raise ValueError('unsupported operator')

    def is_zero(self):
        """
        Return ``True``/``False``, or ``None`` if Calcium cannot decide.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(0).is_zero()
            True
            sage: s = ca_from_rational(2).sqrt()
            sage: s.mul(s).sub(ca_from_rational(2)).is_zero()
            True
        """
        return _truth(ca_check_is_zero(self.x, _ctx.ctx))

    def equal(self, Ca o):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: a = ca_from_rational(2).sqrt().add(ca_from_rational(3).sqrt())
            sage: b = ca_from_rational(3).sqrt().add(ca_from_rational(2).sqrt())
            sage: a.equal(b)
            True
            sage: a.equal(ca_from_rational(2))
            False
        """
        return _truth(ca_check_equal(self.x, o.x, _ctx.ctx))

    def is_real(self):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(2).sqrt().is_real()
            True
            sage: ca_from_rational(-2).sqrt().is_real()
            False
        """
        return _truth(ca_check_is_real(self.x, _ctx.ctx))

    def sign_real(self):
        """
        Sign of the real part: `-1`, `0`, `1`, or ``None`` if undecided.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(-5).sign_real()
            -1
            sage: ca_from_rational(-1).sqrt().sign_real()
            0
        """
        cdef Ca r = self.real()
        return r._sign_of_real_value()

    def sign_imag(self):
        """
        Sign of the imaginary part: `-1`, `0`, `1`, or ``None``.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: ca_from_rational(-1).sqrt().sign_imag()
            1
            sage: ca_from_rational(3).sign_imag()
            0
        """
        cdef Ca r = self.imag()
        return r._sign_of_real_value()

    def _sign_of_real_value(self):
        # sign of self, which the caller promises is real-valued
        z = _truth(ca_check_is_zero(self.x, _ctx.ctx))
        if z is True:
            return 0
        if z is None:
            return None
        cdef Ca zero = ca_from_rational(0)
        lt = _truth(ca_check_lt(self.x, zero.x, _ctx.ctx))
        if lt is True:
            return -1
        if lt is False:
            return 1
        return None

    def cmp_lex(self, Ca o):
        """
        Compare lexicographically by (real, imaginary) parts, QQbar's
        ordering. Return `-1`, `0`, `1`, or ``None`` if undecided.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: i = ca_from_rational(-1).sqrt()
            sage: i.cmp_lex(ca_from_rational(1))
            -1
            sage: i.cmp_lex(i.neg())
            1
            sage: ca_from_rational(2).cmp_lex(ca_from_rational(2))
            0
        """
        cdef Ca d = self.sub(o)
        s = d.sign_real()
        if s is None:
            return None
        if s != 0:
            return s
        return d.sign_imag()

    def get_rational(self):
        """
        Return this value as a Rational, or ``None`` if it is not rational.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: s = ca_from_rational(2).sqrt()
            sage: s.mul(s).get_rational()
            2
            sage: s.get_rational() is None
            True
        """
        cdef fmpq_t t
        cdef Rational res
        fmpq_init(t)
        sig_on()
        ok = ca_get_fmpq(t, self.x, _ctx.ctx)
        sig_off()
        if not ok:
            fmpq_clear(t)
            return None
        res = Rational.__new__(Rational)
        fmpq_get_mpq(res.value, t)
        fmpq_clear(t)
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
