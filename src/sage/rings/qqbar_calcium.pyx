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
    ca_set_fmpq, ca_set_qqbar, ca_get_fmpq, ca_get_qqbar, ca_get_acb,
    ca_add, ca_sub, ca_mul, ca_div, ca_neg, ca_inv, ca_sqrt, ca_abs, ca_re,
    ca_im, ca_conj, ca_pow_fmpq, ca_pow_si, ca_check_is_zero, ca_check_equal,
    ca_check_is_real, ca_check_lt)
from sage.libs.flint.fmpq cimport (
    fmpq_init, fmpq_clear, fmpq_set_mpq, fmpq_get_mpq)
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_get_mpz
from sage.libs.flint.qqbar cimport (
    qqbar_init, qqbar_clear, _qqbar_vec_init, _qqbar_vec_clear,
    qqbar_get_acb, qqbar_roots_fmpz_poly, qqbar_equal)
from sage.libs.flint.acb cimport acb_init, acb_clear
from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si
from sage.rings.complex_arb cimport acb_to_ComplexIntervalFieldElement
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.rings.rational cimport Rational

from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


cdef extern from "flint/qqbar.h":
    const fmpz_poly_struct* QQBAR_POLY(const qqbar_t x)

cdef extern from "flint/fmpz_poly.h":
    slong _minpoly_degree "fmpz_poly_degree"(const fmpz_poly_struct* poly)
    void _minpoly_coeff "fmpz_poly_get_coeff_fmpz"(fmpz_t x, const fmpz_poly_struct* poly, slong n)


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

    def to_qqbar_data(self, long prec=64):
        """
        Return ``(coeffs, enclosure)`` where ``coeffs`` are the integer
        coefficients (constant first) of the minimal polynomial of this
        (necessarily algebraic) value and ``enclosure`` is a complex
        interval isolating it among the roots.

        EXAMPLES::

            sage: from sage.rings.qqbar_calcium import ca_from_rational
            sage: s = ca_from_rational(2).sqrt()
            sage: coeffs, e = s.to_qqbar_data()
            sage: coeffs
            [-2, 0, 1]
            sage: e.real() in RIF(1.414, 1.415)
            True
        """
        cdef qqbar_t q
        cdef acb_t z
        cdef fmpz_t c
        cdef ComplexIntervalFieldElement enc
        cdef Integer zc
        cdef const fmpz_poly_struct* p
        cdef slong i, d
        qqbar_init(q)
        try:
            sig_on()
            ok = ca_get_qqbar(q, self.x, _ctx.ctx)
            sig_off()
            if not ok:
                raise ValueError('unable to express this value as an algebraic number')
            p = QQBAR_POLY(q)
            d = _minpoly_degree(p)
            coeffs = []
            fmpz_init(c)
            try:
                for i in range(d + 1):
                    _minpoly_coeff(c, p, i)
                    zc = Integer.__new__(Integer)
                    fmpz_get_mpz(zc.value, c)
                    coeffs.append(zc)
            finally:
                fmpz_clear(c)
            acb_init(z)
            try:
                sig_on()
                qqbar_get_acb(z, q, prec)
                sig_off()
                enc = ComplexIntervalField(prec)(0)
                acb_to_ComplexIntervalFieldElement(enc, z)
            finally:
                acb_clear(z)
            return coeffs, enc
        finally:
            qqbar_clear(q)

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


def ca_from_minpoly_root(coeffs, enclosure):
    """
    Return the :class:`Ca` root of the integer polynomial with coefficient
    list ``coeffs`` (constant first; a ZZ or QQ polynomial is also accepted)
    isolated by the complex interval ``enclosure``.

    EXAMPLES::

        sage: from sage.rings.qqbar_calcium import ca_from_minpoly_root, ca_from_rational
        sage: r = ca_from_minpoly_root([-2, 0, 1], CIF(1.4142135623730951))
        sage: r.equal(ca_from_rational(2).sqrt())
        True

    A non-isolating enclosure is rejected::

        sage: ca_from_minpoly_root([-2, 0, 1], CIF(RIF(-2, 2)))
        Traceback (most recent call last):
        ...
        ValueError: enclosure does not isolate a single root
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.polynomial.polynomial_element import Polynomial
    R = PolynomialRing(ZZ, 'y')
    if isinstance(coeffs, Polynomial):
        p = coeffs
        p = p * p.denominator() if p.base_ring() is QQ else p
        zp_sage = R(p)
    else:
        zp_sage = R([ZZ(c) for c in coeffs])
    if zp_sage.degree() <= 0:
        raise ValueError('polynomial must be nonconstant')
    # squarefree part, so every root appears exactly once below
    zp_sage = zp_sage // zp_sage.gcd(zp_sage.derivative())
    cdef Polynomial_integer_dense_flint zp = <Polynomial_integer_dense_flint> zp_sage
    cdef slong d = zp_sage.degree()
    cdef qqbar_ptr vec
    cdef acb_t z
    cdef ComplexIntervalFieldElement re
    cdef slong i
    cdef long prec = 64
    cdef Ca res
    from sage.rings.complex_interval_field import ComplexIntervalField as _CIF
    vec = _qqbar_vec_init(d)
    try:
        sig_on()
        qqbar_roots_fmpz_poly(vec, zp._poly, 0)
        sig_off()
        target = _CIF(53)(enclosure)
        while prec <= (1 << 20):
            hits = []
            for i in range(d):
                acb_init(z)
                try:
                    qqbar_get_acb(z, &vec[i], prec)
                    re = _CIF(prec)(0)
                    acb_to_ComplexIntervalFieldElement(re, z)
                finally:
                    acb_clear(z)
                if re.overlaps(_CIF(prec)(target)):
                    hits.append(i)
            if len(hits) == 1:
                res = _new_ca()
                sig_on()
                ca_set_qqbar(res.x, &vec[hits[0]], _ctx.ctx)
                sig_off()
                return res
            if not hits:
                raise ValueError('enclosure does not isolate a single root')
            prec *= 2
        raise ValueError('enclosure does not isolate a single root')
    finally:
        _qqbar_vec_clear(vec, d)


class ANCalciumBase:
    """
    Placeholder base marker for :class:`ANCalcium` (defined after the
    descriptor machinery); only used for ``isinstance`` in
    :func:`ca_from_descr`.

    TESTS::

        sage: from sage.rings.qqbar_calcium import ANCalciumBase
        sage: ANCalciumBase
        <class 'sage.rings.qqbar_calcium.ANCalciumBase'>
    """
    pass


def ca_from_descr(descr):
    """
    Convert a ``qqbar.py`` descriptor to a :class:`Ca`, or return ``None``
    when the descriptor is of a kind this backend cannot handle (the caller
    then falls back to the classical path).

    EXAMPLES::

        sage: from sage.rings.qqbar_calcium import ca_from_descr, ca_from_rational
        sage: ca_from_descr(QQbar(7/3)._descr).get_rational()
        7/3
        sage: rt2 = QQbar(2).sqrt()
        sage: ca_from_descr(rt2._descr).mul(ca_from_descr(rt2._descr)).get_rational()
        2
        sage: s = QQbar(2).sqrt() + QQbar(3).sqrt()
        sage: ca_from_descr(s._descr).equal(
        ....:     ca_from_rational(2).sqrt().add(ca_from_rational(3).sqrt()))
        True
        sage: rt2.exactify()
        sage: type(rt2._descr).__name__
        'ANExtensionElement'
        sage: ca_from_descr(rt2._descr).mul(ca_from_descr(rt2._descr)).get_rational()
        2

    Roots of polynomials with non-rational algebraic coefficients are
    unsupported (the caller falls back to the classical path)::

        sage: x = polygen(QQbar)
        sage: r = (x^2 - QQbar(2).sqrt()).roots(QQbar, multiplicities=False)[0]
        sage: ca_from_descr(r._descr) is None
        True
    """
    from sage.rings.qqbar import (ANRational, ANRoot, ANExtensionElement,
                                  ANUnaryExpr, ANBinaryExpr)
    if isinstance(descr, ANCalciumBase):
        return descr._ca
    if isinstance(descr, ANRational):
        return ca_from_rational(descr._value)
    if isinstance(descr, ANRoot):
        cached = getattr(descr, '_calcium_cache', None)
        if cached is not None:
            return cached
        p = descr._poly._poly
        if p.base_ring() is not QQ:
            # Rational powering constructs its polynomial over QQbar even
            # when all coefficients are rational.
            try:
                p = p.change_ring(QQ)
            except (TypeError, ValueError):
                return None
        try:
            res = ca_from_minpoly_root(p, descr._interval)
        except ValueError:
            return None
        descr._calcium_cache = res
        return res
    if isinstance(descr, ANExtensionElement):
        gen = descr._generator
        ca_gen = getattr(gen, '_calcium_cache', None)
        if ca_gen is None:
            ca_gen = ca_from_descr(gen._root)
            if ca_gen is None:
                return None
            gen._calcium_cache = ca_gen
        res = ca_from_rational(0)
        for c in reversed(descr._value.polynomial().list()):
            res = res.mul(ca_gen).add(ca_from_rational(c))
        return res
    if isinstance(descr, ANUnaryExpr):
        arg = ca_from_descr(descr._arg._descr)
        if arg is None:
            return None
        op = descr._op
        try:
            if op == '-':
                return arg.neg()
            if op == '~':
                return arg.invert()
            if op == 'conjugate':
                return arg.conjugate()
            if op == 'abs':
                return arg.abs()
            if op == 'real':
                return arg.real()
            if op == 'imag':
                return arg.imag()
            if op == 'norm':
                return arg.mul(arg.conjugate())
        except (ZeroDivisionError, ValueError):
            return None
        return None
    if isinstance(descr, ANBinaryExpr):
        left = ca_from_descr(descr._left._descr)
        if left is None:
            return None
        right = ca_from_descr(descr._right._descr)
        if right is None:
            return None
        try:
            return left._binop(right, descr._op)
        except (ZeroDivisionError, ValueError):
            return None
    return None
