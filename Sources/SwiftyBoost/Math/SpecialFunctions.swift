//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//  

import CBoostBridge

// Errors for invalid inputs to special functions
public enum SpecialFunctionError: Error {
    case parameterNotPositive(name: String)
    case parameterOutOfRange(name: String, min: Double, max: Double)
    case parameterNotFinite(name: String)
    case poleAtNonPositiveInteger(name: String)
    case invalidCombination(message: String)
}

// Helper casts
@usableFromInline internal func D<T: BinaryFloatingPoint>(_ x: T) -> Double { Double(x) }

// Boost constants (via CBoostBridge)
@inlinable public var boostE: Double { bs_const_e() }
// If you need π later: @inlinable public var boostPi: Double { bs_const_pi() }

// MARK: - Gamma / Error

@inlinable public func gamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    // Poles at non-positive integers (0, -1, -2, ...)
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_tgamma(dx))
}

@inlinable public func logGamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    // Poles at non-positive integers (0, -1, -2, ...)
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_lgamma(dx))
}

@inlinable public func errorFunction<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    // erf is defined for all real x; guard against non-finite input
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_erf(dx))
}

@inlinable public func complementaryErrorFunction<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    // erfc is defined for all real x; guard against non-finite input
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_erfc(dx))
}

// MARK: - Beta family

/// Euler’s beta function B(a, b).
@inlinable public func beta<T: BinaryFloatingPoint>(_ a: T, _ b: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    return T(bs_beta(D(a), D(b)))
}

/// Unnormalized incomplete beta B(a, b; x).
@inlinable public func incompleteBetaUnnormalized<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return T(bs_fullBeta(D(a), D(b), D(x)))
}

/// Regularized incomplete beta I_x(a, b).
@inlinable public func regularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return T(bs_ibeta(D(a), D(b), D(x)))
}

/// Complement of regularized incomplete beta I_{1−x}(b, a).
@inlinable public func complementaryRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return T(bs_ibetac(D(a), D(b), D(x)))
}

/// Solve for x in I_x(a, b) = p.
@inlinable public func inverseRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, p: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return T(bs_ibeta_inv(D(a), D(b), D(p)))
}

/// Solve for x in I_{1−x}(a, b) = p (complement).
@inlinable public func inverseComplementaryRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, p: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return T(bs_ibetac_inv(D(a), D(b), D(p)))
}

/// Solve for parameter a in I_x(a, b) = p, given b, x, p.
@inlinable public func solveAForRegularizedIncompleteBeta<T: BinaryFloatingPoint>(b: T, x: T, p: T) -> T {
    T(bs_ibeta_inva(D(b), D(x), D(p)))
}

/// Solve for parameter b in I_x(a, b) = p, given a, x, p.
@inlinable public func solveBForRegularizedIncompleteBeta<T: BinaryFloatingPoint>(a: T, x: T, p: T) -> T {
    T(bs_ibeta_invb(D(a), D(x), D(p)))
}

// MARK: - Digamma / Polygamma / Zeta

@inlinable public func digamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    // Domain: poles at non-positive integers (0, -1, -2, ...)
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // Check if dx is an integer and non-positive
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_digamma(dx))
}

@inlinable public func trigamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    // Same poles as digamma at non-positive integers
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_trigamma(dx))
}

/// Polygamma of order n at x.
@inlinable public func polygamma<T: BinaryFloatingPoint>(order n: Int, _ x: T) throws -> T {
    // Order must be >= 0 (n = 0 corresponds to digamma); same poles in x as digamma
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_polygamma(Int32(n), dx))
}

/// Riemann zeta ζ(x).
@inlinable public func riemannZeta<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // Simple pole at x = 1
    guard dx != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
    return T(bs_riemann_zeta(dx))
}

// MARK: - Owen's T

/// Owen’s T(h, a).
@inlinable public func owensT<T: BinaryFloatingPoint>(h: T, a: T) throws -> T {
    let dh = D(h), da = D(a)
    guard dh.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
    guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
    return T(bs_owens_t(dh, da))
}

// MARK: - Exponential integrals and related

/// Generalized exponential integral En(n, x).
@inlinable public func exponentialIntegralEn<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_expint_En(Int32(n), dx))
}

@inlinable public func expm1<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_expm1(dx))
}

@inlinable public func log1p<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // log1p domain: x > -1
    guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
    return T(bs_log1p(dx))
}

@inlinable public func log1pmx<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // log1pmx shares the same domain as log1p
    guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
    return T(bs_log1pmx(dx))
}

@inlinable public func powm1<T: BinaryFloatingPoint>(_ x: T, _ y: T) throws -> T {
    let dx = D(x), dy = D(y)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    // Real power of negative base with non-integer exponent is not real
    let yIsInteger = dy == dy.rounded(.towardZero)
    guard !(dx < 0 && !yIsInteger) else {
        throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
    }
    return T(bs_powm1(dx, dy))
}

@inlinable public func cbrt<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cbrt(dx))
}

// MARK: - Trig helpers

@inlinable public func sinPi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_sin_pi(dx))
}

@inlinable public func cosPi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cos_pi(dx))
}

// MARK: - Airy

@inlinable public func airyAi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_ai(dx))
}

@inlinable public func airyBi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_bi(dx))
}

@inlinable public func airyAiPrime<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_ai_prime(dx))
}

@inlinable public func airyBiPrime<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_bi_prime(dx))
}

// MARK: - Bessel (cylindrical)

/// Cylindrical Bessel J_v(x).
@inlinable public func besselJ<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cyl_bessel_j(dv, dx))
}

/// Cylindrical Neumann Y_v(x).
@inlinable public func besselY<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // Y_v has a singularity at x = 0; principal domain x > 0
    guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return T(bs_cyl_neumann(dv, dx))
}

/// Modified Bessel I_v(x).
@inlinable public func modifiedBesselI<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cyl_bessel_i(dv, dx))
}

/// Modified Bessel K_v(x).
@inlinable public func modifiedBesselK<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // K_v has a singularity at x = 0; principal domain x > 0
    guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return T(bs_cyl_bessel_k(dv, dx))
}

// MARK: - Legendre

/// Legendre polynomial P_n(x).
@inlinable public func legendreP<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_legendre_p(Int32(n), dx))
}

/// Associated Legendre function P_n^m(x) (integer n, m).
@inlinable public func associatedLegendreP<T: BinaryFloatingPoint>(_ n: Int, _ m: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_assoc_legendre_p(Int32(n), Int32(m), dx))
}

// MARK: - Elliptic integrals (Legendre forms)

/// Complete elliptic integral of the first kind K(k).
@inlinable public func completeEllipticIntegralK<T: BinaryFloatingPoint>(_ k: T) throws -> T {
    let dk = D(k)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    // Principal domain for modulus k: |k| <= 1
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_1_complete(dk))
}

/// Incomplete elliptic integral of the first kind F(k, phi).
@inlinable public func incompleteEllipticIntegralF<T: BinaryFloatingPoint>(_ k: T, phi: T) throws -> T {
    let dk = D(k), dphi = D(phi)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_1(dk, dphi))
}

/// Complete elliptic integral of the second kind E(k).
@inlinable public func completeEllipticIntegralE<T: BinaryFloatingPoint>(_ k: T) throws -> T {
    let dk = D(k)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_2_complete(dk))
}

/// Incomplete elliptic integral of the second kind E(k, phi).
@inlinable public func incompleteEllipticIntegralE<T: BinaryFloatingPoint>(_ k: T, phi: T) throws -> T {
    let dk = D(k), dphi = D(phi)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_2(dk, dphi))
}

/// Complete elliptic integral of the third kind Π(n; k).
@inlinable public func completeEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, characteristic n: T) throws -> T {
    let dk = D(k), dn = D(n)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_3_complete(dk, dn))
}

/// Incomplete elliptic integral of the third kind Π(n; k, phi).
@inlinable public func incompleteEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, characteristic n: T, phi: T) throws -> T {
    let dk = D(k), dn = D(n), dphi = D(phi)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_3(dk, dn, dphi))
}

// MARK: - Elliptic integrals (Carlson symmetric forms)

@inlinable public func carlsonRC<T: BinaryFloatingPoint>(_ x: T, _ y: T) throws -> T {
    let dx = D(x), dy = D(y)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    // Principal domain: y > 0 (to avoid singularities)
    guard dy > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y") }
    return T(bs_ellint_rc(dx, dy))
}

@inlinable public func carlsonRF<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    // Typical principal domain: x, y, z >= 0 (non-negative)
    guard dx >= 0, dy >= 0, dz >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rf(dx, dy, dz))
}

@inlinable public func carlsonRD<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dx >= 0, dy >= 0, dz >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rd(dx, dy, dz))
}

@inlinable public func carlsonRJ<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T, _ p: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z), dp = D(p)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dp.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p") }
    guard dx >= 0, dy >= 0, dz >= 0, dp >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rj(dx, dy, dz, dp))
}

@inlinable public func carlsonRG<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dx >= 0, dy >= 0, dz >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rg(dx, dy, dz))
}

// MARK: - Lambert W (real branches)

/// Principal branch W0(x).
@inlinable public func lambertW0<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    // Domain: x >= -1/e
    let minX = -1.0 / boostE
    guard dx >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: Double.infinity) }
    return T(bs_lambert_w0(dx))
}

/// Lower branch W−1(x), defined for x ∈ [−1/e, 0).
@inlinable public func lambertWm1<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 / boostE
    guard dx >= minX && dx < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: 0.0) }
    return T(bs_lambert_wm1(dx))
}
