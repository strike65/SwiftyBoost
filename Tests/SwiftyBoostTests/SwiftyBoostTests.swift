import Testing
@testable import SwiftyBoost
import CBoostBridge

// MARK: - Factorial, Rising Factorial/Pochhammer, Binomial, Double Factorial

@Suite("Factorial, Rising Factorial/Pochhammer, Binomial, Double Factorial")
struct FactorialAndCombinatoricsTests {

    // MARK: - Factorial

    @Test("factorial Double matches backend up to 170")
    func factorialDouble() throws {
        for n in 0...Int(SpecialFunctions.maxFactorialInputDouble) {
            let nn = UInt32(n)
            let expected = bs_factorial(nn)
            let got = try SpecialFunctions.factorial(nn)
            #expect(got == expected, "Mismatch at n=\(n)")
        }
    }

    @Test("factorial Double throws past 170")
    func factorialDoubleOverflow() {
        let n = SpecialFunctions.maxFactorialInputDouble + 1
        #expect(throws: SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(SpecialFunctions.maxFactorialInputDouble))) {
            _ = try SpecialFunctions.factorial(n)
        }
    }

    @Test("factorial Float matches backend up to 34")
    func factorialFloat() throws {
        for n in 0...Int(SpecialFunctions.maxFactorialInputFloat) {
            let nn = UInt32(n)
            let expected = bs_factorial_f(nn)
            let got = try SpecialFunctions.factorial_f(nn)
            #expect(got == expected, "Mismatch at n=\(n)")
        }
    }

    @Test("factorial Float throws past 34")
    func factorialFloatOverflow() {
        let n = SpecialFunctions.maxFactorialInputFloat + 1
        #expect(throws: SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(SpecialFunctions.maxFactorialInputFloat))) {
            _ = try SpecialFunctions.factorial_f(n)
        }
    }

    @Test("factorial generic<T=Float> matches Double-backend converted to Float")
    func factorialGenericFloat() throws {
        typealias T = Float
        for n in 0...Int(SpecialFunctions.maxFactorialInput(for: T.self)) {
            let nn = UInt32(n)
            let expected = T(bs_factorial(nn)) // generic path uses Double backend
            let got: T = try SpecialFunctions.factorial(nn)
            #expect(got == expected, "Mismatch at n=\(n)")
        }
    }

    @Test("factorial generic<T=Double> matches Double backend")
    func factorialGenericDouble() throws {
        typealias T = Double
        for n in 0...Int(SpecialFunctions.maxFactorialInput(for: T.self)) {
            let nn = UInt32(n)
            let expected = T(bs_factorial(nn))
            let got: T = try SpecialFunctions.factorial(nn)
            #expect(got == expected, "Mismatch at n=\(n)")
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("factorial Float80 matches backend up to 1754")
    func factorialFloat80() throws {
        for n in 0...Int(SpecialFunctions.maxFactorialInputFloat80) {
            let nn = UInt32(n)
            let expected = bs_factorial_l(nn)
            let got = try SpecialFunctions.factorial_l(nn)
            #expect(got == expected, "Mismatch at n=\(n)")
        }
    }

    @Test("factorial Float80 throws past 1754")
    func factorialFloat80Overflow() {
        let n = SpecialFunctions.maxFactorialInputFloat80 + 1
        #expect(throws: SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(SpecialFunctions.maxFactorialInputFloat80))) {
            _ = try SpecialFunctions.factorial_l(n)
        }
    }

    @Test("factorial generic<T=Float80> matches Double-backend converted to Float80")
    func factorialGenericFloat80() throws {
        typealias T = Float80
        for n in 0...Int(SpecialFunctions.maxFactorialInput(for: T.self)) {
            let nn = UInt32(n)
            let expected = T(bs_factorial(nn)) // generic path uses Double backend
            let got: T = try SpecialFunctions.factorial(nn)
            #expect(got == expected, "Mismatch at n=\(n)")
        }
    }
    #endif

    // MARK: - Rising factorial / Pochhammer

    @Test("rising_factorial Double matches backend, includes edge cases")
    func risingFactorialDouble() throws {
        let cases: [(Double, UInt32)] = [
            (0.0, 0), (1.0, 0), (-3.0, 0),
            (2.5, 1), (2.5, 5), (-2.5, 5),
            (-3.0, 3), // not exact zero
            (-3.0, 4), // exact zero (k in [0,3], −3 + 3 = 0)
        ]
        for (x, n) in cases {
            let isZero = {
                if n == 0 { return false }
                if x <= 0 {
                    let k = -x
                    let ktrunc = k.rounded(.towardZero)
                    if k == ktrunc {
                        return k >= 0 && k < Double(n)
                    }
                }
                return false
            }()
            if isZero {
                let got = try SpecialFunctions.rising_factorial(x, n)
                #expect(got == 0)
            } else {
                let got = try SpecialFunctions.rising_factorial(x, n)
                let expected = bs_rising_factorial(x, n)
                #expect(got == expected, "Mismatch at x=\(x), n=\(n)")
            }
        }
    }

    @Test("pochhammer Double equals rising_factorial")
    func pochhammerDoubleAlias() throws {
        let xs: [Double] = [-2.5, -1.0, 0, 0.75, 3.0]
        let ns: [UInt32] = [0, 1, 2, 5, 10]
        for x in xs {
            for n in ns {
                let a = try SpecialFunctions.pochhammer(x, n)
                let b = try SpecialFunctions.rising_factorial(x, n)
                #expect(a == b, "Alias mismatch at x=\(x), n=\(n)")
            }
        }
    }

    @Test("rising_factorial Float matches backend for representative cases")
    func risingFactorialFloat() throws {
        let cases: [(Float, UInt32)] = [(0,0), (1,0), (2.5,3), (-2.5,3), (-3,4)]
        for (x, n) in cases {
            let isZero: Bool = {
                if n == 0 { return false }
                if x <= 0 {
                    let k = -Double(x)
                    let ktrunc = k.rounded(.towardZero)
                    if k == ktrunc {
                        return k >= 0 && k < Double(n)
                    }
                }
                return false
            }()
            if isZero {
                let got = try SpecialFunctions.rising_factorial_f(x, n)
                #expect(got == 0)
            } else {
                let got = try SpecialFunctions.rising_factorial_f(x, n)
                let expected = bs_rising_factorial_f(x, n)
                #expect(got == expected, "Mismatch at x=\(x), n=\(n)")
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("rising_factorial Float80 matches backend for representative cases")
    func risingFactorialFloat80() throws {
        let cases: [(Float80, UInt32)] = [(0,0), (1,0), (2.5,3), (-2.5,3), (-3,4)]
        for (x, n) in cases {
            let isZero: Bool = {
                if n == 0 { return false }
                if x <= 0 {
                    let k = -Double(x)
                    let ktrunc = k.rounded(.towardZero)
                    if k == ktrunc {
                        return k >= 0 && k < Double(n)
                    }
                }
                return false
            }()
            if isZero {
                let got = try SpecialFunctions.rising_factorial_l(x, n)
                #expect(got == 0)
            } else {
                let got = try SpecialFunctions.rising_factorial_l(x, n)
                let expected = bs_rising_factorial_l(x, n)
                #expect(got == expected, "Mismatch at x=\(x), n=\(n)")
            }
        }
    }
    #endif

    // MARK: - Binomial coefficient

    @Test("binomial_coeff Double matches backend and special cases")
    func binomialDouble() throws {
        #expect(try SpecialFunctions.binomial_coeff(5, 0) == 1)
        #expect(try SpecialFunctions.binomial_coeff(5, 5) == 1)
        #expect(try SpecialFunctions.binomial_coeff(5, 6) == 0)

        let pairs: [(UInt32, UInt32)] = [(5,1), (5,2), (10,3), (20,10)]
        for (n, k) in pairs {
            let got = try SpecialFunctions.binomial_coeff(n, k)
            let expected = bs_binomial_coefficient(n, k)
            #expect(got == expected, "Mismatch at C(\(n), \(k))")
        }
    }

    @Test("binomial_coeff Float matches backend and special cases")
    func binomialFloat() throws {
        #expect(try SpecialFunctions.binomial_coeff_f(5, 0) == 1)
        #expect(try SpecialFunctions.binomial_coeff_f(5, 5) == 1)
        #expect(try SpecialFunctions.binomial_coeff_f(5, 6) == 0)

        let pairs: [(UInt32, UInt32)] = [(5,1), (5,2), (10,3), (20,10)]
        for (n, k) in pairs {
            let got = try SpecialFunctions.binomial_coeff_f(n, k)
            let expected = bs_binomial_coefficient_f(n, k)
            #expect(got == expected, "Mismatch at C(\(n), \(k))")
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("binomial_coeff Float80 matches backend and special cases")
    func binomialFloat80() throws {
        #expect(try SpecialFunctions.binomial_coeff_l(5, 0) == 1)
        #expect(try SpecialFunctions.binomial_coeff_l(5, 5) == 1)
        #expect(try SpecialFunctions.binomial_coeff_l(5, 6) == 0)

        let pairs: [(UInt32, UInt32)] = [(5,1), (5,2), (10,3), (20,10)]
        for (n, k) in pairs {
            let got = try SpecialFunctions.binomial_coeff_l(n, k)
            let expected = bs_binomial_coefficient_l(n, k)
            #expect(got == expected, "Mismatch at C(\(n), \(k))")
        }
    }
    #endif

    // MARK: - Double factorial

    @Test("double_factorial Double matches backend for representative n")
    func doubleFactorialDouble() throws {
        let ns: [UInt32] = Array(0...25)
        for n in ns {
            let got = try SpecialFunctions.double_factorial(n)
            let expected = bs_double_factorial(n)
            #expect(got == expected, "Mismatch at n!! with n=\(n)")
        }
    }

    @Test("double_factorial Float matches backend for representative n")
    func doubleFactorialFloat() throws {
        let ns: [UInt32] = Array(0...25)
        for n in ns {
            let got = try SpecialFunctions.double_factorial_f(n)
            let expected = bs_double_factorial_f(n)
            #expect(got == expected, "Mismatch at n!! with n=\(n)")
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("double_factorial Float80 matches backend for representative n")
    func doubleFactorialFloat80() throws {
        let ns: [UInt32] = Array(0...25)
        for n in ns {
            let got = try SpecialFunctions.double_factorial_l(n)
            let expected = bs_double_factorial_l(n)
            #expect(got == expected, "Mismatch at n!! with n=\(n)")
        }
    }
    #endif
}

// MARK: - Gamma, LogGamma, Ratios, Incomplete/Regularized, Inverses, Derivative

@Suite("Gamma, LogGamma, Ratios, Incomplete/Regularized, Inverses, Derivative")
struct GammaAndRelatedTests {

    @Test
    func gammaDouble() throws {
        let xs: [Double] = [0.5, 1, 2.5, 5, 10]
        for x in xs {
            let got = try SpecialFunctions.gamma(x)
            let expected = bs_tgamma(x)
            #expect(got == expected, "Γ mismatch at x=\(x)")
        }
        for x in [0.0, -1.0, -2.0] {
            #expect(throws: SpecialFunctionError.poleAtNonPositiveInteger(name: "x")) {
                _ = try SpecialFunctions.gamma(x)
            }
        }
    }

    @Test
    func logGammaDouble() throws {
        let xs: [Double] = [0.5, 1, 2.5, 5, 10]
        for x in xs {
            let got = try SpecialFunctions.logGamma(x)
            let expected = bs_lgamma(x)
            #expect(got == expected, "lnΓ mismatch at x=\(x)")
        }
    }

    @Test
    func gammaRatios() throws {
        let a = 5.5, b = 2.25
        let r = try SpecialFunctions.gammaRatio(a, b)
        #expect(r == bs_tgamma_ratio(a, b))

        let delta = -1.75
        let rd = try SpecialFunctions.gammaDeltaRatio(a, delta: delta)
        #expect(rd == bs_tgamma_delta_ratio(a, delta))
    }

    @Test
    func incompleteAndRegularized() throws {
        let a = 2.5, x = 1.25
        #expect(try SpecialFunctions.incompleteGammaLower(a, x: x) == bs_tgamma_lower(a, x))
        #expect(try SpecialFunctions.incompleteGammaUpper(a, x: x) == bs_tgamma_upper(a, x))
        #expect(try SpecialFunctions.regularizedGammaP(a, x: x) == bs_gamma_p(a, x))
        #expect(try SpecialFunctions.regularizedGammaQ(a, x: x) == bs_gamma_q(a, x))
        #expect(try SpecialFunctions.regularizedGammaPDerivative(a, x: x) == bs_gamma_p_derivative(a, x))
    }

    @Test
    func regularizedInverses() throws {
        let a = 3.5
        for p in [0.1, 0.5, 0.9] {
            #expect(try SpecialFunctions.regularizedGammaPInv(a, p: p) == bs_gamma_p_inv(a, p))
        }
        for q in [0.1, 0.5, 0.9] {
            #expect(try SpecialFunctions.regularizedGammaQInv(a, q: q) == bs_gamma_q_inv(a, q))
        }
    }
}

// MARK: - Beta, Incomplete/Regularized, Inverses, Solvers, Derivative

@Suite("Beta, Incomplete/Regularized, Inverses, Solvers, Derivative")
struct BetaFamilyTests {

    @Test
    func betaDouble() throws {
        let a = 2.5, b = 3.75
        #expect(try SpecialFunctions.beta(a, b) == bs_beta(a, b))
    }

    @Test
    func incompleteAndRegularized() throws {
        let a = 2.5, b = 3.75
        for x in [0.0, 0.25, 0.5, 1.0] {
            #expect(try SpecialFunctions.incompleteBetaUnnormalized(a, b, x: x) == bs_fullBeta(a, b, x))
            #expect(try SpecialFunctions.regularizedIncompleteBeta(a, b, x: x) == bs_ibeta(a, b, x))
            #expect(try SpecialFunctions.complementaryRegularizedIncompleteBeta(a, b, x: x) == bs_ibetac(a, b, x))
            #expect(try SpecialFunctions.regularizedIncompleteBetaDerivative(a, b, x: x) == bs_ibeta_derivative(a, b, x))
        }
    }

    @Test
    func inversesAndSolvers() throws {
        let a = 2.5, b = 3.75
        for p in [0.1, 0.5, 0.9] {
            #expect(try SpecialFunctions.inverseRegularizedIncompleteBeta(a, b, p: p) == bs_ibeta_inv(a, b, p))
            #expect(try SpecialFunctions.inverseComplementaryRegularizedIncompleteBeta(a, b, p: p) == bs_ibetac_inv(a, b, p))
        }
        let x = 0.4, p = 0.6
        #expect(SpecialFunctions.solveAForRegularizedIncompleteBeta(b: b, x: x, p: p) == bs_ibeta_inva(b, x, p))
        #expect(SpecialFunctions.solveBForRegularizedIncompleteBeta(a: a, x: x, p: p) == bs_ibeta_invb(a, x, p))
    }
}

// MARK: - Bessel (J, Y, I, K) and Legendre (P, P^m)

@Suite("Bessel (J, Y, I, K) and Legendre (P, P^m)")
struct BesselAndLegendreTests {

    @Test
    func besselJYIK_Double() throws {
        let v = 2.5
        let xs: [Double] = [0.1, 1.0, 5.0]
        for x in xs {
            #expect(try SpecialFunctions.besselJ(v: v, x: x) == bs_cyl_bessel_j(v, x))
            #expect(try SpecialFunctions.modifiedBesselI(v: v, x: x) == bs_cyl_bessel_i(v, x))
        }
        for x in xs {
            #expect(throws: SpecialFunctionError.parameterNotPositive(name: "x")) {
                _ = try SpecialFunctions.besselY(v: v, x: -abs(x))
            }
            #expect(throws: SpecialFunctionError.parameterNotPositive(name: "x")) {
                _ = try SpecialFunctions.modifiedBesselK(v: v, x: 0.0)
            }
            #expect(try SpecialFunctions.besselY(v: v, x: x) == bs_cyl_neumann(v, x))
            #expect(try SpecialFunctions.modifiedBesselK(v: v, x: x) == bs_cyl_bessel_k(v, x))
        }
    }

    @Test
    func legendreP_and_Associated_Double() throws {
        let xs: [Double] = [-0.75, -0.1, 0.0, 0.25, 0.9]
        for n in 0...6 {
            for x in xs {
                #expect(try SpecialFunctions.legendreP(n, x) == bs_legendre_p(Int32(n), x))
            }
        }
        for n in 0...6 {
            for m in -n...n {
                for x in xs {
                    #expect(try SpecialFunctions.associatedLegendreP(n, m, x) == bs_assoc_legendre_p(Int32(n), Int32(m), x))
                }
            }
        }
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.legendreP(-1, 0.1 as Double)
        }
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-3), max: Double(3))) {
            _ = try SpecialFunctions.associatedLegendreP(3, 4, 0.1 as Double)
        }
    }
}

// MARK: - Legendre Elliptic Integrals and Lambert W

@Suite("Legendre Elliptic Integrals and Lambert W")
struct EllipticAndLambertWTests {

    @Test
    func elliptic_Double() throws {
        let ks: [Double] = [-0.9, 0.0, 0.5, 0.9, 1.0]
        let phis: [Double] = [-1.2, -0.3, 0.0, 0.7, 1.3]
        for k in ks {
            #expect(try SpecialFunctions.completeEllipticIntegralK(k) == bs_ellint_1_complete(k))
            #expect(try SpecialFunctions.completeEllipticIntegralE(k) == bs_ellint_2_complete(k))
            for phi in phis {
                #expect(try SpecialFunctions.incompleteEllipticIntegralF(k, phi: phi) == bs_ellint_1(k, phi))
                #expect(try SpecialFunctions.incompleteEllipticIntegralE(k, phi: phi) == bs_ellint_2(k, phi))
            }
            #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)) {
                _ = try SpecialFunctions.completeEllipticIntegralK(1.1)
            }
        }
    }

    @Test
    func ellipticPi_Double() throws {
        let ks: [Double] = [-0.9, 0.0, 0.5, 0.9]
        let nus: [Double] = [-0.5, 0.0, 0.2, 0.8]
        let phis: [Double] = [-1.2, -0.3, 0.0, 0.7, 1.3]
        for k in ks {
            for nu in nus {
                for phi in phis {
                    let got = SpecialFunctions.incompleteEllipticIntegralPi(k, nu, phi)
                    let expected = bs_ellint_3(k, nu, phi)
                    #expect(got == expected, "Π mismatch at k=\(k), nu=\(nu), phi=\(phi)")
                }
            }
            // Identity check: Π(0; φ | k) = F(φ | k)
            for phi in phis {
                let pi0 = SpecialFunctions.incompleteEllipticIntegralPi(k, 0.0, phi)
                let f = try SpecialFunctions.incompleteEllipticIntegralF(k, phi: phi)
                let diff = abs(pi0 - f)
                #expect(diff <= 1e-15, "Identity Π(0; φ | k) = F(φ|k) failed at k=\(k), phi=\(phi); diff=\(diff)")
            }
        }
    }

    @Test
    func ellipticPiComplete_Double() throws {
        let ks: [Double] = [-0.9, 0.0, 0.5, 0.9]
        let nus: [Double] = [-0.5, 0.0, 0.2, 0.8]
        let halfPi = Double.pi / 2
        for k in ks {
            for nu in nus {
                // Backend comparison
                let got = SpecialFunctions.completeEllipticIntegralPi(k, nu)
                let expected = bs_ellint_3_complete(k, nu)
                #expect(got == expected, "Π_complete mismatch at k=\(k), nu=\(nu)")

                // Identity: Π(n | k) = Π(n; π/2 | k)
                let incAtHalfPi = SpecialFunctions.incompleteEllipticIntegralPi(k, nu, halfPi)
                let diff = abs(got - incAtHalfPi)
                #expect(diff <= 1e-15, "Identity Π(n|k) = Π(n; π/2 | k) failed at k=\(k), nu=\(nu); diff=\(diff)")
            }
        }
    }

    @Test
    func lambertW_Double() throws {
        let e = bs_const_e()
        let minX = -1.0 / e
        let xs0: [Double] = [minX, -0.2, 0.0, 0.5, 2.0]
        for x in xs0 {
            if x >= minX {
                #expect(try SpecialFunctions.lambertW0(x) == bs_lambert_w0(x))
            }
        }
        let xsm1: [Double] = [minX, -0.9, -0.2]
        for x in xsm1 {
            if x >= minX && x < 0 {
                #expect(try SpecialFunctions.lambertWm1(x) == bs_lambert_wm1(x))
            } else {
                #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: 0.0)) {
                    _ = try SpecialFunctions.lambertWm1(x)
                }
            }
        }

        // Do/catch assertions to avoid "Caught error" issue lines:

        // W0 should throw for x < -1/e
        do {
            _ = try SpecialFunctions.lambertW0(minX - 1e-6)
            #expect(Bool(false), "lambertW0 should throw for x < -1/e")
        } catch let e as SpecialFunctionError {
            switch e {
            case .parameterOutOfRange(let name, let min, let max):
                #expect(name == "x")
                #expect(abs(min - minX) <= 1e-15)
                #expect(max == Double.infinity)
            default:
                #expect(Bool(false), "Unexpected error: \(e)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }

        // W−1 should throw for x >= 0
        do {
            _ = try SpecialFunctions.lambertWm1(0.1)
            #expect(Bool(false), "lambertWm1 should throw for x >= 0")
        } catch let e as SpecialFunctionError {
            switch e {
            case .parameterOutOfRange(let name, let min, let max):
                #expect(name == "x")
                #expect(abs(min - minX) <= 1e-15)
                #expect(max == 0.0)
            default:
                #expect(Bool(false), "Unexpected error: \(e)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }
}

// MARK: - Exponential Integrals Ei/En and erf/erfc

@Suite("Exponential Integrals Ei/En and erf/erfc")
struct ExponentialIntegralsAndErrorFunctionTests {

    @Test
    func exponentialIntegralEi_Double() throws {
        let xs: [Double] = [-5, -1, -0.1, 0.1, 1, 5]
        for x in xs {
            #expect(try SpecialFunctions.exponentialIntegralEi(x) == bs_expint_Ei(x))
        }
    }

    @Test
    func exponentialIntegralEn_Double() throws {
        let ns = [0, 1, 2, 5]
        let xs: [Double] = [0.0, 0.1, 1.0, 5.0]
        for n in ns {
            for x in xs {
                #expect(try SpecialFunctions.exponentialIntegralEn(n, x) == bs_expint_En(Int32(n), x))
            }
        }
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.exponentialIntegralEn(-1, 1.0 as Double)
        }
    }

    @Test
    func errorFunctions_Double() throws {
        let xs: [Double] = [-3, -1, -0.1, 0, 0.1, 1, 3]
        for x in xs {
            #expect(try SpecialFunctions.errorFunction(x) == bs_erf(x))
            #expect(try SpecialFunctions.complementaryErrorFunction(x) == bs_erfc(x))
        }
    }
}

// MARK: - Digamma, Trigamma, Polygamma, Riemann Zeta

@Suite("Digamma, Trigamma, Polygamma, Riemann Zeta")
struct DigammaPolygammaZetaTests {

    @Test
    func digamma_trigamma_Double() throws {
        let xs: [Double] = [0.5, 0.75, 1.0, 2.5, 5.0]
        for x in xs {
            #expect(try SpecialFunctions.digamma(x) == bs_digamma(x))
            #expect(try SpecialFunctions.trigamma(x) == bs_trigamma(x))
        }
        for x in [0.0, -1.0, -2.0] {
            #expect(throws: SpecialFunctionError.poleAtNonPositiveInteger(name: "x")) {
                _ = try SpecialFunctions.digamma(x)
            }
            #expect(throws: SpecialFunctionError.poleAtNonPositiveInteger(name: "x")) {
                _ = try SpecialFunctions.trigamma(x)
            }
        }
    }

    @Test
    func polygamma_Double() throws {
        let xs: [Double] = [0.5, 1.0, 2.5]
        for n in 0...5 {
            for x in xs {
                #expect(try SpecialFunctions.polygamma(order: n, x) == bs_polygamma(Int32(n), x))
            }
        }
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "order")) {
            _ = try SpecialFunctions.polygamma(order: -1, 1.0 as Double)
        }
    }

    @Test
    func zeta_Double() throws {
        let xs: [Double] = [-3.0, -2.0, -0.5, 0.0, 2.0, 3.5]
        for x in xs where x != 1.0 {
            #expect(try SpecialFunctions.riemannZeta(x) == bs_riemann_zeta(x))
        }
        #expect(throws: SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1")) {
            _ = try SpecialFunctions.riemannZeta(1.0 as Double)
        }
    }
}

// MARK: - Owen's T and Hypergeometric functions

@Suite("Owen's T and Hypergeometric functions")
struct OwensTAndHypergeometricTests {

    @Test
    func owensT_Double() throws {
        let hs: [Double] = [-2.0, -0.5, 0.0, 0.5, 2.0]
        let as_: [Double] = [-2.0, -0.5, 0.0, 0.5, 2.0]
        for h in hs {
            for a in as_ {
                #expect(try SpecialFunctions.owensT(h: h, a: a) == bs_owens_t(h, a))
            }
        }
    }

    @Test
    func hypergeometric_Scalars_Double() throws {
        let a = 1.25
        for z in [-0.5, -0.1, 0.0, 0.25, 0.75] as [Double] {
            #expect(try SpecialFunctions.hypergeometric1F0(a: a, z: z) == bs_hypergeometric_1F0(a, z))
        }

        let b = 2.75
        for z in [-2.0, -0.5, 0.0, 0.5, 2.0] as [Double] {
            #expect(try SpecialFunctions.hypergeometric0F1(b: b, z: z) == bs_hypergeometric_0F1(b, z))
        }

        let a2 = -1.5, b2 = 0.75
        for z in [-0.5, -0.1, 0.0, 0.1, 0.5] as [Double] {
            #expect(try SpecialFunctions.hypergeometric2F0(a: a2, b: b2, z: z) == bs_hypergeometric_2F0(a2, b2, z))
        }

        let a1 = 1.5, b1 = 2.25
        for z in [-1.0, -0.25, 0.0, 0.5, 2.0] as [Double] {
            #expect(try SpecialFunctions.hypergeometric1F1(a: a1, b: b1, z: z) == bs_hypergeometric_1F1(a1, b1, z))
        }
    }

    @Test
    func hypergeometric_pFq_Double() {
        let a: [Double] = [1.0, 2.0]
        let b: [Double] = [3.0, 4.0]
        let z: Double = 0.25
        let got = SpecialFunctions.hypergeometricPFQ(a: a, b: b, z: z)
        let expected = a.withUnsafeBufferPointer { ap in
            b.withUnsafeBufferPointer { bp in
                bs_hypergeometric_pFq(ap.baseAddress, ap.count, bp.baseAddress, bp.count, z)
            }
        }
        #expect(got == expected)
    }
}

// MARK: - Common helpers: expm1, log1p, log1pmx, powm1, cbrt, sinPi, cosPi

@Suite("Common helpers: expm1, log1p, log1pmx, powm1, cbrt, sinPi, cosPi")
struct CommonHelpersTests {

    @Test
    func expm1_log1p_log1pmx_cbrt_Double() throws {
        let xs: [Double] = [-0.9, -0.1, 0.0, 0.1, 1.0]
        for x in xs {
            if x > -1 {
                #expect(try SpecialFunctions.expm1(x) == bs_expm1(x))
                #expect(try SpecialFunctions.log1p(x) == bs_log1p(x))
                #expect(try SpecialFunctions.log1pmx(x) == bs_log1pmx(x))
            } else {
                #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity)) {
                    _ = try SpecialFunctions.log1p(x)
                }
                #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity)) {
                    _ = try SpecialFunctions.log1pmx(x)
                }
            }
            #expect(try SpecialFunctions.cbrt(x) == bs_cbrt(x))
        }
    }

    @Test
    func powm1_Double() throws {
        let cases: [(Double, Double)] = [(1.000001, 3.0), (2.0, -1.5), (0.5, 2.5), (-2.0, -2.0)]
        for (x, y) in cases {
            #expect(try SpecialFunctions.powm1(x, y) == bs_powm1(x, y))
        }
        #expect(throws: SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")) {
            _ = try SpecialFunctions.powm1(-2.0 as Double, 0.5 as Double)
        }
    }

    @Test
    func sinPi_cosPi_Double() throws {
        let xs: [Double] = [-2.5, -1.0, -0.25, 0.0, 0.25, 1.0, 2.5]
        for x in xs {
            #expect(try SpecialFunctions.sinPi(x) == bs_sin_pi(x))
            #expect(try SpecialFunctions.cosPi(x) == bs_cos_pi(x))
        }
    }
}
