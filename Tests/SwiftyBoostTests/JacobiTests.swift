import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Jacobi Polynomials")
struct JacobiTests {

    @Test("Double evaluations and derivatives match bridge")
    func doubleMatchesBridge() throws {
        let xs: [Double] = [-0.95, -0.2, 0.0, 0.4, 0.9]
        let params: [(Double, Double)] = [
            (0.0, 0.0),
            (0.25, -0.3),
            (1.5, 0.5)
        ]
        for n in 0...8 {
            let nu = UInt32(n)
            for (alpha, beta) in params {
                for x in xs {
                    let expected = bs_jacobi_d(nu, alpha, beta, x)
                    let got = try SpecialFunctions.jacobi(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(got == expected)

                    let expPrime = bs_jacobi_prime_d(nu, alpha, beta, x)
                    let gotPrime = try SpecialFunctions.jacobiPrime(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(gotPrime == expPrime)

                    let expDoublePrime = bs_jacobi_double_prime_d(nu, alpha, beta, x)
                    let gotDoublePrime = try SpecialFunctions.jacobiDoublePrime(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(gotDoublePrime == expDoublePrime)

                    for k in 0...min(3, n) {
                        let expDeriv = bs_jacobi_derivative_d(nu, alpha, beta, x, UInt32(k))
                        let gotDeriv = try SpecialFunctions.jacobiDerivative(n: n, alpha: alpha, beta: beta, x: x, k: k)
                        #expect(gotDeriv == expDeriv)
                        if k == 1 {
                            #expect(abs(gotDeriv - gotPrime) <= 1e-12)
                        } else if k == 2 {
                            #expect(abs(gotDeriv - gotDoublePrime) <= 1e-12)
                        }
                    }
                }
            }
        }
    }

    @Test("Float overloads mirror bridge")
    func floatMatchesBridge() throws {
        let xs: [Float] = [-0.8, -0.3, 0.0, 0.45, 0.85]
        let params: [(Float, Float)] = [
            (0.0, 0.0),
            (0.5, -0.25)
        ]
        for n in 0...6 {
            let nu = UInt32(n)
            for (alpha, beta) in params {
                for x in xs {
                    let expected = bs_jacobi_f(nu, alpha, beta, x)
                    let got = try SpecialFunctions.jacobi(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(got == expected)

                    let expPrime = bs_jacobi_prime_f(nu, alpha, beta, x)
                    let gotPrime = try SpecialFunctions.jacobiPrime(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(gotPrime == expPrime)

                    let expDoublePrime = bs_jacobi_double_prime_f(nu, alpha, beta, x)
                    let gotDoublePrime = try SpecialFunctions.jacobiDoublePrime(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(gotDoublePrime == expDoublePrime)

                    let expDeriv = bs_jacobi_derivative_f(nu, alpha, beta, x, 2)
                    let gotDeriv = try SpecialFunctions.jacobiDerivative(n: n, alpha: alpha, beta: beta, x: x, k: 2)
                    #expect(gotDeriv == expDeriv)
                }
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 overloads mirror bridge")
    func float80MatchesBridge() throws {
        let xs: [Float80] = [-0.75, -0.25, 0.0, 0.5, 0.92]
        let params: [(Float80, Float80)] = [
            (0.0, 0.0),
            (0.75, -0.4)
        ]
        for n in 0...5 {
            let nu = UInt32(n)
            for (alpha, beta) in params {
                for x in xs {
                    let expected = bs_jacobi_l(nu, alpha, beta, x)
                    let got = try SpecialFunctions.jacobi(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(got == expected)

                    let expPrime = bs_jacobi_prime_l(nu, alpha, beta, x)
                    let gotPrime = try SpecialFunctions.jacobiPrime(n: n, alpha: alpha, beta: beta, x: x)
                    #expect(gotPrime == expPrime)

                    let expDerivative = bs_jacobi_derivative_l(nu, alpha, beta, x, 2)
                    let gotDerivative = try SpecialFunctions.jacobiDerivative(n: n, alpha: alpha, beta: beta, x: x, k: 2)
                    #expect(gotDerivative == expDerivative)
                }
            }
        }
    }
    #endif

    @Test("Rejects invalid parameters")
    func rejectsInvalidInputs() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "n", value: -1)) {
            _ = try SpecialFunctions.jacobi(n: -1, alpha: 0.0, beta: 0.0, x: 0.2)
        }
        do {
            _ = try SpecialFunctions.jacobi(n: 2, alpha: Double.nan, beta: 0.1, x: 0.2)
            #expect(Bool(false), "Expected parameterNotFinite for alpha")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actual, value: value) = error {
                #expect(actual == "alpha")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "k", value: -1)) {
            _ = try SpecialFunctions.jacobiDerivative(n: 3, alpha: 0.1, beta: 0.2, x: 0.3, k: -1)
        }
        do {
            _ = try SpecialFunctions.jacobiDerivative(n: 3, alpha: 0.1, beta: Double.infinity, x: 0.3, k: 1)
            #expect(Bool(false), "Expected parameterNotFinite for beta")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actual, value: value) = error {
                #expect(actual == "beta")
                #expect(value.isInfinite)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }
}
