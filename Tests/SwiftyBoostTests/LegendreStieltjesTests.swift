import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Legendre–Stieltjes Polynomials")
struct LegendreStieltjesTests {

    @Test("Double evaluations match bridge")
    func evaluationMatchesBackendDouble() throws {
        let xs: [Double] = [-0.9, -0.25, 0.0, 0.4, 0.95]
        for m in 1...6 {
            let mu = UInt32(m)
            for x in xs {
                let expected = bs_legendre_stieltjes_d(mu, x)
                let got = try SpecialFunctions.legendreStieltjes(m, x)
                #expect(got == expected)

                let expectedPrime = bs_legendre_stieltjes_prime_d(mu, x)
                let gotPrime = try SpecialFunctions.legendreStieltjesPrime(m, x)
                #expect(gotPrime == expectedPrime)
            }
            let expectedNorm = bs_legendre_stieltjes_norm_sq_d(mu)
            let gotNorm = try SpecialFunctions.legendreStieltjesNormSquared(m) as Double
            #expect(gotNorm == expectedNorm)
        }
    }

    @Test("Zeros for Double match bridge and respect parity")
    func zerosMatchBackendDouble() throws {
        for m in 1...7 {
            let mu = UInt32(m)
            let expectedCount = Int(bs_legendre_stieltjes_zeros_d(mu, nil, 0))
            var expected = Array<Double>(repeating: .zero, count: expectedCount)
            expected.withUnsafeMutableBufferPointer { buf in
                _ = bs_legendre_stieltjes_zeros_d(mu, buf.baseAddress, expectedCount)
            }
            let got: [Double] = try SpecialFunctions.legendreStieltjesZeros(order: m)
            #expect(got.count == expected.count)
            for (lhs, rhs) in zip(got, expected) {
                #expect(lhs == rhs)
            }
            if m % 2 == 1 {
                #expect(got.contains(0.0))
            } else {
                #expect(!got.contains(0.0))
            }
        }
    }

    @Test("Float overloads hit native bridge")
    func floatMatchesBackend() throws {
        let xs: [Float] = [-0.7, -0.2, 0.0, 0.45, 0.9]
        for m in 1...5 {
            let mu = UInt32(m)
            for x in xs {
                let expected = bs_legendre_stieltjes_f(mu, x)
                let got = try SpecialFunctions.legendreStieltjes(m, x)
                #expect(got == expected)
                let expectedPrime = bs_legendre_stieltjes_prime_f(mu, x)
                let gotPrime = try SpecialFunctions.legendreStieltjesPrime(m, x)
                #expect(gotPrime == expectedPrime)
            }
            let expectedNorm = bs_legendre_stieltjes_norm_sq_f(mu)
            let gotNorm = try SpecialFunctions.legendreStieltjesNormSquared(m) as Float
            #expect(gotNorm == expectedNorm)

            let expectedCount = Int(bs_legendre_stieltjes_zeros_f(mu, nil, 0))
            var expected = Array<Float>(repeating: .zero, count: expectedCount)
            expected.withUnsafeMutableBufferPointer { buf in
                _ = bs_legendre_stieltjes_zeros_f(mu, buf.baseAddress, expectedCount)
            }
            let got: [Float] = try SpecialFunctions.legendreStieltjesZeros(order: m)
            #expect(got == expected)
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 overloads mirror bridge")
    func float80MatchesBackend() throws {
        let xs: [Float80] = [-0.8, -0.3, 0.0, 0.35, 0.88]
        for m in 1...4 {
            let mu = UInt32(m)
            for x in xs {
                let expected = bs_legendre_stieltjes_l(mu, x)
                let got = try SpecialFunctions.legendreStieltjes(m, x)
                #expect(got == expected)
                let expectedPrime = bs_legendre_stieltjes_prime_l(mu, x)
                let gotPrime = try SpecialFunctions.legendreStieltjesPrime(m, x)
                #expect(gotPrime == expectedPrime)
            }
            let expectedNorm = bs_legendre_stieltjes_norm_sq_l(mu)
            let gotNorm = try SpecialFunctions.legendreStieltjesNormSquared(m) as Float80
            #expect(gotNorm == expectedNorm)

            let expectedCount = Int(bs_legendre_stieltjes_zeros_l(mu, nil, 0))
            var expected = Array<Float80>(repeating: .zero, count: expectedCount)
            expected.withUnsafeMutableBufferPointer { buf in
                _ = bs_legendre_stieltjes_zeros_l(mu, buf.baseAddress, expectedCount)
            }
            let got: [Float80] = try SpecialFunctions.legendreStieltjesZeros(order: m)
            #expect(got == expected)
        }
    }
    #endif

    @Test("Rejects invalid parameters")
    func rejectsInvalidInputs() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "m", value: -1)) {
            _ = try SpecialFunctions.legendreStieltjes(-1, 0.2 as Double)
        }
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "m", value: 0)) {
            _ = try SpecialFunctions.legendreStieltjes(0, 0.1 as Double)
        }
        do {
            _ = try SpecialFunctions.legendreStieltjes(2, Double.nan)
            #expect(Bool(false), "Expected parameterNotFinite for x=nan")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actualName, value: value) = error {
                #expect(actualName == "x")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        do {
            _ = try SpecialFunctions.legendreStieltjes(2, Double.infinity)
            #expect(Bool(false), "Expected parameterNotFinite for x=+∞")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actualName, value: value) = error {
                #expect(actualName == "x")
                #expect(value.isInfinite)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "m", value: 0)) {
            _ = try SpecialFunctions.legendreStieltjesZeros(order: 0) as [Double]
        }
        #expect(throws: SpecialFunctionError<Float>.parameterNotPositive(name: "m", value: 0)) {
            _ = try SpecialFunctions.legendreStieltjesZeros(order: 0) as [Float]
        }
    }
}
