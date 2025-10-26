import Testing
@testable import SwiftyBoost
import Foundation

@Suite("KL Divergence – Dynamic distributions")
struct DynamicKLDivergenceTests {
    @Test("Exponential vs Exponential matches analytic form")
    func dynamicExponential() throws {
        let lambdaP: Double = 1.75
        let lambdaQ: Double = 0.65
        let p = try Distribution.Dynamic<Double>(
            distributionName: "exponential",
            parameters: ["lambda": lambdaP]
        )
        let q = try Distribution.Dynamic<Double>(
            distributionName: "exponential",
            parameters: ["lambda": lambdaQ]
        )
        let expected = Foundation.log(lambdaP / lambdaQ) + (lambdaQ / lambdaP) - 1
        let maybeKL = try p.klDivergence(relativeTo: q, options: .automatic())
        let kl = try #require(maybeKL)
        #expect(abs(kl - expected) <= 5e-10)
    }

    @Test("Bernoulli vs Bernoulli matches analytic form")
    func dynamicBernoulli() throws {
        let p1: Double = 0.27
        let p2: Double = 0.63
        let p = try Distribution.Dynamic<Double>(
            distributionName: "bernoulli",
            parameters: ["p": p1]
        )
        let q = try Distribution.Dynamic<Double>(
            distributionName: "bernoulli",
            parameters: ["p": p2]
        )
        let expected =
            p1 * Foundation.log(p1 / p2)
            + (1 - p1) * Foundation.log((1 - p1) / (1 - p2))
        let maybeKL = try p.klDivergence(relativeTo: q, options: .automatic())
        let kl = try #require(maybeKL)
        #expect(abs(kl - expected) <= 5e-12)
    }
}

@Suite("KL Divergence – Typed wrappers")
struct TypedKLDivergenceTests {
    @Test("Exponential wrapper delegates to dynamic backend")
    func typedExponential() throws {
        let lhs = try Distribution.Exponential<Double>(lambda: 1.2)
        let rhs = try Distribution.Exponential<Double>(lambda: 0.9)
        let expected = Foundation.log(1.2 / 0.9) + (0.9 / 1.2) - 1
        let maybeKL = try lhs.klDivergence(relativeTo: rhs)
        let kl = try #require(maybeKL)
        #expect(abs(kl - expected) <= 1e-9)
    }

    @Test("Bernoulli wrapper delegates to dynamic backend")
    func typedBernoulli() throws {
        let p1: Double = 0.33
        let p2: Double = 0.77
        let lhs = try Distribution.Bernoulli<Double>(p: p1)
        let rhs = try Distribution.Bernoulli<Double>(p: p2)
        let expected =
            p1 * Foundation.log(p1 / p2)
            + (1 - p1) * Foundation.log((1 - p1) / (1 - p2))
        let maybeKL = try lhs.klDivergence(relativeTo: rhs)
        let kl = try #require(maybeKL)
        #expect(abs(kl - expected) <= 1e-12)
    }
}
