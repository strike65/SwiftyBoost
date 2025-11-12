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

    @Test("Typed continuous distributions support cross-family KL divergence")
    func typedCrossFamilyContinuous() throws {
        let gamma = try Distribution.Gamma<Double>(shape: 3.0, scale: 0.75)
        let beta = try Distribution.Beta<Double>(alpha: 2.5, beta: 4.0)
        let maybeTyped = try gamma.klDivergence(relativeTo: beta)
        let typedKL = try #require(maybeTyped)

        let gammaDynamic = try Distribution.Dynamic<Double>(
            distributionName: "gamma",
            parameters: ["shape": 3.0, "scale": 0.75]
        )
        let betaDynamic = try Distribution.Dynamic<Double>(
            distributionName: "beta",
            parameters: ["alpha": 2.5, "beta": 4.0]
        )
        let maybeDynamic = try gammaDynamic.klDivergence(relativeTo: betaDynamic)
        let dynamicKL = try #require(maybeDynamic)
        #expect(abs(typedKL - dynamicKL) <= 5e-9)
    }

    @Test("Typed discrete distributions support cross-family KL divergence")
    func typedCrossFamilyDiscrete() throws {
        let bernoulli = try Distribution.Bernoulli<Double>(p: 0.42)
        let binomial = try Distribution.Binomial<Double>(numberOfTrials: 1, probabibilityOfSuccess: 0.63)
        let maybeKL = try bernoulli.klDivergence(relativeTo: binomial)
        let kl = try #require(maybeKL)

        let p: Double = 0.42
        let q: Double = 0.63
        let expected =
            p * Foundation.log(p / q)
            + (1 - p) * Foundation.log((1 - p) / (1 - q))
        #expect(abs(kl - expected) <= 1e-12)
    }

    @Test("Continuous KL divergence honors integration bounds")
    func typedContinuousIntegrationBounds() throws {
        let lhs = try Distribution.Uniform<Double>(lower: 0, upper: 1)
        let rhs = try Distribution.Uniform<Double>(lower: 0, upper: 2)
        let options = Distribution.KLDivergenceOptions<Double>(
            integrationLowerBound: 0,
            integrationUpperBound: 0.25
        )
        let maybeKL = try lhs.klDivergence(relativeTo: rhs, options: options)
        let kl = try #require(maybeKL)
        #expect(abs(kl - (0.25 * Foundation.log(2))) <= 1e-12)
    }

    @Test("Discrete KL divergence honors integration bounds")
    func typedDiscreteIntegrationBounds() throws {
        let bernoulli = try Distribution.Bernoulli<Double>(p: 0.3)
        let binomial = try Distribution.Binomial<Double>(numberOfTrials: 1, probabibilityOfSuccess: 0.55)
        let options = Distribution.KLDivergenceOptions<Double>(
            integrationLowerBound: 0,
            integrationUpperBound: 0
        )
        let maybeKL = try bernoulli.klDivergence(relativeTo: binomial, options: options)
        let kl = try #require(maybeKL)
        let expected = (1 - 0.3) * Foundation.log((1 - 0.3) / (1 - 0.55))
        #expect(abs(kl - expected) <= 1e-12)
    }
}
