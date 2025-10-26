import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Non-central beta distribution")
struct NonCentralBetaDistributionTests {

    @Test("Skewness and kurtosis from pdf integration (λ > 0)")
    func highOrderMomentsFromIntegration() throws {
        let alpha: Double = 2.75
        let beta: Double = 3.5
        let lambda: Double = 4.2
        let dist = try Distribution.NonCentralBeta<Double>(alpha: alpha, beta: beta, lambda: lambda)

        func rawMoment(_ order: Int) throws -> Double {
            let result = try Quadrature.integrate(
                using: Quadrature.Rule.gaussKronrod(points: 61),
                over: .finite(lower: 0, upper: 1)
            ) { x -> Double in
                let clamped = min(max(x, Double.leastNonzeroMagnitude), 1 - Double.ulpOfOne)
                guard let density = try? dist.pdf(clamped) else { return 0 }
                return Foundation.pow(clamped, Double(order)) * density
            }
            return result.value
        }

        let m1 = try rawMoment(1)
        let m2 = try rawMoment(2)
        let m3 = try rawMoment(3)
        let m4 = try rawMoment(4)

        let central2 = m2 - m1 * m1
        let central3 = m3 - 3 * m1 * m2 + 2 * m1 * m1 * m1
        var central4 = m4 - 4 * m1 * m3
        central4 += 6 * m1 * m1 * m2
        central4 -= 3 * m1 * m1 * m1 * m1

        #expect(central2 > 0)

        let expectedSkew = central3 / Foundation.pow(central2, 1.5)
        let expectedKurtosis = central4 / (central2 * central2)
        let expectedExcess = expectedKurtosis - 3

        let skew = try #require(dist.skewness)
        let kurtosis = try #require(dist.kurtosis)
        let excess = try #require(dist.kurtosisExcess)

        #expect(abs(skew - expectedSkew) <= 5e-8)
        #expect(abs(kurtosis - expectedKurtosis) <= 5e-8)
        #expect(abs(excess - expectedExcess) <= 5e-8)
    }

    @Test("λ → 0 reduces to central beta analytic moments")
    func reducesToCentralBeta() throws {
        let alpha: Double = 4.2
        let beta: Double = 2.8
        let lambda: Double = 0
        let dist = try Distribution.NonCentralBeta<Double>(alpha: alpha, beta: beta, lambda: lambda)

        let expectedVariance = (alpha * beta) /
            Foundation.pow(alpha + beta, 2) /
            (alpha + beta + 1)
        let expectedSkew = 2 * (beta - alpha) *
            Foundation.sqrt(alpha + beta + 1) /
            ((alpha + beta + 2) * Foundation.sqrt(alpha * beta))
        let numerator = 6 * (Foundation.pow(alpha - beta, 2) * (alpha + beta + 1) - alpha * beta * (alpha + beta + 2))
        let denominator = alpha * beta * (alpha + beta + 2) * (alpha + beta + 3)
        let expectedExcess = numerator / denominator
        let expectedKurtosis = expectedExcess + 3

        if let variance = dist.variance {
            #expect(abs(variance - expectedVariance) <= 5e-13)
        } else {
            Issue.record("Expected variance for central beta")
        }

        let skew = try #require(dist.skewness)
        let kurtosis = try #require(dist.kurtosis)
        let excess = try #require(dist.kurtosisExcess)

        #expect(abs(skew - expectedSkew) <= 5e-12)
        #expect(abs(kurtosis - expectedKurtosis) <= 5e-12)
        #expect(abs(excess - expectedExcess) <= 5e-12)
    }
}
