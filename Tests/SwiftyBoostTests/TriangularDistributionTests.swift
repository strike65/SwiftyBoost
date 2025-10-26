import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Triangular distribution (Double)")
struct TriangularDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(
            throws: DistributionError<Double>.parameterOutOfRange(
                name: "upper",
                min: 1,
                max: .infinity
            )
        ) {
            _ = try Distribution.Triangular<Double>(lower: 1, mode: 2, upper: 1)
        }

        #expect(
            throws: DistributionError<Double>.parameterOutOfRange(
                name: "mode",
                min: 0,
                max: 2
            )
        ) {
            _ = try Distribution.Triangular<Double>(lower: 0, mode: 3, upper: 2)
        }
    }

    @Test("Moments match analytic expressions")
    func analyticMoments() throws {
        let lower: Double = -1.5
        let mode: Double = 0.25
        let upper: Double = 2.0

        let dist = try Distribution.Triangular<Double>(lower: lower, mode: mode, upper: upper)

        let meanExpected = (lower + mode + upper) / 3
        let varianceExpected = (
            lower * lower + mode * mode + upper * upper
            - lower * mode - lower * upper - mode * upper
        ) / 18

        if let mean = dist.mean {
            #expect(abs(mean - meanExpected) <= 1e-12)
        } else {
            Issue.record("Mean should be available for triangular distribution")
        }

        if let variance = dist.variance {
            #expect(abs(variance - varianceExpected) <= 1e-12)
        } else {
            Issue.record("Variance should be available for triangular distribution")
        }

        #expect(dist.median >= lower && dist.median <= upper)
        #expect(dist.mode == mode)
    }

    @Test("CDF/quantile round trip and support endpoints")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Triangular<Double>(lower: -2, mode: -0.2, upper: 3.5)

        #expect(try dist.cdf(dist.lower) == 0)
        #expect(try dist.cdf(dist.upper) == 1)
        #expect(try dist.pdf(dist.lower - 10) == 0)
        #expect(try dist.pdf(dist.upper + 10) == 0)

        let ps: [Double] = [1e-8, 1e-5, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round trip mismatch for p=\(p)")
        }
    }
}
