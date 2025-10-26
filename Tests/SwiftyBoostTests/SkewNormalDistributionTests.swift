import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Skew-normal distribution (Double)")
struct SkewNormalDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.SkewNormal<Double>(location: 0, scale: 0, shape: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -0.5)
        ) {
            _ = try Distribution.SkewNormal<Double>(location: 0, scale: -0.5, shape: 1)
        }
    }

    @Test("Reduces to normal when shape is zero")
    func reducesToNormal() throws {
        let skew = try Distribution.SkewNormal<Double>(location: -0.5, scale: 2, shape: 0)
        let normal = try Distribution.Normal<Double>(mean: -0.5, sd: 2)

        let xs: [Double] = [-5, -2, -1, 0, 1, 3, 5]
        for x in xs {
            let skewPdf = try skew.pdf(x)
            let normalPdf = try normal.pdf(x)
            #expect(abs(skewPdf - normalPdf) <= 1e-12)

            let skewCdf = try skew.cdf(x)
            let normalCdf = try normal.cdf(x)
            #expect(abs(skewCdf - normalCdf) <= 1e-12)
        }

        #expect(abs((skew.mean ?? .nan) - (normal.mean ?? .nan)) <= 1e-12)
        #expect(abs((skew.variance ?? .nan) - (normal.variance ?? .nan)) <= 1e-12)
    }

    @Test("Skewness follows the sign of shape")
    func skewnessMatchesShape() throws {
        let positive = try Distribution.SkewNormal<Double>(location: 0, scale: 1, shape: 4)
        let negative = try Distribution.SkewNormal<Double>(location: 0, scale: 1, shape: -4)

        #expect((positive.skewness ?? .nan) > 0)
        #expect((negative.skewness ?? .nan) < 0)
        #expect(abs(positive.median + negative.median) <= 1e-12)
    }

    @Test("Quantile and CDF round trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.SkewNormal<Double>(location: 1.5, scale: 0.8, shape: 3.5)
        let ps: [Double] = [1e-9, 1e-5, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 5e-11)
        }
    }
}
