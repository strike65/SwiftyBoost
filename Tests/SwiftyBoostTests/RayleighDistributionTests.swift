import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Rayleigh distribution (Double)")
struct RayleighDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.Rayleigh<Double>(scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -1)
        ) {
            _ = try Distribution.Rayleigh<Double>(scale: -1)
        }
    }

    @Test("Moments and entropy match analytic forms")
    func closedFormMoments() throws {
        let scale: Double = 1.75
        let dist = try Distribution.Rayleigh<Double>(scale: scale)

        let meanExpected = scale * Foundation.sqrt(Double.pi / 2)
        let varianceExpected = (2 - Double.pi / 2) * scale * scale
        let modeExpected = scale
        let medianExpected = scale * Foundation.sqrt(2 * Foundation.log(2))
        let skewnessExpected = (2 * Foundation.sqrt(Double.pi) * (Double.pi - 3)) / Foundation.pow(4 - Double.pi, 1.5)
        let kurtosisExcessExpected = (-6 * Double.pi * Double.pi + 24 * Double.pi - 16) / Foundation.pow(4 - Double.pi, 2)
        let eulerGamma = 0.5772156649015328606
        let entropyExpected = 1 + Foundation.log(scale / Foundation.sqrt(2)) + 0.5 * eulerGamma

        if let mean = dist.mean {
            #expect(abs(mean - meanExpected) <= 1e-12)
        } else {
            Issue.record("Mean should be available")
        }

        if let variance = dist.variance {
            #expect(abs(variance - varianceExpected) <= 1e-12)
        } else {
            Issue.record("Variance should be available")
        }

        #expect(dist.mode == modeExpected)
        #expect(abs(dist.median - medianExpected) <= 1e-12)
        #expect(abs((dist.skewness ?? .nan) - skewnessExpected) <= 1e-12)
        #expect(abs((dist.kurtosisExcess ?? .nan) - kurtosisExcessExpected) <= 1e-12)
        #expect(abs((dist.kurtosis ?? .nan) - (kurtosisExcessExpected + 3)) <= 1e-12)
        #expect(abs((dist.entropy ?? .nan) - entropyExpected) <= 1e-12)
    }

    @Test("CDF/quantile round trip")
    func cdfQuantileRoundTrip() throws {
        let dist = try Distribution.Rayleigh<Double>(scale: 0.42)
        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-8]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round trip mismatch for p=\(p)")
        }
    }

    @Test("Hazard equals pdf/sf in support")
    func hazardMatchesRatio() throws {
        let dist = try Distribution.Rayleigh<Double>(scale: 0.95)
        let xs: [Double] = [0, 0.5, 1, 2, 4]
        for x in xs {
            let pdf = try dist.pdf(x)
            let sf = try dist.sf(x)
            let hazard = try dist.hazard(x)
            if sf > 0 {
                #expect(abs(hazard - (pdf / sf)) <= 1e-12)
            }
        }
    }
}
