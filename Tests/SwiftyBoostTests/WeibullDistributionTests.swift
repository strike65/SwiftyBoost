import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Weibull distribution (Double)")
struct WeibullDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "shape", value: 0)
        ) {
            _ = try Distribution.Weibull<Double>(shape: 0, scale: 1)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.Weibull<Double>(shape: 1, scale: 0)
        }
    }

    @Test("Moments match gamma-based formula")
    func analyticMoments() throws {
        let shape: Double = 2.3
        let scale: Double = 1.7
        let dist = try Distribution.Weibull<Double>(shape: shape, scale: scale)

        let gamma1 = Foundation.tgamma(1 + 1 / shape)
        let gamma2 = Foundation.tgamma(1 + 2 / shape)

        let meanExpected = scale * gamma1
        let varianceExpected = scale * scale * (gamma2 - gamma1 * gamma1)

        if let mean = dist.mean {
            #expect(abs(mean - meanExpected) <= 1e-12)
        } else {
            Issue.record("Mean should be available for Weibull distribution")
        }

        if let variance = dist.variance {
            #expect(abs(variance - varianceExpected) <= 1e-12)
        } else {
            Issue.record("Variance should be available for Weibull distribution")
        }

        if let mode = dist.mode {
            let modeExpected = scale * Foundation.pow((shape - 1) / shape, 1 / shape)
            #expect(abs(mode - modeExpected) <= 1e-12)
        } else {
            Issue.record("Mode should exist for shape > 1")
        }
    }

    @Test("CDF and quantile round trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Weibull<Double>(shape: 0.9, scale: 1.2)
        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-8]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 5e-11, "Round-trip mismatch p=\(p)")
        }
    }

    @Test("Hazard equals pdf/sf inside support")
    func hazardMatchesRatio() throws {
        let dist = try Distribution.Weibull<Double>(shape: 1.5, scale: 0.75)
        let xs: [Double] = [0, 0.3, 1, 2, 4]
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
