import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Kolmogorov–Smirnov distribution (Double)")
struct KolmogorovSmirnovDistributionTests {

    @Test("Support bounds and hazard relationship")
    func supportAndHazard() throws {
        let n: Double = 25
        let dist = try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: n)
        #expect(dist.supportLowerBound == 0)
        #expect(dist.supportUpperBound.isFinite)

        let xs: [Double] = [0.05, 0.1, 0.2, 0.5]
        for x in xs {
            let pdf = try dist.pdf(x)
            let sf = try dist.sf(x)
            let haz = try dist.hazard(x)
            if sf > 0 {
                #expect(abs(haz - (pdf / sf)) <= 1e-12)
            }
        }
    }

    @Test("CDF/quantile round-trip")
    func quantileRoundTrip() throws {
        let n: Double = 40
        let dist = try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: n)
        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round-trip mismatch for p=\(p)")
        }

        let qs: [Double] = [1e-8, 1e-5, 1e-3, 0.1]
        for q in qs {
            let x = try dist.quantileComplement(q)
            let s = try dist.sf(x)
            #expect(abs(s - q) <= 1e-10, "Round-trip mismatch for q=\(q)")
        }
    }

    @Test("Analytic moments match Boost formulas")
    func analyticMoments() throws {
        let n: Double = 17
        let dist = try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: n)

        let rootHalfPi = Foundation.sqrt(Double.pi / 2)
        let ln2 = Foundation.log(2.0)
        let expectedMean = rootHalfPi * ln2 / Foundation.sqrt(n)
        #expect(abs((dist.mean ?? .nan) - expectedMean) <= 1e-12)

        let piSqOver6 = (Double.pi * Double.pi) / 6.0
        let expectedVariance = (piSqOver6 - Double.pi * ln2 * ln2) / (2 * n)
        #expect(abs((dist.variance ?? .nan) - expectedVariance) <= 1e-12)

        let zeta3 = Constants.zetaThree(Double.self)
        let ex3 = 0.5625 * rootHalfPi * zeta3 / (n * Foundation.sqrt(n))
        let mean = expectedMean
        let variance = expectedVariance
        let expectedSkew = (ex3 - 3 * mean * variance - mean * mean * mean) / (variance * Foundation.sqrt(variance))
        if let skew = dist.skewness {
            #expect(abs(skew - expectedSkew) <= 1e-12)
        } else {
            Issue.record("Expected Kolmogorov–Smirnov skewness to be provided by Boost")
        }

        let ex4 = 7 * piSqOver6 * piSqOver6 / (20 * n * n)
        let expectedKurtosis = (ex4
            - 4 * mean * expectedSkew * variance * Foundation.sqrt(variance)
            - 6 * mean * mean * variance
            - mean * mean * mean * mean) / (variance * variance)
        if let kurtosis = dist.kurtosis {
            #expect(abs(kurtosis - expectedKurtosis) <= 1e-12)
        } else {
            Issue.record("Expected Kolmogorov–Smirnov kurtosis to be provided by Boost")
        }
        if let kurtosisExcess = dist.kurtosisExcess {
            #expect(abs(kurtosisExcess - (expectedKurtosis - 3)) <= 1e-12)
        } else {
            Issue.record("Expected Kolmogorov–Smirnov kurtosisExcess to be provided by Boost")
        }
    }

    @Test("Median aligns with 0.5 quantile")
    func medianMatchesQuantile() throws {
        let dist = try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: 8)
        let median = dist.median
        let quantile = try dist.quantile(0.5)
        #expect(abs(median - quantile) <= 1e-12)
    }
}
