import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Log-normal distribution (Double)")
struct LogNormalDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.LogNormal<Double>(location: 0, scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -0.5)
        ) {
            _ = try Distribution.LogNormal<Double>(location: 0, scale: -0.5)
        }
    }

    @Test("Support is strictly positive")
    func supportBounds() throws {
        let dist = try Distribution.LogNormal<Double>(location: -0.5, scale: 0.9)
        #expect(dist.supportLowerBound >= 0)
        let xNegative = -0.1
        let c = try dist.cdf(xNegative)
        #expect(c.isNaN)
    }

    @Test("CDF/quantile round-trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.LogNormal<Double>(location: 0.3, scale: 0.65)
        let ps: [Double] = [1e-12, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round-trip mismatch for p=\(p)")
        }

        let qs: [Double] = [1e-8, 1e-4, 0.05]
        for q in qs {
            let x = try dist.quantileComplement(q)
            let s = try dist.sf(x)
            #expect(abs(s - q) <= 1e-10)
        }
    }

    @Test("Analytic moments and entropy")
    func analyticMoments() throws {
        let mu: Double = 0.7
        let sigma: Double = 0.8
        let dist = try Distribution.LogNormal<Double>(location: mu, scale: sigma)

        let expSigmaSq = Foundation.exp(sigma * sigma)
        let expectedMean = Foundation.exp(mu + (sigma * sigma) / 2)
        let expectedVar = (expSigmaSq - 1) * Foundation.exp(2 * mu + sigma * sigma)
        let expectedMode = Foundation.exp(mu - sigma * sigma)
        let expectedMedian = Foundation.exp(mu)
        let expectedSkew = (expSigmaSq + 2) * Foundation.sqrt(expSigmaSq - 1)

        if let mean = dist.mean {
            #expect(abs(mean - expectedMean) <= 1e-12)
        } else {
            Issue.record("Expected log-normal mean to be available")
        }
        if let variance = dist.variance {
            #expect(abs(variance - expectedVar) <= 1e-12)
        } else {
            Issue.record("Expected log-normal variance to be available")
        }
        #expect(abs((dist.mode ?? .nan) - expectedMode) <= 1e-12)
        #expect(abs(dist.median - expectedMedian) <= 1e-12)
        if let skew = dist.skewness {
            #expect(abs(skew - expectedSkew) <= 1e-12)
        } else {
            Issue.record("Expected log-normal skewness to be available")
        }
        let expectedKurtosisExcess = Foundation.exp(4 * sigma * sigma) + 2 * Foundation.exp(3 * sigma * sigma) + 3 * Foundation.exp(2 * sigma * sigma) - 6
        if let kurtosis = dist.kurtosis {
            let expectedKurtosis = expectedKurtosisExcess + 3
            #expect(abs(kurtosis - expectedKurtosis) <= 1e-12)
        } else {
            Issue.record("Expected log-normal kurtosis to be available")
        }
        if let excess = dist.kurtosisExcess {
            #expect(abs(excess - expectedKurtosisExcess) <= 1e-12)
        } else {
            Issue.record("Expected log-normal kurtosisExcess to be available")
        }

        if let entropy = dist.entropy {
            let expectedEntropy = mu + 0.5 + Foundation.log(sigma * Foundation.sqrt(2 * Double.pi))
            #expect(abs(entropy - expectedEntropy) <= 1e-12)
        } else {
            Issue.record("Expected log-normal entropy to be available")
        }
    }
}
