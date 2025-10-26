import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Hyperexponential distribution (Double)")
struct HyperexponentialDistributionTests {

    @Test("PDF/CDF agree with analytical mixture")
    func pdfCdfMatchAnalytical() throws {
        let rawProbabilities: [Double] = [0.2, 0.5, 0.8]
        let rates: [Double] = [0.5, 1.75, 4.2]
        let dist = try Distribution.Hyperexponential<Double>(
            probabilities: rawProbabilities,
            rates: rates
        )

        let probabilities = dist.probabilities
        let xs: [Double] = [0.0, 0.1, 0.25, 1.0, 2.5]
        for x in xs {
            let pdf = try dist.pdf(x)
            let expectedPdf = zip(probabilities, rates).reduce(0.0) { acc, pair in
                let (p, lambda) = pair
                return acc + p * lambda * Foundation.exp(-lambda * x)
            }
            #expect(abs(pdf - expectedPdf) <= 1e-12)

            let cdf = try dist.cdf(x)
            let expectedCdf = 1.0 - zip(probabilities, rates).reduce(0.0) { acc, pair in
                let (p, lambda) = pair
                return acc + p * Foundation.exp(-lambda * x)
            }
            #expect(abs(cdf - expectedCdf) <= 1e-12)
        }
    }

    @Test("Quantile/CDF round-trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Hyperexponential<Double>(
            probabilities: [0.1, 0.4, 0.3, 0.6],
            rates: [0.75, 1.25, 2.5, 3.5]
        )

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 5e-11, "Round-trip mismatch for p=\(p)")
        }
    }

    @Test("Moments and hazard derived values")
    func momentsAndHazard() throws {
        let dist = try Distribution.Hyperexponential<Double>(
            probabilities: [1.0, 3.0],
            rates: [0.8, 2.7]
        )

        let probabilities = dist.probabilities
        let rates = dist.rates

        let s1 = zip(probabilities, rates).reduce(0.0) { partial, pair in
            let (p, lambda) = pair
            return partial + p / lambda
        }
        let s2 = zip(probabilities, rates).reduce(0.0) { partial, pair in
            let (p, lambda) = pair
            let denom = lambda * lambda
            return partial + p / denom
        }
        let s3 = zip(probabilities, rates).reduce(0.0) { partial, pair in
            let (p, lambda) = pair
            return partial + p / Foundation.pow(lambda, 3)
        }
        let s4 = zip(probabilities, rates).reduce(0.0) { partial, pair in
            let (p, lambda) = pair
            return partial + p / Foundation.pow(lambda, 4)
        }

        let meanExpected = s1
        let varianceExpected = 2 * s2 - s1 * s1
        let skewnessExpected = {
            let numerator = 6 * s3 - (3 * (2 * s2 - s1 * s1) + s1 * s1) * s1
            let denominator = 2 * s2 - s1 * s1
            return numerator / Foundation.pow(denominator, 1.5)
        }()
        let kurtosisExpected = {
            let s1s1 = s1 * s1
            let numerator = 24 * s4 - 24 * s3 * s1 + 3 * (2 * (2 * s2 - s1s1) + s1s1) * s1s1
            let denominator = 2 * s2 - s1s1
            return numerator / (denominator * denominator)
        }()

        #expect(abs((dist.mean ?? .nan) - meanExpected) <= 1e-12)
        #expect(abs((dist.variance ?? .nan) - varianceExpected) <= 1e-12)
        #expect(abs((dist.skewness ?? .nan) - skewnessExpected) <= 1e-11)
        #expect(abs((dist.kurtosis ?? .nan) - kurtosisExpected) <= 1e-11)
        #expect(abs((dist.kurtosisExcess ?? .nan) - (kurtosisExpected - 3)) <= 1e-11)
        #expect((dist.mode ?? .nan) == 0.0)

        let x: Double = 0.85
        let hazard = try dist.hazard(x)
        let pdf = try dist.pdf(x)
        let sf = try dist.sf(x)
        #expect(abs(hazard - pdf / sf) <= 1e-12)

        let chf = try dist.chf(x)
        #expect(abs(chf + Foundation.log(sf)) <= 1e-12)
    }

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Hyperexponential<Double>(rates: [])
        }
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "rate[0]", value: 0)) {
            _ = try Distribution.Hyperexponential<Double>(rates: [0])
        }
        #expect(throws: DistributionError<Double>.invalidCombination(message: "Probabilities must match the number of rates", value: nil)) {
            _ = try Distribution.Hyperexponential<Double>(probabilities: [1, 2], rates: [1])
        }
        #expect(throws: DistributionError<Double>.invalidCombination(message: "Probabilities must contain at least one positive entry", value: nil)) {
            _ = try Distribution.Hyperexponential<Double>(probabilities: [0, 0], rates: [1, 2])
        }
    }
}
