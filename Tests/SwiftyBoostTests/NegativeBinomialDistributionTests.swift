import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Negative binomial distribution (Double)")
struct NegativeBinomialDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "successes", value: 0)
        ) {
            _ = try Distribution.NegativeBinomial<Double>(successes: 0, probabilityOfSuccess: 0.5)
        }
        #expect(
            throws: DistributionError<Double>.parameterOutOfRange(name: "probabilityOfSuccess", min: 0, max: 1)
        ) {
            _ = try Distribution.NegativeBinomial<Double>(successes: 2, probabilityOfSuccess: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterOutOfRange(name: "probabilityOfSuccess", min: 0, max: 1)
        ) {
            _ = try Distribution.NegativeBinomial<Double>(successes: 2, probabilityOfSuccess: 1.5)
        }
    }

    @Test("CDF/quantile round-trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.NegativeBinomial<Double>(successes: 3.5, probabilityOfSuccess: 0.35)
        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            var k = Foundation.floor(x)
            if k < 0 { k = 0 }
            var c = try dist.cdf(k)
            var iterations = 0
            while c < p - 1e-12 && iterations < 1_000_000 {
                k += 1
                c = try dist.cdf(k)
                iterations += 1
            }
            #expect(c >= p - 1e-10)
            if k > 0 {
                let lower = try dist.cdf(k - 1)
                #expect(lower <= p + 1e-10)
            }
        }

        let qs: [Double] = [1e-5, 1e-3, 0.05]
        for q in qs {
            let x = try dist.quantileComplement(q)
            var k = Foundation.floor(x)
            if k < 0 { k = 0 }
            var s = try dist.sf(k)
            var iterations = 0
            while s < q - 1e-12 && k > 0 && iterations < 1_000_000 {
                k -= 1
                s = try dist.sf(k)
                iterations += 1
            }
            #expect(s >= q - 1e-10)
            let upper = try dist.sf(k + 1)
            #expect(upper <= q + 1e-10)
        }
    }

    @Test("Analytic moments")
    func analyticMoments() throws {
        let successes: Double = 4.5
        let p: Double = 0.42
        let dist = try Distribution.NegativeBinomial<Double>(successes: successes, probabilityOfSuccess: p)

        if let mean = dist.mean {
            let expectedMean = successes * (1 - p) / p
            #expect(abs(mean - expectedMean) <= 1e-12)
        } else {
            Issue.record("Expected negative binomial mean to be available")
        }
        if let variance = dist.variance {
            let expectedVar = successes * (1 - p) / (p * p)
            #expect(abs(variance - expectedVar) <= 1e-12)
        } else {
            Issue.record("Expected negative binomial variance to be available")
        }
        if let skewness = dist.skewness {
            let expectedSkew = (2 - p) / Foundation.sqrt(successes * (1 - p))
            #expect(abs(skewness - expectedSkew) <= 1e-12)
        } else {
            Issue.record("Expected negative binomial skewness to be available")
        }
        if let kurtosis = dist.kurtosis {
            let expectedKurtosis = 3 + (6 / successes) + ((p * p) / (successes * (1 - p)))
            #expect(abs(kurtosis - expectedKurtosis) <= 1e-12)
        } else {
            Issue.record("Expected negative binomial kurtosis to be available")
        }
        if let excess = dist.kurtosisExcess {
            let expectedExcess = (6 - p * (6 - p)) / (successes * (1 - p))
            #expect(abs(excess - expectedExcess) <= 1e-12)
        } else {
            Issue.record("Expected negative binomial kurtosisExcess to be available")
        }
    }

    @Test("Hazard equals pdf/sf where defined")
    func hazardMatchesRatio() throws {
        let dist = try Distribution.NegativeBinomial<Double>(successes: 2.5, probabilityOfSuccess: 0.3)
        let xs: [Double] = [0, 1, 2, 5, 10]
        for x in xs {
            let pdf = try dist.pdf(Double(x))
            let sf = try dist.sf(Double(x))
            let hazard = try dist.hazard(Double(x))
            if sf > 0 {
                #expect(abs(hazard - (pdf / sf)) <= 1e-12)
            }
        }
    }

    @Test("Entropy matches analytic summation and handles degeneracies")
    func entropyMatchesAnalyticValues() throws {
        func analyticEntropy(successes r: Double, probability p: Double) throws -> Double {
            if p == 1 {
                return 0
            }
            let q = 1 - p
            let logGammaR = try SpecialFunctions.logGamma(r)
            let logP0 = r * Foundation.log(p)

            var entropy = 0.0
            var correction = 0.0
            var cumulative = 0.0

            let tailTolerance = 1e-12
            let termTolerance = 1e-18
            let maxIterations = 1_000_000
            var k = 0

            while k < maxIterations {
                let logGammaRK = try SpecialFunctions.logGamma(r + Double(k))
                let logGammaK1 = try SpecialFunctions.logGamma(Double(k + 1))
                let logProbability = logGammaRK - logGammaR - logGammaK1
                    + logP0 + Double(k) * Foundation.log(q)
                let probability: Double = Foundation.exp(logProbability)
                if probability <= 0 || !probability.isFinite {
                    break
                }
                let term: Double = -probability * Foundation.log(probability)
                let sum = entropy + term
                if abs(entropy) >= abs(term) {
                    correction += (entropy - sum) + term
                } else {
                    correction += (term - sum) + entropy
                }
                entropy = sum
                cumulative = min(1.0, cumulative + probability)
                let tail = max(0.0, 1.0 - cumulative)
                if tail <= tailTolerance || probability <= termTolerance && tail <= 1e-9 {
                    break
                }
                k += 1
            }

            return entropy + correction
        }

        let parameterSets: [(successes: Double, probability: Double, tolerance: Double)] = [
            (successes: 2.75, probability: 0.35, tolerance: 5e-12),
            (successes: 5.5, probability: 0.8, tolerance: 5e-12),
            (successes: 1.25, probability: 0.62, tolerance: 5e-12),
            (successes: 3.0, probability: 0.999, tolerance: 5e-11),
            (successes: 7.0, probability: 0.3, tolerance: 5e-12)
        ]

        for (successes, probability, tolerance) in parameterSets {
            let dist = try Distribution.NegativeBinomial<Double>(
                successes: successes,
                probabilityOfSuccess: probability
            )
            guard let entropy = dist.entropy else {
                Issue.record("Expected entropy to be available for NegativeBinomial(r: \(successes), p: \(probability)).")
                continue
            }
            let expected = try analyticEntropy(successes: successes, probability: probability)
            #expect(abs(entropy - expected) <= tolerance)
        }

        let degenerate = try Distribution.NegativeBinomial<Double>(successes: 4.0, probabilityOfSuccess: 1)
        if let entropy = degenerate.entropy {
            #expect(abs(entropy) <= 1e-12)
        } else {
            Issue.record("Expected entropy for degenerate NegativeBinomial with p = 1.")
        }
    }
}
