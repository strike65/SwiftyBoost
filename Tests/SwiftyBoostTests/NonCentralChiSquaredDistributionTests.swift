import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Non-central chi-squared distribution (Double)")
struct NonCentralChiSquaredDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "degreesOfFreedom", value: 0)
        ) {
            _ = try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: 0, nonCentrality: 1)
        }
        #expect(
            throws: DistributionError<Double>.parameterOutOfRange(name: "nonCentrality", min: 0, max: .infinity)
        ) {
            _ = try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: 3, nonCentrality: -0.1)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: .nan, nonCentrality: 0.5)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: 3, nonCentrality: .infinity)
        }
    }

    @Test("Analytic moments align with cumulant formulas")
    func analyticMoments() throws {
        let df: Double = 6.5
        let lambda: Double = 2.75
        let dist = try Distribution.NonCentralChiSquared<Double>(
            degreesOfFreedom: df,
            nonCentrality: lambda
        )

        let expectedMean = df + lambda
        let variance = 2 * (df + 2 * lambda)
        let mu3 = 8 * (df + 3 * lambda)
        let base = df + 2 * lambda
        let mu4 = 48 * (df + 4 * lambda) + 12 * base * base

        let sigma = Foundation.sqrt(variance)
        let expectedSkewness = mu3 / (variance * sigma)
        let expectedKurtosis = mu4 / (variance * variance)
        let expectedExcess = expectedKurtosis - 3

        if let mean = dist.mean {
            #expect(abs(mean - expectedMean) <= 1e-12)
        } else {
            Issue.record("Expected mean to be available for non-central chi-squared")
        }
        if let varianceValue = dist.variance {
            #expect(abs(varianceValue - variance) <= 1e-12)
        } else {
            Issue.record("Expected variance to be available for non-central chi-squared")
        }
        if let skew = dist.skewness {
            #expect(abs(skew - expectedSkewness) <= 1e-12)
        } else {
            Issue.record("Expected skewness to be available for non-central chi-squared")
        }
        if let kurtosis = dist.kurtosis {
            #expect(abs(kurtosis - expectedKurtosis) <= 1e-12)
        } else {
            Issue.record("Expected kurtosis to be available for non-central chi-squared")
        }
        if let excess = dist.kurtosisExcess {
            #expect(abs(excess - expectedExcess) <= 1e-12)
        } else {
            Issue.record("Expected kurtosisExcess to be available for non-central chi-squared")
        }
    }

    @Test("Typed vs dynamic parity")
    func parityWithDynamic() throws {
        let df: Double = 5.5
        let lambda: Double = 1.75
        let typed = try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: df, nonCentrality: lambda)
        let dynamic = try Distribution.Dynamic<Double>(
            distributionName: "non_central_chi_squared",
            parameters: ["df": df, "lambda": lambda]
        )

        let xs: [Double] = [0.1, 0.5, 1.0, 3.0, 6.0, 10.0]
        for x in xs {
            let pdfTyped = try typed.pdf(x)
            let pdfDyn = try dynamic.pdf(x)
            #expect(abs(pdfTyped - pdfDyn) <= 1e-12)

            let cdfTyped = try typed.cdf(x)
            let cdfDyn = try dynamic.cdf(x)
            #expect(abs(cdfTyped - cdfDyn) <= 1e-12)

            let sfTyped = try typed.sf(x)
            let sfDyn = try dynamic.sf(x)
            #expect(abs(sfTyped - sfDyn) <= 1e-12)
        }

        let probabilities: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in probabilities {
            let qt = try typed.quantile(p)
            let qd = try dynamic.quantile(p)
            #expect(abs(qt - qd) <= 1e-10)

            let ct = try typed.cdf(qt)
            #expect(ct >= p - 1e-12)
        }

        let complementProbs: [Double] = [1e-6, 1e-3, 0.1, 0.5]
        for q in complementProbs {
            let qt = try typed.quantileComplement(q)
            let qd = try dynamic.quantileComplement(q)
            #expect(abs(qt - qd) <= 1e-10)

            let st = try typed.sf(qt)
            #expect(abs(st - q) <= 1e-12)
        }
    }

    @Test("Lambda -> 0 reduces to central chi-squared")
    func reducesToCentralCase() throws {
        let df: Double = 9
        let lambda: Double = 0
        let nonCentral = try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: df, nonCentrality: lambda)
        let central = try Distribution.ChiSquared<Double>(degreesOfFreedom: df)

        let x: Double = 3.5
        #expect(abs((try nonCentral.pdf(x)) - (try central.pdf(x))) <= 1e-12)
        #expect(abs((try nonCentral.cdf(x)) - (try central.cdf(x))) <= 1e-12)
        #expect(abs((try nonCentral.sf(x)) - (try central.sf(x))) <= 1e-12)

        if let mean = nonCentral.mean, let centralMean = central.mean {
            #expect(abs(mean - centralMean) <= 1e-12)
        }
        if let variance = nonCentral.variance, let centralVariance = central.variance {
            #expect(abs(variance - centralVariance) <= 1e-12)
        }
        if let skew = nonCentral.skewness, let centralSkew = central.skewness {
            #expect(abs(skew - centralSkew) <= 1e-12)
        }
        if let kurt = nonCentral.kurtosis, let centralKurt = central.kurtosis {
            #expect(abs(kurt - centralKurt) <= 1e-12)
        }
    }
}
