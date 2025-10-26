import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Non-central F distribution (Double)")
struct NonCentralFDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "degreesOfFreedom1", value: 0)) {
            _ = try Distribution.NonCentralF<Double>(degreesOfFreedom1: 0, degreesOfFreedom2: 5, nonCentrality: 1)
        }
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "degreesOfFreedom2", value: 0)) {
            _ = try Distribution.NonCentralF<Double>(degreesOfFreedom1: 4, degreesOfFreedom2: 0, nonCentrality: 1)
        }
        #expect(throws: DistributionError<Double>.parameterOutOfRange(name: "nonCentrality", min: 0, max: .infinity)) {
            _ = try Distribution.NonCentralF<Double>(degreesOfFreedom1: 4, degreesOfFreedom2: 5, nonCentrality: -0.5)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.NonCentralF<Double>(degreesOfFreedom1: .nan, degreesOfFreedom2: 5, nonCentrality: 1)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.NonCentralF<Double>(degreesOfFreedom1: 4, degreesOfFreedom2: 5, nonCentrality: .infinity)
        }
    }

    @Test("Analytic mean and variance")
    func analyticMoments() throws {
        let df1: Double = 6.0
        let df2: Double = 8.5
        let lambda: Double = 3.0
        let dist = try Distribution.NonCentralF<Double>(
            degreesOfFreedom1: df1,
            degreesOfFreedom2: df2,
            nonCentrality: lambda
        )

        let expectedMean = (df2 * (df1 + lambda)) / (df1 * (df2 - 2))
        let numeratorVariance = 2 * pow(df2 / df1, 2) *
            (pow(df1 + lambda, 2) + (df1 + 2 * lambda) * (df2 - 2))
        let expectedVariance = numeratorVariance /
            (pow(df2 - 2, 2) * (df2 - 4))

        if let mean = dist.mean {
            #expect(abs(mean - expectedMean) <= 1e-12)
        } else {
            Issue.record("Expected mean to be available for non-central F")
        }
        if let variance = dist.variance {
            #expect(abs(variance - expectedVariance) <= 1e-12)
        } else {
            Issue.record("Expected variance to be available for non-central F")
        }
    }

    @Test("Typed vs dynamic parity")
    func parityWithDynamic() throws {
        let df1: Double = 5.5
        let df2: Double = 7.0
        let lambda: Double = 2.25
        let typed = try Distribution.NonCentralF<Double>(
            degreesOfFreedom1: df1,
            degreesOfFreedom2: df2,
            nonCentrality: lambda
        )
        let dynamic = try Distribution.Dynamic<Double>(
            distributionName: "non_central_f",
            parameters: ["df1": df1, "df2": df2, "lambda": lambda]
        )

        let xs: [Double] = [0.05, 0.3, 0.75, 1.2, 2.5, 4.0]
        for x in xs {
            let pdfTyped = try typed.pdf(x)
            let pdfDynamic = try dynamic.pdf(x)
            #expect(abs(pdfTyped - pdfDynamic) <= 1e-12)

            let cdfTyped = try typed.cdf(x)
            let cdfDynamic = try dynamic.cdf(x)
            #expect(abs(cdfTyped - cdfDynamic) <= 1e-12)

            let sfTyped = try typed.sf(x)
            let sfDynamic = try dynamic.sf(x)
            #expect(abs(sfTyped - sfDynamic) <= 1e-12)
        }

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let qt = try typed.quantile(p)
            let qd = try dynamic.quantile(p)
            #expect(abs(qt - qd) <= 1e-10)
        }

        if let meanTyped = typed.mean, let meanDynamic = dynamic.mean {
            #expect(abs(meanTyped - meanDynamic) <= 1e-12)
        }
        if let varianceTyped = typed.variance, let varianceDynamic = dynamic.variance {
            #expect(abs(varianceTyped - varianceDynamic) <= 1e-12)
        }
    }

    @Test("Lambda -> 0 reduces to central FisherF")
    func reducesToCentralCase() throws {
        let df1: Double = 4
        let df2: Double = 9
        let nonCentral = try Distribution.NonCentralF<Double>(
            degreesOfFreedom1: df1,
            degreesOfFreedom2: df2,
            nonCentrality: 0
        )
        let central = try Distribution.FisherF<Double>(
            degreesOfFreedom1: df1,
            degreesOfFreedom2: df2
        )

        let xs: [Double] = [0.1, 0.5, 1.0, 2.0]
        let tolerance = 1e-10
        for x in xs {
            #expect(abs((try nonCentral.pdf(x)) - (try central.pdf(x))) <= tolerance)
            #expect(abs((try nonCentral.cdf(x)) - (try central.cdf(x))) <= tolerance)
        }

        if let mean = nonCentral.mean, let centralMean = central.mean {
            #expect(abs(mean - centralMean) <= 1e-12)
        }
        if let variance = nonCentral.variance, let centralVariance = central.variance {
            #expect(abs(variance - centralVariance) <= 1e-12)
        }
    }
}
