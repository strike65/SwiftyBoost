import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Non-central StudentT distribution (Double)")
struct NonCentralStudentTDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "degreesOfFreedom", value: 0)) {
            _ = try Distribution.NonCentralStudentT<Double>(degreesOfFreedom: 0, nonCentrality: 0.5)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.NonCentralStudentT<Double>(degreesOfFreedom: .nan, nonCentrality: 0.5)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.NonCentralStudentT<Double>(degreesOfFreedom: 5, nonCentrality: .infinity)
        }
    }

    @Test("Typed vs dynamic parity")
    func parityWithDynamic() throws {
        let df: Double = 7.5
        let lambda: Double = -1.25
        let typed = try Distribution.NonCentralStudentT<Double>(
            degreesOfFreedom: df,
            nonCentrality: lambda
        )
        let dynamic = try Distribution.Dynamic<Double>(
            distributionName: "non_central_t",
            parameters: ["df": df, "lambda": lambda]
        )

        let xs: [Double] = [-3.0, -1.0, 0.0, 1.0, 2.5]
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

    @Test("Lambda -> 0 reduces to central StudentT")
    func reducesToCentralCase() throws {
        let df: Double = 9
        let nonCentral = try Distribution.NonCentralStudentT<Double>(
            degreesOfFreedom: df,
            nonCentrality: 0
        )
        let central = try Distribution.StudentT<Double>(degreesOfFreedom: df)

        let xs: [Double] = [-2.0, -0.5, 0.0, 0.5, 2.0]
        for x in xs {
            #expect(abs((try nonCentral.pdf(x)) - (try central.pdf(x))) <= 1e-12)
            #expect(abs((try nonCentral.cdf(x)) - (try central.cdf(x))) <= 1e-12)
            #expect(abs((try nonCentral.sf(x)) - (try central.sf(x))) <= 1e-12)
        }

        if let mean = nonCentral.mean, let centralMean = central.mean {
            #expect(abs(mean - centralMean) <= 1e-12)
        }
        if let variance = nonCentral.variance, let centralVariance = central.variance {
            #expect(abs(variance - centralVariance) <= 1e-12)
        }
    }
}
