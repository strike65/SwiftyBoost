import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Inverse Chi-Squared distribution (Double)")
struct InverseChiSquaredDistributionTests {

    @Test("Typed wrapper matches dynamic backend")
    func typedMatchesDynamic() throws {
        let df: Double = 7.5
        let scale: Double = 0.35
        let typed = try Distribution.InverseChiSquared<Double>(degreesOfFreedom: df, scale: scale)
        let dynamic = try Distribution.Dynamic<Double>(
            distributionName: "inverse_chi_squared",
            parameters: ["df": df, "scale": scale]
        )

        let xs: [Double] = [0.05, 0.2, 0.5, 1.0, 3.0]
        for x in xs {
            let pt = try typed.pdf(x)
            let pd = try dynamic.pdf(x)
            #expect(abs(pt - pd) <= 1e-13)

            let ct = try typed.cdf(x)
            let cd = try dynamic.cdf(x)
            #expect(abs(ct - cd) <= 1e-13)

            let st = try typed.sf(x)
            let sd = try dynamic.sf(x)
            #expect(abs(st - sd) <= 1e-13)
        }
    }

    @Test("Quantile/CDF round-trip")
    func quantileRoundTrip() throws {
        let df: Double = 9.5
        let typed = try Distribution.InverseChiSquared<Double>(degreesOfFreedom: df)

        let ps: [Double] = [1e-6, 1e-4, 0.01, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try typed.quantile(p)
            let c = try typed.cdf(x)
            #expect(abs(c - p) <= 5e-11, "Round-trip mismatch for p = \(p)")
        }
    }

    @Test("Analytic moments where defined")
    func analyticMoments() throws {
        let df: Double = 10
        let scale: Double = 0.75
        let typed = try Distribution.InverseChiSquared<Double>(degreesOfFreedom: df, scale: scale)

        // Mean exists for ν > 2.
        let expectedMean = scale * df / (df - 2)
        if let mean = typed.mean {
            #expect(abs(mean - expectedMean) <= 1e-12)
        } else {
            Issue.record("Expected mean to exist for df = \(df)")
        }

        // Variance exists for ν > 4.
        let expectedVariance = (2 * pow(scale, 2) * pow(df, 2)) / (pow(df - 2, 2) * (df - 4))
        if let variance = typed.variance {
            #expect(abs(variance - expectedVariance) <= 1e-11)
        } else {
            Issue.record("Expected variance to exist for df = \(df)")
        }

        // Hazard should agree with PDF / SF.
        let x: Double = 0.6
        let hazard = try typed.hazard(x)
        let pdf = try typed.pdf(x)
        let sf = try typed.sf(x)
        #expect(abs(hazard - pdf / sf) <= 1e-12)
    }

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "degreesOfFreedom", value: -1)) {
            _ = try Distribution.InverseChiSquared<Double>(degreesOfFreedom: -1)
        }
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "degreesOfFreedom", value: 0)) {
            _ = try Distribution.InverseChiSquared<Double>(degreesOfFreedom: 0)
        }
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -0.1)) {
            _ = try Distribution.InverseChiSquared<Double>(degreesOfFreedom: 5, scale: -0.1)
        }
    }
}
