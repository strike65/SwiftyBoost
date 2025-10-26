import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Pareto distribution (Double)")
struct ParetoDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -1)) {
            _ = try Distribution.Pareto<Double>(scale: -1, shape: 2)
        }
        #expect(throws: DistributionError<Double>.parameterNotPositive(name: "shape", value: 0)) {
            _ = try Distribution.Pareto<Double>(scale: 1, shape: 0)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Pareto<Double>(scale: .nan, shape: 2)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Pareto<Double>(scale: 1, shape: .infinity)
        }
    }

    @Test("Analytic mean and variance")
    func analyticMoments() throws {
        let scale: Double = 1.5
        let shape: Double = 3.25
        let pareto = try Distribution.Pareto<Double>(scale: scale, shape: shape)

        let expectedMean = (shape * scale) / (shape - 1)
        let expectedVariance = (shape * scale * scale) /
            ((shape - 1) * (shape - 1) * (shape - 2))

        if let mean = pareto.mean {
            #expect(abs(mean - expectedMean) <= 1e-12)
        } else {
            Issue.record("Expected Pareto mean to be available")
        }
        if let variance = pareto.variance {
            #expect(abs(variance - expectedVariance) <= 1e-12)
        } else {
            Issue.record("Expected Pareto variance to be available")
        }
        #expect(pareto.mode == scale)
    }

    @Test("Typed vs dynamic parity")
    func parityWithDynamic() throws {
        let scale: Double = 2.0
        let shape: Double = 4.5
        let typed = try Distribution.Pareto<Double>(scale: scale, shape: shape)
        let dynamic = try Distribution.Dynamic<Double>(
            distributionName: "pareto",
            parameters: ["scale": scale, "shape": shape]
        )

        let xs: [Double] = [2.0, 2.5, 3.0, 5.0, 10.0]
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

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let qt = try typed.quantile(p)
            let qd = try dynamic.quantile(p)
            #expect(abs(qt - qd) <= 1e-10)
        }
    }
}
