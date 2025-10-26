import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Poisson distribution (Double)")
struct PoissonDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(throws: DistributionError<Double>.parameterOutOfRange(name: "lambda", min: 0, max: .infinity)) {
            _ = try Distribution.Poisson<Double>(lambda: -0.1)
        }
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Poisson<Double>(lambda: .nan)
        }
    }

    @Test("Analytic moments")
    func analyticMoments() throws {
        let lambda: Double = 3.75
        let p = try Distribution.Poisson<Double>(lambda: lambda)

        #expect(abs((p.mean ?? .nan) - lambda) <= 1e-12)
        #expect(abs((p.variance ?? .nan) - lambda) <= 1e-12)

        if let skewness = p.skewness {
            #expect(abs(skewness - (1 / Foundation.sqrt(lambda))) <= 1e-12)
        } else {
            Issue.record("Expected skewness to be available for Poisson")
        }
        if let kurtosisExcess = p.kurtosisExcess {
            #expect(abs(kurtosisExcess - (1 / lambda)) <= 1e-12)
        } else {
            Issue.record("Expected kurtosisExcess to be available for Poisson")
        }
    }

    @Test("Typed vs dynamic parity")
    func parityWithDynamic() throws {
        let lambda: Double = 4.2
        let typed = try Distribution.Poisson<Double>(lambda: lambda)
        let dynamic = try Distribution.Dynamic<Double>(
            distributionName: "poisson",
            parameters: ["lambda": lambda]
        )

        let ks: [Double] = [0, 1, 2, 3, 5, 10]
        for k in ks {
            let pdfTyped = try typed.pdf(k)
            let pdfDynamic = try dynamic.pdf(k)
            #expect(abs(pdfTyped - pdfDynamic) <= 1e-12)

            let cdfTyped = try typed.cdf(k)
            let cdfDynamic = try dynamic.cdf(k)
            #expect(abs(cdfTyped - cdfDynamic) <= 1e-12)

            let sfTyped = try typed.sf(k)
            let sfDynamic = try dynamic.sf(k)
            #expect(abs(sfTyped - sfDynamic) <= 1e-12)
        }

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9]
        for p in ps {
            let qt = try typed.quantile(p)
            let qd = try dynamic.quantile(p)
            #expect(abs(qt - qd) <= 1e-10)
        }
    }
}
