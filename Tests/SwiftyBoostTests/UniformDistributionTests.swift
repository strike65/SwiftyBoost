import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Uniform distribution (Double)")
struct UniformDistributionTests {

    @Test("Parameter validation")
    func parameterValidation() {
        #expect(
            throws: DistributionError<Double>.parameterOutOfRange(
                name: "upper",
                min: 1,
                max: .infinity
            )
        ) {
            _ = try Distribution.Uniform<Double>(lower: 1, upper: 1)
        }
    }

    @Test("PDF and CDF match closed form")
    func pdfCdfMatch() throws {
        let lower: Double = -0.8
        let upper: Double = 3.2
        let dist = try Distribution.Uniform<Double>(lower: lower, upper: upper)

        let width = upper - lower
        let density = 1.0 / width

        let xs: [Double] = [lower - 1, lower, (lower + upper) / 2, upper, upper + 1]
        for x in xs {
            let pdf = try dist.pdf(x)
            if x < lower || x > upper {
                #expect(pdf == 0)
            } else {
                #expect(abs(pdf - density) <= 1e-12)
            }

            let cdf = try dist.cdf(x)
            if x <= lower {
                #expect(cdf == 0)
            } else if x >= upper {
                #expect(cdf == 1)
            } else {
                let expected = (x - lower) / width
                #expect(abs(cdf - expected) <= 1e-12)
            }
        }
    }

    @Test("Quantile round trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Uniform<Double>(lower: -5, upper: 7)
        let ps: [Double] = [0, 1e-6, 0.1, 0.5, 0.9, 1 - 1e-6, 1]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-12)
        }
    }
}
