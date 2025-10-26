import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Laplace distribution (Double)")
struct LaplaceDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.Laplace<Double>(location: 0, scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -0.2)
        ) {
            _ = try Distribution.Laplace<Double>(location: 0, scale: -0.2)
        }
    }

    @Test("Closed-form moments and entropy")
    func closedFormMoments() throws {
        let location: Double = -1.2
        let scale: Double = 0.8
        let dist = try Distribution.Laplace<Double>(location: location, scale: scale)

        #expect(dist.mean == location)
        #expect(abs((dist.variance ?? .nan) - (2 * scale * scale)) <= 1e-12)
        #expect(dist.mode == location)
        #expect(dist.median == location)
        #expect((dist.skewness ?? .nan) == 0)
        #expect((dist.kurtosis ?? .nan) == 6)
        #expect((dist.kurtosisExcess ?? .nan) == 3)

        if let entropy = dist.entropy {
            let expected = Foundation.log(2 * Foundation.exp(1.0) * scale)
            #expect(abs(entropy - expected) <= 1e-12)
        } else {
            Issue.record("Expected Laplace entropy to be available")
        }
    }

    @Test("CDF/quantile round-trip and symmetry")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Laplace<Double>(location: 0.25, scale: 1.5)
        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round-trip mismatch for p=\(p)")
        }

        let offsets: [Double] = [0, 0.5, 1.0, 2.0]
        for offset in offsets {
            let xPlus = dist.location + offset
            let xMinus = dist.location - offset
            let pdfPlus = try dist.pdf(xPlus)
            let pdfMinus = try dist.pdf(xMinus)
            #expect(abs(pdfPlus - pdfMinus) <= 1e-12)

            let cdfPlus = try dist.cdf(xPlus)
            let cdfMinus = try dist.cdf(xMinus)
            #expect(abs((cdfPlus + cdfMinus) - 1) <= 1e-12)
        }
    }

    @Test("Hazard equals pdf/sf where defined")
    func hazardMatchesRatio() throws {
        let dist = try Distribution.Laplace<Double>(location: -0.5, scale: 0.9)
        let xs: [Double] = [-2.0, -0.5, 0.5, 2.0]
        for x in xs {
            let pdf = try dist.pdf(x)
            let sf = try dist.sf(x)
            let hazard = try dist.hazard(x)
            if sf > 0 {
                #expect(abs(hazard - (pdf / sf)) <= 1e-12)
            }
        }
    }
}
