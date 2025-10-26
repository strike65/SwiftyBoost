import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Logistic distribution (Double)")
struct LogisticDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.Logistic<Double>(location: 0, scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -1)
        ) {
            _ = try Distribution.Logistic<Double>(location: 0, scale: -1)
        }
    }

    @Test("CDF/quantile round-trip and symmetry")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Logistic<Double>(location: 0.25, scale: 1.8)
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
            #expect(abs((cdfPlus + cdfMinus) - 1.0) <= 1e-12)
        }
    }

    @Test("Closed-form moments and entropy")
    func analyticMoments() throws {
        let location: Double = -1.0
        let scale: Double = 0.75
        let dist = try Distribution.Logistic<Double>(location: location, scale: scale)

        #expect(dist.mean == location)
        if let variance = dist.variance {
            let expectedVar = (Double.pi * Double.pi * scale * scale) / 3.0
            #expect(abs(variance - expectedVar) <= 1e-12)
        } else {
            Issue.record("Expected logistic variance to be available")
        }

        #expect(dist.mode == location)
        #expect(dist.median == location)
        #expect((dist.skewness ?? .nan) == 0)
        if let kurtosis = dist.kurtosis {
            let expectedKurtosis = 21.0 / 5.0 // 4.2
            #expect(abs(kurtosis - expectedKurtosis) <= 1e-12)
        } else {
            Issue.record("Expected logistic kurtosis to be available")
        }
        if let kurtosisExcess = dist.kurtosisExcess {
            let expectedExcess = 21.0 / 5.0 - 3.0
            #expect(abs(kurtosisExcess - expectedExcess) <= 1e-12)
        } else {
            Issue.record("Expected logistic kurtosisExcess to be available")
        }

        if let entropy = dist.entropy {
            let expectedEntropy = Foundation.log(scale) + 2.0
            #expect(abs(entropy - expectedEntropy) <= 1e-12)
        } else {
            Issue.record("Expected logistic entropy to be available")
        }
    }

    @Test("Hazard equals pdf/sf where defined")
    func hazardMatchesRatio() throws {
        let dist = try Distribution.Logistic<Double>(location: -0.2, scale: 1.4)
        let xs: [Double] = [-3.0, -1.0, 0.0, 1.0, 3.0]
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
