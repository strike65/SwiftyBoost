import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Holtsmark distribution (Double)")
struct HoltsmarkDistributionTests {

    @Test("PDF/CDF symmetry and shift invariance")
    func symmetryAndShift() throws {
        let holts = try Distribution.Holtsmark<Double>(loc: 0.4, scale: 1.2)

        let offsets: [Double] = [0.0, 0.5, 1.0, 1.5]
        for offset in offsets {
            let xPlus = holts.location + offset
            let xMinus = holts.location - offset

            let pdfPlus = try holts.pdf(xPlus)
            let pdfMinus = try holts.pdf(xMinus)
            #expect(abs(pdfPlus - pdfMinus) <= 1e-12)

            let cdfPlus = try holts.cdf(xPlus)
            let cdfMinus = try holts.cdf(xMinus)
            #expect(abs((cdfPlus + cdfMinus) - 1.0) <= 1e-12)
        }
    }

    @Test("CDF/quantile round-trip")
    func quantileRoundTrip() throws {
        let holts = try Distribution.Holtsmark<Double>(loc: 0, scale: 1)
        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try holts.quantile(p)
            let c = try holts.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round-trip mismatch for p=\(p)")
        }
    }

    @Test("Moments and entropy availability")
    func momentsAndEntropy() throws {
        let scale: Double = 2.5
        let holts = try Distribution.Holtsmark<Double>(loc: -0.75, scale: scale)

        // Mean/median/mode coincide with the location parameter.
        #expect(abs((holts.mean ?? .nan) - holts.location) <= 1e-12)
        #expect(holts.median == holts.location)
        #expect((holts.mode ?? .nan) == holts.location)

        // Variance and higher moments diverge for Holtsmark.
        #expect(holts.variance == nil)
        #expect(holts.skewness == nil)
        #expect(holts.kurtosis == nil)
        #expect(holts.kurtosisExcess == nil)

        // Differential entropy is provided analytically: H = log(scale) + H_standard.
        if let entropy = holts.entropy {
            let expected = Foundation.log(scale) + Constants.holtsmarkEntropy()
            #expect(abs(entropy - expected) <= 1e-12)
        } else {
            Issue.record("Expected Holtsmark entropy to be available")
        }
    }
}
