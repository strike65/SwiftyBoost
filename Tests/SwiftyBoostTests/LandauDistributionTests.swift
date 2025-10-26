import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Landau distribution (Double)")
struct LandauDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.Landau<Double>(location: 0, scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -1)
        ) {
            _ = try Distribution.Landau<Double>(location: 0, scale: -1)
        }
    }

    @Test("CDF/quantile round-trip and tail complement")
    func quantileRoundTrip() throws {
        let dist = try Distribution.Landau<Double>(location: 0.5, scale: 1.7)
        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10, "Round-trip mismatch for p=\(p)")
        }

        let qs: [Double] = [1e-5, 1e-3, 0.1]
        for q in qs {
            let x = try dist.quantileComplement(q)
            let s = try dist.sf(x)
            #expect(abs(s - q) <= 1e-10, "Round-trip mismatch for q=\(q)")
        }
    }

    @Test("Location-shift invariance for mode and median")
    func locationShiftInvariance() throws {
        let base = try Distribution.Landau<Double>(location: 0, scale: 1.1)
        let shifted = try Distribution.Landau<Double>(location: 1.25, scale: 1.1)

        if let baseMode = base.mode, let shiftedMode = shifted.mode {
            #expect(abs((shiftedMode - baseMode) - 1.25) <= 1e-12)
        } else {
            Issue.record("Expected Landau mode to be available")
        }

        let baseMedian = base.median
        let shiftedMedian = shifted.median
        #expect(abs((shiftedMedian - baseMedian) - 1.25) <= 1e-12)
    }

    @Test("Heavy-tailed moments and entropy scaling")
    func momentsAndEntropy() throws {
        let scale1: Double = 0.75
        let scale2: Double = 2.0
        let d1 = try Distribution.Landau<Double>(location: -0.3, scale: scale1)
        let d2 = try Distribution.Landau<Double>(location: -0.3, scale: scale2)

        #expect(d1.mean == nil)
        #expect(d1.variance == nil)
        #expect(d1.skewness == nil)
        #expect(d1.kurtosis == nil)
        #expect(d1.kurtosisExcess == nil)

        if let e1 = d1.entropy, let e2 = d2.entropy {
            let diff = e2 - e1
            let expectedDiff = Foundation.log(scale2 / scale1)
            #expect(abs(diff - expectedDiff) <= 1e-10)
        } else {
            Issue.record("Expected Landau entropy to be available")
        }
    }
}
