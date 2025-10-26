import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Map-Airy distribution (Double)")
struct MapAiryDistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.MapAiry<Double>(location: 0, scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -0.1)
        ) {
            _ = try Distribution.MapAiry<Double>(location: 0, scale: -0.1)
        }
    }

    @Test("CDF/quantile round-trip")
    func quantileRoundTrip() throws {
        let dist = try Distribution.MapAiry<Double>(location: -0.25, scale: 1.2)
        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-8, "Round-trip mismatch for p=\(p)")
        }
    }

    @Test("Moments and entropy behaviour")
    func momentBehaviour() throws {
        let scale1: Double = 0.75
        let scale2: Double = 2.25
        let d1 = try Distribution.MapAiry<Double>(location: 0.6, scale: scale1)
        let d2 = try Distribution.MapAiry<Double>(location: 0.6, scale: scale2)

        // Mean exists and equals the location.
        #expect(abs((d1.mean ?? .nan) - d1.location) <= 1e-12)
        #expect(abs((d2.mean ?? .nan) - d2.location) <= 1e-12)

        // Variance and higher moments diverge.
        #expect(d1.variance == nil)
        #expect(d1.skewness == nil)
        #expect(d1.kurtosis == nil)
        #expect(d1.kurtosisExcess == nil)

        // Entropy scales with log(scale).
        if let e1 = d1.entropy, let e2 = d2.entropy {
            let expected = Foundation.log(scale2 / scale1)
            #expect(abs((e2 - e1) - expected) <= 1e-10)
        } else {
            Issue.record("Expected Map-Airy entropy to be available")
        }

        // Mode/median finite.
        #expect((d1.mode ?? .nan).isFinite)
        #expect(d1.median.isFinite)
    }
}
