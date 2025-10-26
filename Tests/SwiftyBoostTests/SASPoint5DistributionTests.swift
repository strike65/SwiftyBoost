import Testing
@testable import SwiftyBoost

@Suite("SAS Point5 distribution (Double)")
struct SASPoint5DistributionTests {

    @Test("Scale must be positive")
    func invalidScaleThrows() {
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: 0)
        ) {
            _ = try Distribution.SASPoint5<Double>(scale: 0)
        }
        #expect(
            throws: DistributionError<Double>.parameterNotPositive(name: "scale", value: -2)
        ) {
            _ = try Distribution.SASPoint5<Double>(scale: -2)
        }
    }

    @Test("Undefined moments propagate as nil")
    func undefinedMoments() throws {
        let dist = try Distribution.SASPoint5<Double>(location: 1.2, scale: 0.75)
        #expect(dist.mean == nil)
        #expect(dist.variance == nil)
        #expect(dist.skewness == nil)
        #expect(dist.kurtosis == nil)
        #expect(dist.kurtosisExcess == nil)
        #expect(dist.mode == dist.location)
        #expect(abs(dist.median - dist.location) <= 1e-12)
        #expect(dist.entropy != nil)
    }

    @Test("Symmetry around the location parameter")
    func symmetry() throws {
        let dist = try Distribution.SASPoint5<Double>(location: -0.4, scale: 1.3)
        let offsets: [Double] = [0, 0.25, 0.5, 1.0, 2.5]
        for offset in offsets {
            let plus = dist.location + offset
            let minus = dist.location - offset
            let pdfPlus = try dist.pdf(plus)
            let pdfMinus = try dist.pdf(minus)
            #expect(abs(pdfPlus - pdfMinus) <= 1e-12)

            let cdfPlus = try dist.cdf(plus)
            let cdfMinus = try dist.cdf(minus)
            #expect(abs((cdfPlus + cdfMinus) - 1) <= 1e-12)
        }
    }

    @Test("CDF/quantile round trip")
    func cdfQuantileRoundTrip() throws {
        let dist = try Distribution.SASPoint5<Double>(location: 0.1, scale: 0.9)
        let ps: [Double] = [1e-9, 1e-5, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try dist.quantile(p)
            let c = try dist.cdf(x)
            #expect(abs(c - p) <= 1e-10)
        }
    }
}
