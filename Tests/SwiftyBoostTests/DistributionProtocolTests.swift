import Testing
@testable import SwiftyBoost
import Foundation

@Suite("DistributionProtocol – Gamma")
struct GammaProtocolTests {
    @Test("support and range endpoints")
    func supportAndRange() throws {
        let g = try Distribution.Gamma<Double>(shape: 2.5, scale: 1.2)
        #expect(g.supportLowerBound >= 0)
        #expect(g.supportUpperBound > g.supportLowerBound)
        let r = g.range
        // Boost may report an open lower support bound ~ DBL_MIN, while range may report 0.
        #expect(r.lower <= g.supportLowerBound + 1e-300)
        #expect(r.upper <= g.supportUpperBound || r.upper.isFinite)
    }

    @Test("cdf/quantile round-trip with tolerance")
    func cdfQuantileRoundTrip() throws {
        let g = try Distribution.Gamma<Double>(shape: 3.2, scale: 0.8)
        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 0.999, 1-1e-9]
        for p in ps {
            let x = try g.quantile(p)
            let c = try g.cdf(x)
            #expect(abs(c - p) <= 1e-10, "p=\(p), got cdf(quantile)=\(c)")
        }
    }

    @Test("hazard and CHF relationships")
    func hazardAndChf() throws {
        let g = try Distribution.Gamma<Double>(shape: 4.0, scale: 1.5)
        let xs: [Double] = [0.5, 1.0, 1.5, 2.0]
        for x in xs {
            let h = try g.hazard(x)
            let p = try g.pdf(x)
            let s = try g.sf(x)
            let expected = p / s
            if s > 0 && p.isFinite && s.isFinite && expected.isFinite && h.isFinite {
                #expect(abs(h - expected) <= 1e-10)
                let H = try g.chf(x)
                if s > 1e-300 {
                    let expH = -Foundation.log(s)
                    if expH.isFinite && H.isFinite {
                        #expect(abs(H - expH) <= 1e-10)
                    }
                }
            }
        }
    }
}

@Suite("DistributionProtocol – StudentT")
struct StudentTProtocolTests {
    @Test("support and range infinite endpoints")
    func supportAndRange() throws {
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: 5)
        #expect(t.supportLowerBound.isInfinite && t.supportLowerBound < 0)
        #expect(t.supportUpperBound.isInfinite && t.supportUpperBound > 0)
        let r = t.range
        #expect(r.lower.isInfinite && r.lower < 0)
        #expect(r.upper.isInfinite && r.upper > 0)
    }

    @Test("moments for v=5")
    func momentsV5() throws {
        let v: Double = 5
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: v)
        #expect(t.mode == 0)
        #expect(t.median == 0)
        // mean exists and is 0
        #expect(t.mean == 0)
        // variance = v/(v-2)
        #expect(abs((t.variance ?? .nan) - (v/(v-2))) <= 1e-12)
    }

    @Test("cdf/quantile round-trip with tolerance")
    func cdfQuantileRoundTrip() throws {
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: 7)
        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 0.999, 1-1e-9]
        for p in ps {
            let x = try t.quantile(p)
            let c = try t.cdf(x)
            #expect(abs(c - p) <= 1e-10, "p=\(p), got cdf(quantile)=\(c)")
        }
    }

    @Test("hazard and CHF relationships")
    func hazardAndChf() throws {
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: 9)
        let xs: [Double] = [-1.5, -0.5, 0.0, 0.5, 1.5]
        for x in xs {
            let h = try t.hazard(x)
            let p = try t.pdf(x)
            let s = try t.sf(x)
            let expected = p / s
            if s > 0 && p.isFinite && s.isFinite && expected.isFinite && h.isFinite {
                #expect(abs(h - expected) <= 1e-10)
                let H = try t.chf(x)
                if s > 1e-300 {
                    let expH = -Foundation.log(s)
                    if expH.isFinite && H.isFinite {
                        #expect(abs(H - expH) <= 1e-10)
                    }
                }
            }
        }
    }

    @Test("undefined moments become nil")
    func undefinedMoments() throws {
        // v <= 1: mean undefined
        let t1 = try Distribution.StudentT<Double>(degreesOfFreedom: 0.75)
        #expect(t1.mean == nil)
        // 1 < v <= 2: variance diverges (∞) – we report nil
        let t2 = try Distribution.StudentT<Double>(degreesOfFreedom: 1.5)
        #expect(t2.variance == nil)
    }

    @Test("StudentT<Double> extreme tail quantiles (1e-15)")
    func studentTExtremeTails() throws {
        let v: Double = 10
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: v)
        // Upper tail q = 1e-15
        let q: Double = 1e-15
        let pUpper = 1 - q
        let xu = try t.quantile(pUpper)
        let su = try t.sf(xu)
        #expect(xu.isFinite)
        #expect(abs(su - q) <= 1e-12)
        // Compare with direct upper-tail quantile
        let xu2 = try t.quantileComplement(q)
        let su2 = try t.sf(xu2)
        #expect(abs(su2 - q) <= 1e-12)
        let rel = abs(xu2 - xu) / max(abs(xu), 1.0)
        #expect(rel <= 1e-4)

        // Lower tail p = 1e-15
        let pLow: Double = 1e-15
        let xl = try t.quantile(pLow)
        let cl = try t.cdf(xl)
        #expect(abs(cl - pLow) <= 1e-12)
    }
}
