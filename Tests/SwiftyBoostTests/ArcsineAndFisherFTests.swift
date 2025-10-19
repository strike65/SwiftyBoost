import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Arcsine distribution tests")
struct ArcsineDistributionTests {

    @Test("Arcsine[0,1] basic properties")
    func arcsineUnitInterval() throws {
        let a: Double = 0.0, b: Double = 1.0
        let d = try Distribution.Arcsine<Double>(minX: a, maxX: b)

        // Support and range
        #expect(d.supportLowerBound == a)
        #expect(d.supportUpperBound == b)
        #expect(d.range.lower == a)
        #expect(d.range.upper == b)

        // Mean and variance
        #expect(abs((d.mean ?? .nan) - 0.5) <= 1e-12)
        #expect(abs((d.variance ?? .nan) - (1.0/8.0)) <= 1e-12)

        // Midpoint density = 2/pi
        let m = 0.5
        let pdfMid = try d.pdf(m)
        #expect(abs(pdfMid - 2.0/Double.pi) <= 1e-12)

        // Symmetry at midpoint
        #expect(abs((try d.cdf(m)) - 0.5) <= 1e-12)
        #expect(abs((try d.quantile(0.5)) - m) <= 1e-12)
    }

    @Test("Arcsine[a,b] midpoint PDF and round-trip")
    func arcsineGeneralInterval() throws {
        let a: Double = 2.0, b: Double = 6.0
        let d = try Distribution.Arcsine<Double>(minX: a, maxX: b)
        let mid = (a + b) / 2
        // PDF(mid) = 2 / (pi * (b - a))
        let expected = 2.0 / (Double.pi * (b - a))
        #expect(abs((try d.pdf(mid)) - expected) <= 1e-12)

        // Quantile/CDF round-trip
        let ps: [Double] = [1e-6, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try d.quantile(p)
            let c = try d.cdf(x)
            #expect(abs(c - p) <= 1e-10)
        }
    }
}

@Suite("Arcsine distribution Float/Float80 tests")
struct ArcsineFloatAndFloat80Tests {

    @Test("Arcsine<Float>[0,1] endpoints and round-trip")
    func arcsineFloatUnit() throws {
        let d = try Distribution.Arcsine<Float>(minX: 0, maxX: 1)
        #expect(d.supportLowerBound == 0)
        #expect(d.supportUpperBound == 1)
        #expect(try d.cdf(0) == 0)
        #expect(try d.cdf(1) == 1)
        #expect(try d.quantile(0) == 0)
        #expect(try d.quantile(1) == 1)
        let ps: [Float] = [1e-4, 0.1, 0.5, 0.9, 1 - 1e-4]
        for p in ps {
            let x = try d.quantile(p)
            let c = try d.cdf(x)
            #expect(abs(c - p) <= 2e-4)
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Arcsine<Float80>[2,6] midpoint and round-trip")
    func arcsineFloat80General() throws {
        let a: Float80 = 2, b: Float80 = 6
        let d = try Distribution.Arcsine<Float80>(minX: a, maxX: b)
        let mid = (a + b)/2
        let expected = 2.0 / (Float80.pi * (b - a))
        #expect(abs((try d.pdf(mid)) - expected) <= 1e-14)
        let ps: [Float80] = [1e-10, 0.1, 0.5, 0.9, 1 - 1e-10]
        for p in ps {
            let x = try d.quantile(p)
            let c = try d.cdf(x)
            #expect(abs(c - p) <= 1e-12)
        }
    }
    #endif
}

@Suite("Fisher–Snedecor F distribution tests")
struct FisherFDistributionTests {

    private func fisherMean(df1: Double, df2: Double) -> Double? {
        guard df2 > 2 else { return nil }
        return df2 / (df2 - 2)
    }

    private func fisherVariance(df1: Double, df2: Double) -> Double? {
        guard df2 > 4 else { return nil }
        let num = 2 * df2 * df2 * (df1 + df2 - 2)
        let den = df1 * (df2 - 2) * (df2 - 2) * (df2 - 4)
        return num / den
    }

    private func fisherMode(df1: Double, df2: Double) -> Double? {
        guard df1 > 2 else { return nil }
        return ((df1 - 2.0) / df1) * (df2 / (df2 + 2.0))
    }

    @Test("FisherF(df1=5, df2=10) moments and round-trip")
    func fisherFMomentsAndRoundTrip() throws {
        let df1 = 5.0, df2 = 10.0
        let f = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: df2)

        // Support and range
        #expect(f.supportLowerBound >= 0)
        #expect(f.supportUpperBound > f.supportLowerBound)

        // Mean and variance (exist for df2>2 and df2>4)
        if let m = fisherMean(df1: df1, df2: df2) {
            #expect(abs((f.mean ?? .nan) - m) <= 1e-12)
        }
        if let v = fisherVariance(df1: df1, df2: df2) {
            #expect(abs((f.variance ?? .nan) - v) <= 1e-10)
        }
        if let mo = fisherMode(df1: df1, df2: df2), let fm = f.mode, fm.isFinite {
            #expect(abs(fm - mo) <= 1e-10)
        }

        // CDF/Quantile round-trip
        let ps: [Double] = [1e-6, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let x = try f.quantile(p)
            let c = try f.cdf(x)
            #expect(abs(c - p) <= 1e-10)
        }

        // Hazard identity at safe interior points
        let xs: [Double] = [0.5, 1.0, 2.0, 3.0]
        for x in xs {
            let p = try f.pdf(x)
            let s = try f.sf(x)
            let h = try f.hazard(x)
            let expected = p / s
            if s > 0 && p.isFinite && s.isFinite && expected.isFinite && h.isFinite {
                #expect(abs(h - expected) <= 1e-8)
            }
        }
    }

    @Test("FisherF<Double> existence boundaries for moments")
    func fisherFExistenceBoundaries() throws {
        let df1 = 5.0
        // Mean undefined at df2 ≤ 2, defined just above
        do {
            let fEdge = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 2.0)
            #expect(fEdge.mean == nil)
        }
        do {
            let fNear = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 2.0000001)
            #expect((fNear.mean ?? .nan).isFinite)
        }
        // Variance undefined at df2 ≤ 4, defined just above
        do {
            let fEdge = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 4.0)
            #expect(fEdge.variance == nil)
        }
        do {
            let fNear = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 4.0000001)
            #expect((fNear.variance ?? .nan).isFinite)
        }
        // Skewness undefined at df2 ≤ 6, defined just above
        do {
            let fEdge = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 6.0)
            #expect(fEdge.skewness == nil)
        }
        do {
            let fNear = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 6.0000001)
            #expect((fNear.skewness ?? .nan).isFinite)
        }
        // Kurtosis undefined at df2 ≤ 8, defined just above
        do {
            let fEdge = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 8.0)
            #expect(fEdge.kurtosis == nil)
            #expect(fEdge.kurtosisExcess == nil)
        }
        do {
            let fNear = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: 8.0000001)
            #expect((fNear.kurtosis ?? .nan).isFinite)
            #expect((fNear.kurtosisExcess ?? .nan).isFinite)
        }
    }

    @Test("FisherF<Double> extreme tail quantiles (1e-12)")
    func fisherFDoubleExtremeTails() throws {
        let df1 = 5.0, df2 = 10.0
        let f = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: df2)

        // Upper tail q = 1e-12
        let q: Double = 1e-12
        let pUpper = 1 - q
        let xu1 = try f.quantile(pUpper)
        let su1 = try f.sf(xu1)
        #expect(xu1.isFinite && xu1 > 0)
        #expect(abs(su1 - q) <= 1e-12)

        // Direct upper-tail quantile should match closely
        let xu2 = try f.quantileComplement(q)
        let su2 = try f.sf(xu2)
        #expect(abs(su2 - q) <= 1e-12)
        // Relative agreement of the two x’s within a tight bound
        let rel = abs(xu2 - xu1) / max(xu1, 1.0)
        #expect(rel <= 1e-5)

        // Lower tail p = 1e-12
        let pLow: Double = 1e-12
        let xl = try f.quantile(pLow)
        let cl = try f.cdf(xl)
        #expect(xl >= 0)
        #expect(abs(cl - pLow) <= 1e-12)
    }
}

@Suite("Fisher F Float/Float80 edge cases")
struct FisherFFloatAndFloat80Tests {

    @Test("FisherF<Float> endpoints and round-trip")
    func fisherFFloat() throws {
        let f = try Distribution.FisherF<Float>(degreesOfFreedom1: 5, degreesOfFreedom2: 10)
        #expect(f.supportLowerBound >= 0)
        #expect(try f.cdf(0) == 0)
        #expect(try f.quantile(0) == 0)
        // upper tail: quantile(1) should be +inf
        let q1 = try f.quantile(1)
        #expect(q1.isInfinite && q1 > 0)
        // round-trip
        let ps: [Float] = [1e-5, 0.1, 0.5, 0.9, 1 - 1e-5]
        for p in ps {
            let x = try f.quantile(p)
            let c = try f.cdf(x)
            #expect(abs(c - p) <= 1e-5)
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("FisherF<Float80> endpoints and round-trip")
    func fisherFFloat80() throws {
        let f = try Distribution.FisherF<Float80>(degreesOfFreedom1: 6, degreesOfFreedom2: 12)
        #expect(f.supportLowerBound >= 0)
        #expect(try f.cdf(0) == 0)
        #expect(try f.quantile(0) == 0)
        let q1 = try f.quantile(1)
        #expect(q1.isInfinite && q1 > 0)
        let ps: [Float80] = [1e-12, 0.1, 0.5, 0.9, 1 - 1e-12]
        for p in ps {
            let x = try f.quantile(p)
            let c = try f.cdf(x)
            #expect(abs(c - p) <= 1e-12)
        }
    }
    #endif
}
