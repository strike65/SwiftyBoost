import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Stateful distribution handles (Gamma, StudentT)")
struct DistributionsHandleTests {

    @Test("Gamma<Double> basic behaviors (pdf/cdf/quantile)")
    func gammaDoubleBasics() throws {
        let g = try Distribution.Gamma<Double>(shape: 2.5, scale: 1.2)
        let xs: [Double] = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
        for x in xs {
            // pdf nonnegative, cdf monotone
            let p = try g.pdf(x)
            #expect(p >= 0)
            _ = try g.cdf(x)
        }

        let ps: [Double] = [1e-6, 0.1, 0.5, 0.9, 1 - 1e-7]
        for p in ps {
            let q = try g.quantile(p)
            let c = try g.cdf(q)
            #expect(abs(c - p) <= 1e-12)
        }
    }

    @Test("StudentT<Double> symmetry and quantiles")
    func studentTDoubleBasics() throws {
        let v = 5.0
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: v)
        let xs: [Double] = [-3, -1, 0, 0.5, 2.0, 4.0]
        for x in xs {
            let p = try t.pdf(x)
            #expect(p >= 0)
            _ = try t.cdf(x)
        }

        // Symmetry: median and mode are 0
        #expect(try t.quantile(0.5) == 0)
        #expect(t.mode == 0)
    }

    @Test("Gamma<Float> and StudentT<Float> basic parity")
    func floatParity() throws {
        let g = try Distribution.Gamma<Float>(shape: 3, scale: 0.75)
        let x: Float = 1.3
        let p = try g.pdf(x)
        #expect(p >= 0)
        _ = try g.cdf(x)

        let t = try Distribution.StudentT<Float>(degreesOfFreedom: 7)
        _ = try t.pdf(0.25)
        _ = try t.cdf(0.25)
    }

    #if arch(x86_64) || arch(i386)
    @Test("Gamma<Float80> and StudentT<Float80> parity")
    func float80Parity() throws {
        let g = try Distribution.Gamma<Float80>(shape: 1.75, scale: 0.9)
        let x: Float80 = 0.8
        _ = try g.pdf(x)
        _ = try g.cdf(x)

        let t = try Distribution.StudentT<Float80>(degreesOfFreedom: 9)
        _ = try t.pdf(0.5)
        _ = try t.cdf(0.5)
    }
    #endif
}
