import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Dynamic distributions vs typed wrappers")
struct DynamicDistributionTests {

    @Test("Gamma<Double> dynamic matches typed for pdf/cdf/quantile")
    func gammaDynamicMatchesTypedDouble() throws {
        let shape: Double = 3.7, scale: Double = 0.8
        let dyn = try Distribution.Dynamic<Double>(distributionName: "gamma", parameters: ["shape": shape, "scale": scale])
        let typ = try Distribution.Gamma<Double>(shape: shape, scale: scale)

        let xs: [Double] = [1e-9, 1e-3, 0.1, 0.5, 1.0, 2.0]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)
            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)
        }

        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let xd = try dyn.quantile(p)
            let xt = try typ.quantile(p)
            #expect(abs(xd - xt) / max(xt, 1.0) <= 1e-12)
            let c = try dyn.cdf(xd)
            #expect(abs(c - p) <= 1e-10)
        }
    }

    @Test("StudentT<Float> dynamic matches typed for pdf/cdf/quantile")
    func studentTDynamicMatchesTypedFloat() throws {
        let df: Float = 12
        let dyn = try Distribution.Dynamic<Float>(distributionName: "student_t", parameters: ["df": df])
        let typ = try Distribution.StudentT<Float>(degreesOfFreedom: df)

        let xs: [Float] = [-2, -1, 0, 0.5, 1, 2]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-6)
            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-6)
        }

        let ps: [Float] = [1e-6, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let xd = try dyn.quantile(p)
            let xt = try typ.quantile(p)
            #expect(abs(xd - xt) / max(abs(xt), 1) <= 1e-5)
        }
    }

    @Test("FisherF<Double> dynamic matches typed")
    func fisherFDynamicMatchesTypedDouble() throws {
        let df1: Double = 5, df2: Double = 10
        let dyn = try Distribution.Dynamic<Double>(distributionName: "fisherf", parameters: ["df1": df1, "df2": df2])
        let typ = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: df2)

        let xs: [Double] = [0.01, 0.1, 0.5, 1.0, 2.0]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)
            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)
        }

        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let xd = try dyn.quantile(p)
            let xt = try typ.quantile(p)
            #expect(abs(xd - xt) / max(xt, 1.0) <= 1e-12)
        }
    }

    @Test("Arcsine<Float> dynamic matches typed")
    func arcsineDynamicMatchesTypedFloat() throws {
        let a: Float = 0, b: Float = 1
        let dyn = try Distribution.Dynamic<Float>(distributionName: "arcsine", parameters: ["minX": a, "maxX": b])
        let typ = try Distribution.Arcsine<Float>(minX: a, maxX: b)

        let xs: [Float] = [0.01, 0.1, 0.5, 0.9, 0.99]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            // arcsine has singularities at endpoints; avoid exact 0/1
            #expect(abs(pd - pt) <= 1e-5)
            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-5)
        }

        let ps: [Float] = [1e-6, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let xd = try dyn.quantile(p)
            let xt = try typ.quantile(p)
            #expect(abs(xd - xt) <= 1e-5)
        }
    }

    @Test("Geometric<Double> dynamic matches typed and closed forms")
    func geometricDynamicMatchesTypedDouble() throws {
        let p: Double = 0.37
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "geometric",
            parameters: ["theta": p] // alias for probability
        )
        let typ = try Distribution.Geometric<Double>(probabibilityOfSuccess: p)

        let ks: [Double] = [0, 1, 2, 5, 10]
        for k in ks {
            let pd = try dyn.pdf(k)
            let pt = try typ.pdf(k)
            #expect(abs(pd - pt) <= 1e-12)
            let expectedPdf = pow(1 - p, k) * p
            #expect(abs(pd - expectedPdf) <= 1e-12)

            let cd = try dyn.cdf(k)
            let ct = try typ.cdf(k)
            #expect(abs(cd - ct) <= 1e-12)
            let expectedCdf = 1 - pow(1 - p, k + 1)
            #expect(abs(cd - expectedCdf) <= 1e-12)

            let sd = try dyn.sf(k)
            let st = try typ.sf(k)
            #expect(abs(sd - st) <= 1e-12)
        }

        let probs: [Double] = [1e-9, 1e-6, 0.05, 0.5, 0.9, 1 - 1e-9]
        for prob in probs {
            let qd = try dyn.quantile(prob)
            let qt = try typ.quantile(prob)
            #expect(abs(qd - qt) <= 1e-10)
            let cd = try dyn.cdf(qd)
            #expect(cd >= prob - 1e-12)
        }

        let expectedMean = (1 - p) / p
        let expectedVariance = (1 - p) / (p * p)
        let expectedEntropy = -((1 - p) / p) * Foundation.log(1 - p) - Foundation.log(p)

        if let md = dyn.mean {
            #expect(abs(md - expectedMean) <= 1e-12)
        }
        if let vd = dyn.variance {
            #expect(abs(vd - expectedVariance) <= 1e-12)
        }
        if let entropy = dyn.entropy {
            #expect(abs(entropy - expectedEntropy) <= 1e-12)
        }
    }

    @Test("Holtsmark<Double> dynamic matches typed for core metrics")
    func holtsmarkDynamicMatchesTypedDouble() throws {
        let loc: Double = 0.25
        let scale: Double = 1.5
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "holtsmark",
            parameters: ["mu": loc, "sigma": scale]
        )
        let typ = try Distribution.Holtsmark<Double>(location: loc, scale: scale)

        let xs: [Double] = [-2.5, -0.75, 0.0, 0.75, 2.5]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)

            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)

            let sd = try dyn.sf(x)
            let st = try typ.sf(x)
            #expect(abs(sd - st) <= 1e-12)
        }

        let probs: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in probs {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-10)
        }

        // Location-scale properties
        if let mean = typ.mean {
            #expect(abs(mean - loc) <= 1e-12)
        }
        #expect(typ.median == loc)
        #expect((typ.mode ?? .nan) == loc)

        // Variance and higher moments diverge for Holtsmark (α = 3/2)
        #expect(typ.variance == nil)
        #expect(typ.skewness == nil)
        #expect(typ.kurtosis == nil)
        #expect(typ.kurtosisExcess == nil)

        // Differential entropy supplied in Swift using Constants.holtsmarkEntropy()
        if let entropy = typ.entropy {
            let expected = Foundation.log(scale) + Constants.holtsmarkEntropy()
            #expect(abs(entropy - expected) <= 1e-12)
        }

        // Hazard and CHF parity (only where survival probability is well-behaved)
        let interior: [Double] = [-0.5, 0.0, 0.5]
        for x in interior {
            let hd = try dyn.hazard(x)
            let ht = try typ.hazard(x)
            #expect(abs(hd - ht) <= 1e-12)

            let Hd = try dyn.chf(x)
            let Ht = try typ.chf(x)
            #expect(abs(Hd - Ht) <= 1e-12)
        }
    }

    @Test("Gamma<Double> moments and hazards")
    func gammaDynamicMomentsAndHazard() throws {
        let shape: Double = 5.0, scale: Double = 2.0
        let d = try Distribution.Dynamic<Double>(distributionName: "gamma", parameters: ["k": shape, "theta": scale])
        let t = try Distribution.Gamma<Double>(shape: shape, scale: scale)

        // Support and range
        #expect(d.supportLowerBound >= 0)
        #expect(d.supportUpperBound >= d.supportLowerBound)
        let r = d.range
        #expect(r.lower >= 0)
        #expect(r.upper >= r.lower)
        // Moments
        #expect(abs((d.mean ?? .nan) - (t.mean ?? .nan)) <= 1e-12)
        #expect(abs((d.variance ?? .nan) - (t.variance ?? .nan)) <= 1e-12)
        #expect(abs((d.kurtosis ?? .nan) - (t.kurtosis ?? .nan)) <= 1e-6)
        #expect(abs((d.kurtosisExcess ?? .nan) - (t.kurtosisExcess ?? .nan)) <= 1e-6)
        #expect(abs((d.skewness ?? .nan) - (t.skewness ?? .nan)) <= 1e-12)
        // Hazards
        let xs: [Double] = [0.5, 1, 3, 5]
        for x in xs {
            let hd = try d.hazard(x)
            let ht = try t.hazard(x)
            #expect(abs(hd - ht) <= 1e-12)
            let Hd = try d.chf(x)
            let Ht = try t.chf(x)
            #expect(abs(Hd - Ht) <= 1e-12)
        }
        // Lattice fields must be nil for continuous
        #expect(d.latticeStep == nil)
        #expect(d.latticeOrigin == nil)
    }

    @Test("StudentT alias handling and tails")
    func studentTAliasAndTails() throws {
        // Use alias "nu"
        let t = try Distribution.Dynamic<Double>(distributionName: "students_t", parameters: ["nu": 7.0])
        // Extreme tails
        let q = 1e-12
        let xu = try t.quantile(1 - q)
        let su = try t.sf(xu)
        #expect(abs(su - q) <= 1e-9)
        let xl = try t.quantile(q)
        let cl = try t.cdf(xl)
        #expect(abs(cl - q) <= 1e-9)
    }

    @Test("FisherF: entropy fallback and aliasing")
    func fisherFFallbackEntropyAndAliases() throws {
        // Provide df via aliases m/n
        let d = try Distribution.Dynamic<Double>(distributionName: "f", parameters: ["m": 10.0, "n": 20.0])
        let t = try Distribution.FisherF<Double>(degreesOfFreedom1: 10, degreesOfFreedom2: 20)
        // Entropy is Swift-side fallback in Dynamic; verify numerical agreement with typed wrapper
        if let ed = d.entropy, let et = t.entropy {
            #expect(abs(ed - et) <= 1e-10)
        }
        // Mode/Median present
        #expect((d.mode ?? .nan).isFinite)
        #expect(d.median.isFinite)
    }

    @Test("Arcsine: entropy fallback and aliases")
    func arcsineFallbackEntropyAndAliases() throws {
        let d = try Distribution.Dynamic<Float>(distributionName: "arcsine_distribution", parameters: ["lower": 0, "upper": 1])
        let t = try Distribution.Arcsine<Float>(minX: 0, maxX: 1)
        // Entropy fallback equals typed formula
        if let ed = d.entropy, let et = t.entropy {
            #expect(abs(ed - et) <= 1e-6)
        }
        // Hazard/CHF relationship for a few interior points
        let xs: [Float] = [0.1, 0.3, 0.7, 0.9]
        for x in xs {
            let hd = try d.hazard(x)
            let ht = try t.hazard(x)
            #expect(abs(hd - ht) <= 1e-5)
            let Hd = try d.chf(x)
            let Ht = try t.chf(x)
            #expect(abs(Hd - Ht) <= 1e-5)
        }
    }

    @Test("Factory errors for unknown name and missing parameters")
    func factoryErrors() throws {
        // Unknown distribution name
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Dynamic<Double>(distributionName: "not_a_dist", parameters: [:])
        }
        // Missing required parameter (gamma.shape)
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Dynamic<Double>(distributionName: "gamma", parameters: ["scale": 1.0])
        }
        // Invalid arcsine interval (min >= max)
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Dynamic<Double>(distributionName: "arcsine", parameters: ["minX": 1.0, "maxX": 0.5])
        }
        // Geometric probability outside [0, 1]
        #expect(throws: DistributionError<Double>.self) {
            _ = try Distribution.Dynamic<Double>(distributionName: "geometric", parameters: ["p": 1.5])
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 path smoke test")
    func float80Smoke() throws {
        let g = try Distribution.Dynamic<Float80>(distributionName: "gamma", parameters: ["k": 2.5, "theta": 1.1])
        let x: Float80 = 1.7
        let p = try g.pdf(x)
        let c = try g.cdf(x)
        #expect(p.isFinite && c >= 0 && c <= 1)
    }
    #endif
}
