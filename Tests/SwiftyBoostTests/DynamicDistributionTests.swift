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
        let typ = try Distribution.Holtsmark<Double>(loc: loc, scale: scale)

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

    @Test("Kolmogorov–Smirnov<Double> dynamic matches typed")
    func kolmogorovSmirnovDynamicMatchesTypedDouble() throws {
        let n: Double = 32
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "kolmogorov_smirnov",
            parameters: ["n": n]
        )
        let typ = try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: n)

        let xs: [Double] = [0.05, 0.1, 0.2, 0.4, 0.8]
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

        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-10)
        }

        if let meanDyn = dyn.mean, let meanTyp = typ.mean {
            #expect(abs(meanDyn - meanTyp) <= 1e-12)
        }
        if let varianceDyn = dyn.variance, let varianceTyp = typ.variance {
            #expect(abs(varianceDyn - varianceTyp) <= 1e-12)
        }
        if let skewDyn = dyn.skewness, let skewTyp = typ.skewness {
            #expect(abs(skewDyn - skewTyp) <= 1e-12)
        }
        if let kurtDyn = dyn.kurtosis, let kurtTyp = typ.kurtosis {
            #expect(abs(kurtDyn - kurtTyp) <= 1e-12)
        }
    }

    @Test("Landau<Double> dynamic matches typed")
    func landauDynamicMatchesTypedDouble() throws {
        let location: Double = -0.75
        let scale: Double = 1.6
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "landau",
            parameters: ["mu": location, "c": scale]
        )
        let typ = try Distribution.Landau<Double>(location: location, scale: scale)

        let xs: [Double] = [-3.0, -1.0, 0.0, 1.0, 3.0]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)

            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)
        }

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-10)
        }

        #expect(dyn.mean == nil && typ.mean == nil)
        #expect(dyn.variance == nil && typ.variance == nil)
        #expect(dyn.skewness == nil && typ.skewness == nil)
        #expect(dyn.kurtosis == nil && typ.kurtosis == nil)
        if let entropyDyn = dyn.entropy, let entropyTyp = typ.entropy {
            #expect(abs(entropyDyn - entropyTyp) <= 1e-12)
        }
    }

    @Test("Laplace<Double> dynamic matches typed")
    func laplaceDynamicMatchesTypedDouble() throws {
        let location: Double = 0.3
        let scale: Double = 0.9
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "laplace",
            parameters: ["location": location, "scale": scale]
        )
        let typ = try Distribution.Laplace<Double>(location: location, scale: scale)

        let xs: [Double] = [-1.5, -0.3, 0.0, 0.3, 1.5]
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
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-12)
        }

        #expect(dyn.mean == typ.mean)
        #expect(dyn.variance == typ.variance)
        #expect(dyn.skewness == typ.skewness)
        #expect(dyn.kurtosis == typ.kurtosis)
        #expect(dyn.kurtosisExcess == typ.kurtosisExcess)
        if let entropyDyn = dyn.entropy, let entropyTyp = typ.entropy {
            #expect(abs(entropyDyn - entropyTyp) <= 1e-12)
        }
    }

    @Test("Logistic<Double> dynamic matches typed")
    func logisticDynamicMatchesTypedDouble() throws {
        let location: Double = -0.4
        let scale: Double = 1.25
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "logistic",
            parameters: ["location": location, "scale": scale]
        )
        let typ = try Distribution.Logistic<Double>(location: location, scale: scale)

        let xs: [Double] = [-3.0, -1.0, 0.0, 1.0, 3.0]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)

            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)

            let hd = try dyn.hazard(x)
            let ht = try typ.hazard(x)
            #expect(abs(hd - ht) <= 1e-12)
        }

        let ps: [Double] = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-12)
        }

        #expect(dyn.mean == typ.mean)
        #expect(dyn.variance == typ.variance)
        #expect(dyn.skewness == typ.skewness)
        #expect(dyn.kurtosis == typ.kurtosis)
        #expect(dyn.kurtosisExcess == typ.kurtosisExcess)
        if let eDyn = dyn.entropy, let eTyp = typ.entropy {
            #expect(abs(eDyn - eTyp) <= 1e-12)
        }
    }

    @Test("LogNormal<Double> dynamic matches typed")
    func logNormalDynamicMatchesTypedDouble() throws {
        let location: Double = 0.5
        let scale: Double = 0.75
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "lognormal",
            parameters: ["location": location, "scale": scale]
        )
        let typ = try Distribution.LogNormal<Double>(location: location, scale: scale)

        let xs: [Double] = [0.0001, 0.01, 0.1, 0.5, 2.0]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)

            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)
        }

        let ps: [Double] = [1e-12, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-10)
        }

        #expect(dyn.mean == typ.mean)
        #expect(dyn.variance == typ.variance)
        #expect(dyn.skewness == typ.skewness)
        #expect(dyn.kurtosis == typ.kurtosis)
        #expect(dyn.kurtosisExcess == typ.kurtosisExcess)
        if let eDyn = dyn.entropy, let eTyp = typ.entropy {
            #expect(abs(eDyn - eTyp) <= 1e-12)
        }
    }

    @Test("Map-Airy<Double> dynamic matches typed")
    func mapAiryDynamicMatchesTypedDouble() throws {
        let location: Double = 0.35
        let scale: Double = 1.4
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "mapairy",
            parameters: ["location": location, "scale": scale]
        )
        let typ = try Distribution.MapAiry<Double>(location: location, scale: scale)

        let xs: [Double] = [-3.0, -1.0, 0.0, 1.0, 3.0]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-10)

            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-10)
        }

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-8)
        }

        #expect(abs((dyn.mean ?? .nan) - (typ.mean ?? .nan)) <= 1e-12)
        #expect(dyn.variance == nil && typ.variance == nil)
        #expect(dyn.skewness == nil && typ.skewness == nil)
        #expect(dyn.kurtosis == nil && typ.kurtosis == nil)
        if let eDyn = dyn.entropy, let eTyp = typ.entropy {
            #expect(abs(eDyn - eTyp) <= 1e-10)
        }
    }

    @Test("NegativeBinomial<Double> dynamic matches typed")
    func negativeBinomialDynamicMatchesTypedDouble() throws {
        let successes: Double = 3.2
        let p: Double = 0.37
        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "negative_binomial",
            parameters: ["successes": successes, "p": p]
        )
        let typ = try Distribution.NegativeBinomial<Double>(successes: successes, probabilityOfSuccess: p)

        let xs: [Double] = [0, 1, 2, 5, 10]
        for x in xs {
            let pd = try dyn.pdf(x)
            let pt = try typ.pdf(x)
            #expect(abs(pd - pt) <= 1e-12)

            let cd = try dyn.cdf(x)
            let ct = try typ.cdf(x)
            #expect(abs(cd - ct) <= 1e-12)

            let hd = try dyn.hazard(Double(x))
            let ht = try typ.hazard(Double(x))
            if hd.isFinite && ht.isFinite {
                #expect(abs(hd - ht) <= 1e-12)
            }
        }

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for prob in ps {
            let qd = try dyn.quantile(prob)
            let qt = try typ.quantile(prob)
            #expect(abs(qd - qt) <= 1e-10)
        }

        #expect(abs((dyn.mean ?? .nan) - (typ.mean ?? .nan)) <= 1e-12)
        #expect(abs((dyn.variance ?? .nan) - (typ.variance ?? .nan)) <= 1e-12)
        #expect(abs((dyn.skewness ?? .nan) - (typ.skewness ?? .nan)) <= 1e-12)
        #expect(abs((dyn.kurtosis ?? .nan) - (typ.kurtosis ?? .nan)) <= 1e-12)
        if let eDyn = dyn.entropy, let eTyp = typ.entropy {
            #expect(abs(eDyn - eTyp) <= 1e-12)
        }
    }

    @Test("Hyperexponential<Double> dynamic matches typed")
    func hyperexponentialDynamicMatchesTypedDouble() throws {
        let probabilities: [Double] = [0.15, 0.55, 0.3]
        let rates: [Double] = [0.65, 1.45, 3.25]

        var params: [String: Double] = [:]
        for (idx, rate) in rates.enumerated() {
            params["rate\(idx)"] = rate
        }
        for (idx, probability) in probabilities.enumerated() {
            params["prob\(idx)"] = probability
        }

        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "hyperexponential",
            parameters: params
        )
        let typ = try Distribution.Hyperexponential<Double>(
            probabilities: probabilities,
            rates: rates
        )

        let xs: [Double] = [0.0, 0.05, 0.2, 0.75, 1.2]
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

        let ps: [Double] = [1e-9, 1e-6, 0.01, 0.5, 0.9, 1 - 1e-9]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-11)
        }

        let x: Double = 0.8
        let hd = try dyn.hazard(x)
        let ht = try typ.hazard(x)
        #expect(abs(hd - ht) <= 1e-12)

        let cd = try dyn.chf(x)
        let ct = try typ.chf(x)
        #expect(abs(cd - ct) <= 1e-12)
    }

    @Test("InverseChiSquared<Double> dynamic matches typed")
    func inverseChiSquaredDynamicMatchesTypedDouble() throws {
        let df: Double = 8.5
        let scale: Double = 0.45

        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "inverse_chi_squared",
            parameters: [
                "df": df,
                "scale": scale
            ]
        )
        let typ = try Distribution.InverseChiSquared<Double>(
            degreesOfFreedom: df,
            scale: scale
        )

        let xs: [Double] = [0.05, 0.1, 0.25, 0.75, 1.5]
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

        let ps: [Double] = [1e-9, 1e-6, 0.05, 0.5, 0.95, 1 - 1e-9]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-11)
        }

        let x: Double = 0.4
        let hd = try dyn.hazard(x)
        let ht = try typ.hazard(x)
        #expect(abs(hd - ht) <= 1e-12)

        let cd = try dyn.chf(x)
        let ct = try typ.chf(x)
        #expect(abs(cd - ct) <= 1e-12)
    }

    @Test("NonCentralChiSquared<Double> dynamic matches typed")
    func nonCentralChiSquaredDynamicMatchesTypedDouble() throws {
        let df: Double = 6.0
        let lambda: Double = 3.5

        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "non_central_chi_squared",
            parameters: [
                "df": df,
                "lambda": lambda
            ]
        )
        let typ = try Distribution.NonCentralChiSquared<Double>(
            degreesOfFreedom: df,
            nonCentrality: lambda
        )

        let xs: [Double] = [0.1, 0.5, 1.0, 2.5, 5.0, 8.0]
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

        let ps: [Double] = [1e-6, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-6]
        for p in ps {
            let qd = try dyn.quantile(p)
            let qt = try typ.quantile(p)
            #expect(abs(qd - qt) <= 1e-10)
        }

        if let meanDyn = dyn.mean, let meanTyp = typ.mean {
            #expect(abs(meanDyn - meanTyp) <= 1e-12)
        }
        if let varDyn = dyn.variance, let varTyp = typ.variance {
            #expect(abs(varDyn - varTyp) <= 1e-12)
        }
        if let skewDyn = dyn.skewness, let skewTyp = typ.skewness {
            #expect(abs(skewDyn - skewTyp) <= 1e-12)
        }
        if let kurtDyn = dyn.kurtosis, let kurtTyp = typ.kurtosis {
            #expect(abs(kurtDyn - kurtTyp) <= 1e-12)
        }
    }

    @Test("InverseGamma<Double> dynamic matches typed wrapper")
    func inverseGammaDynamicMatchesTypedDouble() throws {
        let shape: Double = 6.5
        let scale: Double = 1.3

        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "inverse_gamma",
            parameters: [
                "shape": shape,
                "scale": scale
            ]
        )
        let typ = try Distribution.InverseGamma<Double>(shape: shape, scale: scale)

        let xs: [Double] = [0.01, 0.05, 0.1, 0.5, 1.5, 4.0]
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

        let ps: [Double] = [1e-8, 1e-5, 1e-3, 0.1, 0.5, 0.9, 0.999]
        for p in ps {
            let xd = try dyn.quantile(p)
            let xt = try typ.quantile(p)
            let denom = max(abs(xt), 1.0)
            #expect(abs(xd - xt) / denom <= 1e-12)
        }

        let qs: [Double] = [1e-6, 1e-4, 0.05, 0.5]
        for q in qs {
            let xd = try dyn.quantileComplement(q)
            let xt = try typ.quantileComplement(q)
            let denom = max(abs(xt), 1.0)
            #expect(abs(xd - xt) / denom <= 1e-12)
        }

        let expectedMean = scale / (shape - 1)
        let expectedVariance = (scale * scale) / ((shape - 1) * (shape - 1) * (shape - 2))
        let expectedMode = scale / (shape + 1)
        let expectedSkewness = 4 * Foundation.sqrt(shape - 2) / (shape - 3)
        let expectedKurtosisExcess = (30 * shape - 66) / ((shape - 3) * (shape - 4))

        #expect(abs((dyn.mean ?? .nan) - expectedMean) <= 1e-12)
        #expect(abs((dyn.variance ?? .nan) - expectedVariance) <= 1e-12)
        #expect(abs((dyn.mode ?? .nan) - expectedMode) <= 1e-12)
        #expect(abs((dyn.skewness ?? .nan) - expectedSkewness) <= 1e-12)
        #expect(abs((dyn.kurtosisExcess ?? .nan) - expectedKurtosisExcess) <= 1e-12)
        if let kurtosis = dyn.kurtosis {
            #expect(abs(kurtosis - (expectedKurtosisExcess + 3)) <= 1e-12)
        } else {
            Issue.record("Expected kurtosis for inverse gamma with shape \(shape)")
        }
    }

    @Test("InverseNormal<Double> dynamic matches typed wrapper")
    func inverseNormalDynamicMatchesTypedDouble() throws {
        let mu: Double = 1.7
        let lambda: Double = 2.2

        let dyn = try Distribution.Dynamic<Double>(
            distributionName: "inverse_gaussian",
            parameters: [
                "mean": mu,
                "lambda": lambda
            ]
        )
        let typ = try Distribution.InverseNormal<Double>(mean: mu, shape: lambda)

        let xs: [Double] = [0.02, 0.1, 0.5, 1.0, 2.5]
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

        let ps: [Double] = [1e-6, 1e-3, 0.05, 0.5, 0.9, 0.999]
        for p in ps {
            let xd = try dyn.quantile(p)
            let xt = try typ.quantile(p)
            let denom = max(abs(xt), 1.0)
            #expect(abs(xd - xt) / denom <= 1e-12)
        }

        let qs: [Double] = [1e-5, 1e-3, 0.01, 0.25]
        for q in qs {
            let xd = try dyn.quantileComplement(q)
            let xt = try typ.quantileComplement(q)
            let denom = max(abs(xt), 1.0)
            #expect(abs(xd - xt) / denom <= 1e-12)
        }

        let expectedVariance = (mu * mu * mu) / lambda
        let expectedSkewness = 3 * Foundation.sqrt(mu / lambda)
        let expectedKurtosisExcess = 15 * mu / lambda
        let expectedMode = mu * (
            Foundation.sqrt(1 + (9 * mu * mu) / (4 * lambda * lambda))
            - (3 * mu) / (2 * lambda)
        )

        #expect(dyn.mean == mu)
        #expect(abs((dyn.variance ?? .nan) - expectedVariance) <= 1e-12)
        #expect(abs((dyn.mode ?? .nan) - expectedMode) <= 1e-12)
        #expect(abs((dyn.skewness ?? .nan) - expectedSkewness) <= 1e-12)
        #expect(abs((dyn.kurtosisExcess ?? .nan) - expectedKurtosisExcess) <= 1e-12)
        if let kurtosis = dyn.kurtosis {
            let expectedBoostKurtosis = expectedKurtosisExcess - 3
            #expect(abs(kurtosis - expectedBoostKurtosis) <= 1e-12)
        } else {
            Issue.record("Expected kurtosis for inverse normal with λ \(lambda)")
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

    @Test("Dynamic aliases for newly added continuous distributions")
    func dynamicAliasesForNewDistributions() throws {
        // Rayleigh with sigma alias
        do {
            let scale: Double = 1.3
            let dyn = try Distribution.Dynamic<Double>(distributionName: "rayleigh", parameters: ["sigma": scale])
            let typ = try Distribution.Rayleigh<Double>(scale: scale)
            let x: Double = 0.85
            #expect(abs(try dyn.pdf(x) - typ.pdf(x)) <= 1e-12)
            #expect(abs(try dyn.cdf(x) - typ.cdf(x)) <= 1e-12)
            let p: Double = 0.42
            #expect(abs(try dyn.quantile(p) - typ.quantile(p)) <= 1e-12)
        }

        // SAS Point5 with name alias and scale key "c"
        do {
            let location: Double = -0.3
            let scale: Double = 0.9
            let dyn = try Distribution.Dynamic<Double>(distributionName: "sas_point5", parameters: ["loc": location, "c": scale])
            let typ = try Distribution.SASPoint5<Double>(location: location, scale: scale)
            let x: Double = 1.1
            #expect(abs(try dyn.pdf(x) - typ.pdf(x)) <= 1e-12)
            #expect(abs(try dyn.cdf(x) - typ.cdf(x)) <= 1e-12)
        }

        // Skew-normal with alias and omega/alpha parameters
        do {
            let dyn = try Distribution.Dynamic<Double>(
                distributionName: "skewnormal",
                parameters: ["mu": 0.5, "omega": 1.2, "alpha": -3.0]
            )
            let typ = try Distribution.SkewNormal<Double>(location: 0.5, scale: 1.2, shape: -3)
            let x: Double = -0.4
            #expect(abs(try dyn.pdf(x) - typ.pdf(x)) <= 1e-12)
            #expect(abs(try dyn.cdf(x) - typ.cdf(x)) <= 1e-12)
            let p: Double = 0.73
            #expect(abs(try dyn.quantile(p) - typ.quantile(p)) <= 5e-11)
        }

        // Triangular with shorthand parameters a,b,c
        do {
            let dyn = try Distribution.Dynamic<Double>(
                distributionName: "triangle",
                parameters: ["a": -1.0, "c": 0.2, "b": 2.5]
            )
            let typ = try Distribution.Triangular<Double>(lower: -1.0, mode: 0.2, upper: 2.5)
            let x: Double = 0.9
            #expect(abs(try dyn.pdf(x) - typ.pdf(x)) <= 1e-12)
            #expect(abs(try dyn.cdf(x) - typ.cdf(x)) <= 1e-12)
        }

        // Uniform with rectangular alias and endpoints a/b
        do {
            let dyn = try Distribution.Dynamic<Double>(
                distributionName: "rectangular",
                parameters: ["a": -2.0, "b": 3.0]
            )
            let typ = try Distribution.Uniform<Double>(lower: -2.0, upper: 3.0)
            let x: Double = 1.1
            #expect(abs(try dyn.pdf(x) - typ.pdf(x)) <= 1e-12)
            #expect(abs(try dyn.cdf(x) - typ.cdf(x)) <= 1e-12)
            let p: Double = 0.6
            #expect(abs(try dyn.quantile(p) - typ.quantile(p)) <= 1e-12)
        }

        // Weibull with aliases for shape/scale
        do {
            let dyn = try Distribution.Dynamic<Double>(
                distributionName: "weibull",
                parameters: ["alpha": 1.7, "lambda": 0.8]
            )
            let typ = try Distribution.Weibull<Double>(shape: 1.7, scale: 0.8)
            let x: Double = 0.9
            #expect(abs(try dyn.pdf(x) - typ.pdf(x)) <= 1e-12)
            #expect(abs(try dyn.cdf(x) - typ.cdf(x)) <= 1e-12)
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
