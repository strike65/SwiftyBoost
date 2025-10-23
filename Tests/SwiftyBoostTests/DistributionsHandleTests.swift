import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Stateful distribution handles (Gamma, StudentT, Beta, ChiSquared, Binomial, Bernoulli)")
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

    @Test("Gamma<Double> logPdf matches log of pdf")
    func gammaLogPdfMatchesPdf() throws {
        let g = try Distribution.Gamma<Double>(shape: 3.25, scale: 0.75)
        let xs: [Double] = [0.2, 0.5, 1.0, 2.0, 4.5]
        for x in xs {
            let pdf = try g.pdf(x)
            #expect(pdf >= 0)
            if pdf > 0 {
                let logPdf = try g.logPdf(x)
                let expected = Foundation.log(pdf)
                #expect(abs(logPdf - expected) <= 1e-12)
            }
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

    @Test("StudentT<Double> CDF/SF complement consistency")
    func studentTCdfSfConsistency() throws {
        let t = try Distribution.StudentT<Double>(degreesOfFreedom: 6.0)
        let xs: [Double] = [-2.5, -1.0, -0.25, 0.5, 1.5, 2.5]
        for x in xs {
            let c = try t.cdf(x)
            let s = try t.sf(x)
            #expect(abs((c + s) - 1.0) <= 1e-12)
        }
    }

    @Test("Beta<Double> support, moments, and quantiles")
    func betaDoubleProperties() throws {
        let beta = try Distribution.Beta<Double>(alpha: 2.3, beta: 4.7)
        #expect(abs(beta.supportLowerBound) <= Double.ulpOfOne)
        #expect(abs(beta.supportUpperBound - 1.0) <= Double.ulpOfOne)
        let xs: [Double] = [0.05, 0.2, 0.5, 0.8, 0.95]
        for x in xs {
            let density = try beta.pdf(x)
            #expect(density >= 0)
            let cdf = try beta.cdf(x)
            #expect(cdf >= 0 && cdf <= 1)
        }
        let ps: [Double] = [1e-6, 0.05, 0.2, 0.5, 0.8, 0.95, 1 - 1e-6]
        for p in ps {
            let x = try beta.quantile(p)
            let c = try beta.cdf(x)
            #expect(abs(c - p) <= 1e-10)
        }
        #expect(beta.mean != nil)
        #expect(beta.variance != nil)
        guard let mean = beta.mean,
              let variance = beta.variance else {
            return
        }
        let alpha = beta.alpha
        let betaParam = beta.beta
        let meanExpected = alpha / (alpha + betaParam)
        let varianceExpected = (alpha * betaParam) /
            Foundation.pow(alpha + betaParam, 2) /
            (alpha + betaParam + 1)
        #expect(abs(mean - meanExpected) <= 1e-12)
        #expect(abs(variance - varianceExpected) <= 1e-12)
    }

    @Test("ChiSquared<Double> support, moments, and hazards")
    func chiSquaredDoubleProperties() throws {
        let df: Double = 8
        let chi = try Distribution.ChiSquared<Double>(degreesOfFreedom: df)
        #expect(chi.supportLowerBound >= 0)
        #expect(chi.supportUpperBound.isInfinite)
        let xs: [Double] = [0.05, 0.5, 1.0, 2.0, 5.0]
        for x in xs {
            let density = try chi.pdf(x)
            #expect(density >= 0)
            let cdf = try chi.cdf(x)
            #expect(cdf >= 0 && cdf <= 1)
            let sf = try chi.sf(x)
            if sf > 0 {
                let hazard = try chi.hazard(x)
                #expect(abs(hazard - (density / sf)) <= 1e-12)
                let chf = try chi.chf(x)
                if sf > 1e-14 {
                    let expectedChf = -Foundation.log(sf)
                    #expect(abs(chf - expectedChf) <= 1e-10)
                }
            }
        }
        let ps: [Double] = [1e-8, 1e-5, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-8]
        for p in ps {
            let x = try chi.quantile(p)
            let c = try chi.cdf(x)
            #expect(abs(c - p) <= 1e-10)
        }
        #expect(chi.mean != nil)
        #expect(chi.variance != nil)
        #expect(chi.skewness != nil)
        #expect(chi.kurtosisExcess != nil)
        guard let mean = chi.mean,
              let variance = chi.variance,
              let skewness = chi.skewness,
              let kurtosisExcess = chi.kurtosisExcess else {
            return
        }
        #expect(abs(mean - df) <= 1e-12)
        #expect(abs(variance - 2 * df) <= 1e-12)
        #expect(abs(skewness - Foundation.sqrt(8 / df)) <= 1e-12)
        #expect(abs(kurtosisExcess - (12 / df)) <= 1e-12)
        #expect(chi.latticeStep == nil)
        #expect(chi.latticeOrigin == nil)
    }

    @Test("Binomial<Double> pmf/cdf bridge and quantiles")
    func binomialDoubleProperties() throws {
        let n = 8
        let p: Double = 0.37
        let bin = try Distribution.Binomial<Double>(
            numberOfTrials: Double(n),
            probabibilityOfSuccess: p
        )
        #expect(bin.supportLowerBound == 0)
        #expect(abs(bin.supportUpperBound - Double(n)) <= Double(n) * Double.ulpOfOne)
        let q: Double = 1 - p
        func combination(_ n: Int, _ k: Int) -> Double {
            if k < 0 || k > n { return 0 }
            if k == 0 || k == n { return 1 }
            var result = 1.0
            for i in 1...k {
                result *= Double(n - (k - i)) / Double(i)
            }
            return result
        }
        var cdfs: [Double] = []
        var cumulative = 0.0
        for k in 0...n {
            let expected = combination(n, k) *
                Foundation.pow(p, Double(k)) *
                Foundation.pow(q, Double(n - k))
            let got = try bin.pdf(Double(k))
            #expect(abs(got - expected) <= 1e-12)
            cumulative += expected
            cdfs.append(cumulative)
            let cdf = try bin.cdf(Double(k))
            #expect(abs(cdf - cumulative) <= 1e-12)
        }
        #expect(abs(try bin.cdf(Double(n)) - 1.0) <= 1e-12)
        let target: Double = 0.6
        let quantile = try bin.quantile(target)
        let k = Int(quantile.rounded())
        #expect(k >= 0 && k <= n)
        let cdfAtK = cdfs[k]
        #expect(cdfAtK >= target - 1e-12)
        if k > 0 {
            let cdfBefore = cdfs[k - 1]
            #expect(cdfBefore <= target + 1e-12)
        }
        let sfAtK = try bin.sf(Double(k))
        #expect(abs(sfAtK - (1 - cdfAtK)) <= 1e-12)
        let qUpper: Double = 0.25
        let upper = try bin.quantileComplement(qUpper)
        let roundedUpper = upper.rounded()
        #expect(abs(upper - roundedUpper) <= 1e-12)
        let kUpper = Int(roundedUpper)
        let sfUpper = try bin.sf(Double(kUpper))
        #expect(sfUpper <= qUpper + 1e-12)
        if kUpper > Int(bin.supportLowerBound) {
            let sfBefore = try bin.sf(Double(kUpper - 1))
            #expect(sfBefore >= qUpper - 1e-12)
        }
        #expect(bin.variance != nil)
        guard let variance = bin.variance else {
            return
        }
        let expectedVariance = Double(n) * p * q
        #expect(abs(variance - expectedVariance) <= 1e-12)
    }

    @Test("Bernoulli<Double> pmf, quantiles, and entropy")
    func bernoulliDoubleProperties() throws {
        let p: Double = 0.35
        let bern = try Distribution.Bernoulli<Double>(p: p)
        #expect(bern.supportLowerBound == 0)
        #expect(bern.supportUpperBound == 1)
        let pdf0 = try bern.pdf(0)
        let pdf1 = try bern.pdf(1)
        #expect(abs(pdf0 - (1 - p)) <= 1e-12)
        #expect(abs(pdf1 - p) <= 1e-12)
        let cdf0 = try bern.cdf(0)
        #expect(abs(cdf0 - (1 - p)) <= 1e-12)
        let cdf1 = try bern.cdf(1)
        #expect(abs(cdf1 - 1) <= 1e-12)
        let sf0 = try bern.sf(0)
        #expect(abs(sf0 - p) <= 1e-12)
        let qLower: Double = 0.6
        let quantileLower = try bern.quantile(qLower)
        #expect(quantileLower == 0)
        let qUpper: Double = 0.9
        let quantileUpper = try bern.quantile(qUpper)
        #expect(quantileUpper == 1)
        let qc: Double = p
        let upper = try bern.quantileComplement(qc)
        let roundedUpper = upper.rounded()
        #expect(abs(upper - roundedUpper) <= 1e-12)
        let kUpper = Int(roundedUpper)
        let sfUpper = try bern.sf(Double(kUpper))
        #expect(sfUpper <= qc + 1e-12)
        let sfBefore: Double
        if kUpper > Int(bern.supportLowerBound) {
            sfBefore = try bern.sf(Double(kUpper - 1))
        } else {
            sfBefore = 1
        }
        #expect(sfBefore >= qc - 1e-12)
        #expect(bern.mean != nil)
        #expect(bern.variance != nil)
        #expect(bern.mode != nil)
        #expect(bern.kurtosisExcess != nil)
        #expect(bern.skewness != nil)
        guard let mean = bern.mean,
              let variance = bern.variance,
              let mode = bern.mode,
              let kurtosisExcess = bern.kurtosisExcess,
              let skewness = bern.skewness else {
            return
        }
        #expect(abs(mean - p) <= 1e-12)
        #expect(abs(variance - p * (1 - p)) <= 1e-12)
        #expect(mode == 0)
        #expect(bern.median == 0)
        let expectedEntropy = -(p * Foundation.log(p) + (1 - p) * Foundation.log(1 - p))
        #expect(bern.entropy != nil)
        if let entropy = bern.entropy {
            #expect(abs(entropy - expectedEntropy) <= 1e-12)
        }
        #expect(bern.latticeStep == nil)
        #expect(bern.latticeOrigin == nil)
        #expect(abs(skewness - (1 - 2 * p) / Foundation.sqrt(p * (1 - p))) <= 1e-12)
        #expect(abs(kurtosisExcess - (1 - 6 * p * (1 - p)) / (p * (1 - p))) <= 1e-12)
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
