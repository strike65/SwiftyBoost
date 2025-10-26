import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Empirical Distribution – Discrete & Continuous Detection")
struct EmpiricalDistributionTests {
    @Test("discrete sample uses lattice detection and smoothed entropy")
    func discreteSampleProperties() throws {
        let samples: [Double] = [1, 2, 2, 4]
        let dist = try Distribution.Empirical(samples: samples)

        #expect(dist.isDiscrete)
        #expect(dist.latticeStep == 1)
        #expect(dist.latticeOrigin == 1)

        let pdf2 = try dist.pdf(2)
        #expect(abs(pdf2 - 0.4545454545) <= 1e-6)
        #expect(try dist.pdf(3) == 0)

        let cdfAt2 = try dist.cdf(2)
        #expect(abs(cdfAt2 - 0.7272727272) <= 1e-6)
        #expect(try dist.sf(2) == 1 - cdfAt2)

        #expect(abs((try dist.hazard(2)) - (pdf2 / (1 - cdfAt2))) <= 1e-10)
        #expect(abs((try dist.chf(2)) - (-Foundation.log(1 - cdfAt2))) <= 1e-10)

        #expect(try dist.quantile(0) == 1)
        #expect(try dist.quantile(0.5) == 2)
        #expect(try dist.quantile(0.75) == 4)
        #expect(try dist.quantileComplement(0.25) == 4)

        #expect(dist.mode == 2)
        #expect(dist.median == 2)

        let entropy = try #require(dist.entropy)
        #expect(abs(entropy - 1.316140) <= 1e-3)

        let entropyEstimate = try dist.entropyEstimate(bootstrapSamples: 5, confidenceLevel: 0.8)
        #expect(entropyEstimate.confidenceInterval != nil)
    }

    @Test("discrete KL divergence with smoothing is finite")
    func discreteKLEstimate() throws {
        let a = try Distribution.Empirical(samples: [0.0, 0.0, 1.0, 2.0])
        let b = try Distribution.Empirical(samples: [0.0, 1.0, 2.0, 2.0])

        let maybeKL = try a.klDivergence(relativeTo: b)
        let kl = try #require(maybeKL)
        #expect(kl.isFinite)

        let maybeBootstrap = try a.klDivergenceEstimate(relativeTo: b, bootstrapSamples: 5, confidenceLevel: 0.8)
        let bootstrap = try #require(maybeBootstrap)
        #expect(bootstrap.confidenceInterval != nil)
    }

    @Test("continuous sample falls back to KDE/KNN estimators")
    func continuousSampleProperties() throws {
        let samples: [Double] = [-1.2, -0.6, -0.2, 0.0, 0.35, 0.8, 1.2, 1.6, 2.1]
        let dist = try Distribution.Empirical(samples: samples)

        #expect(!dist.isDiscrete)
        #expect(dist.latticeStep == nil)
        #expect(dist.latticeOrigin == nil)

        let pdf0 = try dist.pdf(0)
        #expect(pdf0 > 0)
        let entropy = dist.entropy
        #expect(entropy != nil)

        let mode = try #require(dist.mode)
        #expect(mode >= samples.min()! && mode <= samples.max()!)

        let q25 = try dist.quantile(0.25)
        let q75 = try dist.quantile(0.75)
        #expect(q25 <= q75)

        let entropyEstimate = try dist.entropyEstimate(estimator: .knn(k: 2), bootstrapSamples: 5, confidenceLevel: 0.8)
        #expect(entropyEstimate.confidenceInterval != nil)
    }

    @Test("continuous KL divergence via KNN remains finite")
    func continuousKLEstimate() throws {
        let p = try Distribution.Empirical(samples: [-1.0, -0.5, -0.2, 0.0, 0.3, 0.8, 1.1])
        let q = try Distribution.Empirical(samples: [-0.8, -0.4, -0.1, 0.2, 0.5, 1.0, 1.4])

        let maybeContinuousKL = try p.klDivergence(relativeTo: q, estimator: .knn(k: 2))
        let kl = try #require(maybeContinuousKL)
        #expect(kl.isFinite)

        let maybeKLBootstrap = try p.klDivergenceEstimate(
            relativeTo: q,
            estimator: .knn(k: 2),
            bootstrapSamples: 5,
            confidenceLevel: 0.8
        )
        let bootstrap = try #require(maybeKLBootstrap)
        #expect(bootstrap.confidenceInterval != nil)
    }

    @Test("multimodality detection via Hartigan dip heuristic")
    func multimodalityHeuristic() throws {
        let unimodal = try Distribution.Empirical(samples: stride(from: -1.0, through: 1.0, by: 0.1).map { $0 })
        #expect(!unimodal.isLikelyMultimodal)

        let multimodal = try Distribution.Empirical(samples: [-3.0, -2.6, -2.5, -0.2, 0, 0.2, 2.4, 2.8, 3.2])
        #expect(multimodal.isLikelyMultimodal)
    }
}
