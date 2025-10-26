import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Empirical Distribution – Discrete & Continuous Detection")
struct EmpiricalDistributionTests {
    private let discreteFile = "empirical_discrete_n200.csv"
    private let continuousFile = "empirical_continuous_n200.csv"

    @Test("discrete sample uses lattice detection and smoothed entropy")
    func discreteSampleProperties() throws {
        let samples = try loadSamples(discreteFile)
        let dist = try Distribution.Empirical(samples: samples)

        #expect(dist.isDiscrete)
        #expect(dist.latticeStep == 1)
        #expect(dist.latticeOrigin == -3)

        // Expected smoothed probabilities with alpha = 0.5
        let counts = frequencyTable(samples)
        let uniqueValues = Array(counts.keys).sorted()
        let alpha = 0.5
        let total = Double(samples.count)
        let uniqueCount = Double(uniqueValues.count)
        let denominator = total + alpha * uniqueCount

        func smoothedProbability(_ value: Int) -> Double {
            let count = Double(counts[value] ?? 0)
            return (count + alpha) / denominator
        }

        let expectedPDF2 = smoothedProbability(2)
        let pdf2 = try dist.pdf(2)
        #expect(abs(pdf2 - expectedPDF2) <= 1e-12)

        let expectedCDF2 = uniqueValues
            .filter { $0 <= 2 }
            .reduce(0.0) { $0 + smoothedProbability($1) }
        let cdf2 = try dist.cdf(2)
        #expect(abs(cdf2 - expectedCDF2) <= 1e-12)
        #expect(abs((try dist.sf(2)) - (1 - cdf2)) <= 1e-12)

        let hazard2 = try dist.hazard(2)
        let expectedHazard2 = expectedPDF2 / max(1 - expectedCDF2, Double.leastNonzeroMagnitude)
        #expect(abs(hazard2 - expectedHazard2) <= 1e-9)

        let chf2 = try dist.chf(2)
        #expect(abs(chf2 - (-log(max(1 - expectedCDF2, Double.leastNonzeroMagnitude)))) <= 1e-9)

        let quantile0 = try dist.quantile(0)
        #expect(quantile0 == -3)
        let complement1 = try dist.quantileComplement(0)
        #expect(complement1 == dist.supportUpperBound)

        let expectedMean = uniqueValues.reduce(0.0) { $0 + Double($1) * smoothedProbability($1) }
        #expect(abs((dist.mean ?? .nan) - expectedMean) <= 1e-9)

        let entropy = try #require(dist.entropy)
        let probabilities = uniqueValues.map(smoothedProbability)
        let entropyCore = -probabilities.reduce(0.0) { sum, p in
            guard p > 0 else { return sum }
            return sum + p * log(p)
        }
        let correction = (uniqueCount - 1) / (2 * total)
        let expectedEntropy = entropyCore + correction
        #expect(abs(entropy - expectedEntropy) <= 1e-9)

        let entropyEstimate = try dist.entropyEstimate(bootstrapSamples: 5, confidenceLevel: 0.8)
        #expect(entropyEstimate.confidenceInterval != nil)
    }

    @Test("discrete KL divergence with smoothing is finite")
    func discreteKLEstimate() throws {
        let samples = try loadSamples(discreteFile)
        let otherSamples = samples.map { value -> Double in
            if value >= 2 { return min(3, value + 1) }
            if value <= -2 { return max(-3, value + 1) }
            return value
        }

        let p = try Distribution.Empirical(samples: samples)
        let q = try Distribution.Empirical(samples: otherSamples)

        let maybeKL = try p.klDivergence(relativeTo: q)
        let kl = try #require(maybeKL)
        #expect(kl.isFinite)

        let expected = discreteKLEstimate(samplesP: samples, samplesQ: otherSamples)
        #expect(abs(Double(kl) - expected) <= 1e-6)

        let maybeBootstrap = try p.klDivergenceEstimate(relativeTo: q, bootstrapSamples: 5, confidenceLevel: 0.8)
        let bootstrap = try #require(maybeBootstrap)
        #expect(bootstrap.confidenceInterval != nil)
    }

    @Test("continuous sample falls back to KDE/KNN estimators")
    func continuousSampleProperties() throws {
        let samples = try loadSamples(continuousFile)
        let dist = try Distribution.Empirical(samples: samples)

        #expect(!dist.isDiscrete)
        #expect(dist.latticeStep == nil)
        #expect(dist.latticeOrigin == nil)

        let pdf0 = try dist.pdf(0)
        #expect(pdf0 > 0)

        let q50 = try dist.quantile(0.5)
        let cdfAtQ50 = try dist.cdf(q50)
        #expect(cdfAtQ50 >= 0.5 - 1e-6 && cdfAtQ50 <= 0.5 + 1e-6)

        #expect(dist.entropy != nil)

        let entropyEstimate = try dist.entropyEstimate(estimator: .knn(k: 3), bootstrapSamples: 5, confidenceLevel: 0.8)
        #expect(entropyEstimate.confidenceInterval != nil)
    }

    @Test("continuous KL divergence via KNN remains finite")
    func continuousKLEstimate() throws {
        let samples = try loadSamples(continuousFile)
        let shifted = samples.map { $0 + 0.2 }

        let p = try Distribution.Empirical(samples: samples)
        let q = try Distribution.Empirical(samples: shifted)

        let maybeContinuousKL = try p.klDivergence(relativeTo: q, estimator: .knn(k: 3))
        let kl = try #require(maybeContinuousKL)
        #expect(kl.isFinite)

        let maybeKLBootstrap = try p.klDivergenceEstimate(
            relativeTo: q,
            estimator: .knn(k: 3),
            bootstrapSamples: 5,
            confidenceLevel: 0.8
        )
        let bootstrap = try #require(maybeKLBootstrap)
        #expect(bootstrap.confidenceInterval != nil)
    }

    @Test("multimodality detection heuristic reacts to bimodal data")
    func multimodalityHeuristic() throws {
        let unimodal = try Distribution.Empirical(samples: stride(from: -1.0, through: 1.0, by: 0.1).map { $0 })
        #expect(!unimodal.isLikelyMultimodal)

        let continuous = try loadSamples(continuousFile)
        let bimodal = try Distribution.Empirical(samples: continuous)
        #expect(bimodal.isLikelyMultimodal)
    }

    // MARK: Helpers

    private func loadSamples(_ fileName: String) throws -> [Double] {
        let testsDirectory = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent() // EmpiricalDistributionTests.swift
            .appendingPathComponent("TestData")
        let url = testsDirectory.appendingPathComponent(fileName)
        let data = try String(contentsOf: url)
        return data
            .split(whereSeparator: \.isNewline)
            .dropFirst() // header
            .compactMap { Double($0.trimmingCharacters(in: .whitespacesAndNewlines)) }
    }

    private func frequencyTable(_ samples: [Double]) -> [Int: Int] {
        var table: [Int: Int] = [:]
        for value in samples {
            table[Int(value), default: 0] += 1
        }
        return table
    }

    private func discreteKLEstimate(samplesP: [Double], samplesQ: [Double]) -> Double {
        let alpha = 0.5

        func smoothedDistribution(_ samples: [Double]) -> (counts: [Int: Int], denominator: Double, unique: Set<Int>) {
            let counts = frequencyTable(samples)
            let unique = Set(counts.keys)
            let denominator = Double(samples.count) + alpha * Double(unique.count)
            return (counts, denominator, unique)
        }

        let distP = smoothedDistribution(samplesP)
        let distQ = smoothedDistribution(samplesQ)
        let support = distP.unique.union(distQ.unique).sorted()

        func probability(_ value: Int, dist: (counts: [Int: Int], denominator: Double, unique: Set<Int>)) -> Double {
            let count = Double(dist.counts[value] ?? 0)
            return (count + alpha) / dist.denominator
        }

        var kl: Double = 0
        for value in support {
            let p = probability(value, dist: distP)
            let q = probability(value, dist: distQ)
            kl += p * (log(p) - log(q))
        }
        return kl
    }
}
