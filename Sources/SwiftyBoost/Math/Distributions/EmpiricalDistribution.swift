//
//  EmpiricalDistribution.swift
//  SwiftyBoost
//
//  Created by Codex Agent (GPT-5) for Volker Thieme.
//

import SwiftyBoostPrelude
import Foundation

extension Distribution {
    /// Empirical distribution backed by a finite set of observations.
    ///
    /// The implementation distinguishes between lattice/discrete and continuous samples,
    /// exposes the full ``DistributionProtocol`` surface, and provides bootstrap-based
    /// uncertainty estimates for entropy and KL divergence.
    public struct Empirical<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        // MARK: Nested types

        /// Bootstrap confidence-interval construction.
        /// Bootstrap resampling strategy used when constructing confidence intervals.
        ///
        /// - ``BootstrapMethod/percentile``: Uses the percentile method and returns quantiles of the empirical
        ///   bootstrap distribution without jackknife acceleration.
        /// - ``BootstrapMethod/bca``: Applies bias-corrected and accelerated adjustments when jackknife estimates
        ///   are available. Falls back to percentile behaviour when jackknife data cannot be produced.
        public enum BootstrapMethod: Sendable {
            case percentile
            case bca
        }

        /// Result bundle containing a point estimate and (optionally) a confidence interval.
        public struct Estimate<Value: Sendable>: Sendable {
            /// Point estimate returned by the requested statistic.
            public let value: Value
            /// Optional confidence interval covering the `confidenceLevel` supplied to the bootstrap call.
            public let confidenceInterval: (lower: Value, upper: Value)?
            /// Number of bootstrap replicates used to form the estimate.
            public let replicates: Int
            /// Resampling strategy that produced the interval.
            public let method: BootstrapMethod

            @inlinable
            public init(
                value: Value,
                confidenceInterval: (lower: Value, upper: Value)?,
                replicates: Int,
                method: BootstrapMethod
            ) {
                self.value = value
                self.confidenceInterval = confidenceInterval
                self.replicates = replicates
                self.method = method
            }
        }

        /// Density estimator preference for differential entropy / KL divergence.
        public enum DensityEstimator: Sendable {
            /// Choose between KNN and KDE automatically based on sample characteristics.
            case automatic
            /// Explicitly request the k-nearest neighbour estimator.
            case knn(k: Int)
            /// Use a Gaussian kernel density estimate. When `bandwidth` is nil a plug-in bandwidth is computed.
            case kdeGaussian(bandwidth: T?)
        }

        // MARK: Stored properties

        private let samplesSorted: [T]
        private let samplesDouble: [Double]
        private let uniqueValues: [T]
        private let uniqueDoubles: [Double]
        private let counts: [Int]
        private let totalCount: Int
        private let smoothingAlpha: T = 0.5
        private let knnDefaultK: Int = 3

        private let meanCache: T
        private let varianceCache: T
        private let skewnessCache: T?
        private let kurtosisCache: T?
        private let modeCache: T?
        private let medianCache: T
        private let entropyCache: T?
        private let isDiscreteDistribution: Bool
        private let latticeStepValue: T?
        private let latticeOriginValue: T?
        private let rangeLower: T
        private let rangeUpper: T
        private let likelyMultimodal: Bool

        // MARK: Initialisation

        /// Creates an empirical distribution from a finite sample.
        ///
        /// The constructor validates and sorts the supplied observations, automatically detects
        /// whether the data is best treated as discrete (lattice) or continuous, and precomputes
        /// summary statistics as well as smoothed probability masses. Repeated observations are
        /// collapsed to unique support points for efficiency, while the original multiplicities
        /// remain available for bootstrap resampling.
        ///
        /// - Parameter samples: Sequence of finite observations. The sequence is consumed eagerly.
        /// - Throws:
        ///   - ``DistributionError/invalidCombination(message:value:)`` when the sequence is empty.
        ///   - ``DistributionError/parameterNotFinite(name:value:)`` if any sample is NaN or ±infinity.
        public init<S: Sequence>(samples: S) throws where S.Element == T {
            var collected: [T] = []
            collected.reserveCapacity(64)
            for value in samples {
                guard value.isFinite else {
                    throw DistributionError<T>.parameterNotFinite(name: "samples", value: value)
                }
                collected.append(value)
            }
            guard !collected.isEmpty else {
                throw DistributionError<T>.invalidCombination(
                    message: "Empirical distribution requires at least one observation.",
                    value: nil
                )
            }

            collected.sort()
            self.samplesSorted = collected
            self.samplesDouble = collected.map { Double($0) }
            self.totalCount = collected.count
            self.rangeLower = collected.first!
            self.rangeUpper = collected.last!

            // Collapse to unique support points with multiplicities.
            var values: [T] = []
            var doubles: [Double] = []
            var multiplicities: [Int] = []
            values.reserveCapacity(collected.count)
            multiplicities.reserveCapacity(collected.count)

            var iterator = collected.makeIterator()
            if let first = iterator.next() {
                var current = first
                var currentCount = 1
                while let next = iterator.next() {
                    if next == current {
                        currentCount += 1
                    } else {
                        values.append(current)
                        doubles.append(Double(current))
                        multiplicities.append(currentCount)
                        current = next
                        currentCount = 1
                    }
                }
                values.append(current)
                doubles.append(Double(current))
                multiplicities.append(currentCount)
            }
            self.uniqueValues = values
            self.uniqueDoubles = doubles
            self.counts = multiplicities

            let detection = Self.detectDiscreteSupport(values: values, totalCount: collected.count, toleranceScale: T(1e-6))
            self.isDiscreteDistribution = detection.isDiscrete
            self.latticeStepValue = detection.latticeStep
            self.latticeOriginValue = detection.latticeOrigin

            let probabilities = Self.normalisedFrequencies(
                counts: multiplicities,
                total: collected.count,
                alpha: smoothingAlpha
            )

            // Mean/variance/skewness/kurtosis
            var mean: T = 0
            for idx in values.indices {
                mean += values[idx] * probabilities[idx]
            }
            self.meanCache = mean

            var centeredSecond: T = 0
            var centeredThird: T = 0
            var centeredFourth: T = 0
            for idx in values.indices {
                let diff = values[idx] - mean
                let prob = probabilities[idx]
                centeredSecond += prob * diff * diff
                centeredThird += prob * diff * diff * diff
                centeredFourth += prob * diff * diff * diff * diff
            }
            self.varianceCache = centeredSecond
            if centeredSecond > 0 {
                let sigma = centeredSecond.squareRoot()
                self.skewnessCache = centeredThird / (sigma * sigma * sigma)
                self.kurtosisCache = centeredFourth / (sigma * sigma * sigma * sigma)
            } else {
                self.skewnessCache = nil
                self.kurtosisCache = nil
            }

            // Mode: discrete -> highest smoothed probability; continuous -> KDE max.
            if detection.isDiscrete {
                var bestValue = values.first
                var bestProbability = probabilities.first ?? 1
                for idx in values.indices.dropFirst() {
                    let prob = probabilities[idx]
                    if prob > bestProbability || (prob == bestProbability && values[idx] < bestValue!) {
                        bestProbability = prob
                        bestValue = values[idx]
                    }
                }
                self.modeCache = bestValue
            } else {
                self.modeCache = Self.estimateContinuousMode(samples: self.samplesDouble)
            }

            // Median via quantile 0.5
            self.medianCache = Self.quantile(samples: collected, probabilities: probabilities, p: 0.5)

            if detection.isDiscrete {
                self.entropyCache = Self.discreteEntropy(probabilities: probabilities, uniqueCount: values.count, sampleSize: collected.count)
            } else if let knnEntropy = Self.continuousEntropyKNN(samples: self.samplesDouble, k: knnDefaultK) {
                self.entropyCache = knnEntropy
            } else {
                self.entropyCache = Self.continuousEntropyKDE(samples: self.samplesDouble, bandwidth: nil)
            }

            self.likelyMultimodal = Self.detectMultimodality(
                samples: self.samplesDouble,
                isDiscrete: detection.isDiscrete
            )
        }

        // MARK: DistributionProtocol – Support

        /// Lowest observed value in the sample (minimum of the support).
        public var supportLowerBound: T { rangeLower }
        /// Highest observed value in the sample (maximum of the support).
        public var supportUpperBound: T { rangeUpper }
        /// Convenience tuple exposing both support endpoints.
        public var range: (lower: T, upper: T) { (rangeLower, rangeUpper) }

        // MARK: - Core functions

        /// Evaluates the empirical PDF/PMF or a KDE density depending on the detected support type.
        ///
        /// - Parameter x: Query position.
        /// - Returns: Smoothed probability mass for discrete samples or a continuous density estimate.
        /// - Throws: ``DistributionError/parameterNotFinite(name:value:)`` when `x` is not finite.
        public func pdf(_ x: T) throws -> T {
            guard x.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "x", value: x)
            }
            if isDiscreteDistribution {
                guard let index = Self.binarySearch(values: uniqueValues, target: x) else {
                    return 0
                }
                let probs = Self.normalisedFrequencies(
                    counts: counts,
                    total: totalCount,
                    alpha: smoothingAlpha
                )
                return probs[index]
            }
            let density = Self.kdeDensityGaussian(
                samples: samplesDouble,
                x: Double(x),
                bandwidth: nil
            )
            return T(density)
        }

        /// Natural logarithm of ``pdf(_:)``. Returns `-infinity` for zero mass.
        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            if p > 0 {
                return T.log(p)
            }
            return -.infinity
        }

        /// Empirical CDF evaluated at `x`.
        ///
        /// Discrete samples accumulate smoothed point probabilities; continuous samples use the
        /// empirical cumulative function based on order statistics.
        public func cdf(_ x: T) throws -> T {
            guard x.isFinite else {
                if x == T.infinity { return 1 }
                if x == -T.infinity { return 0 }
                throw DistributionError<T>.parameterNotFinite(name: "x", value: x)
            }
            if x < rangeLower { return 0 }
            if x >= rangeUpper { return 1 }
            if isDiscreteDistribution {
                let probs = Self.normalisedFrequencies(
                    counts: counts,
                    total: totalCount,
                    alpha: smoothingAlpha
                )
                var cumulative: T = 0
                for idx in uniqueValues.indices {
                    if uniqueValues[idx] > x { break }
                    cumulative += probs[idx]
                }
                return min(max(cumulative, 0), 1)
            } else {
                // Continuous CDF via empirical cumulative (midpoint correction).
                let idx = Self.upperBound(in: samplesSorted, value: x)
                return T(idx) / T(totalCount)
            }
        }

        /// Survival function `1 - CDF(x)` with guarding against negative underflow.
        public func sf(_ x: T) throws -> T {
            let c = try cdf(x)
            return max(0, 1 - c)
        }

        /// Hazard rate. For discrete samples this is defined using the smoothed PMF.
        public func hazard(_ x: T) throws -> T {
            let density = try pdf(x)
            if density == 0 { return 0 }
            let survival = try sf(x)
            if survival > 0 {
                return density / survival
            }
            return .infinity
        }

        /// Cumulative hazard `-log(SF(x))` guarded for degenerate tails.
        public func chf(_ x: T) throws -> T {
            let survival = try sf(x)
            if survival <= 0 { return .infinity }
            return -T.log(survival)
        }

        // MARK: - Inverses

        /// Lower-tail quantile. For discrete samples the first point with cumulative mass ≥ `p` is returned.
        public func quantile(_ p: T) throws -> T {
            guard p >= 0, p <= 1 else {
                throw DistributionError<T>.parameterOutOfRange(name: "p", min: 0, max: 1)
            }
            if p == 0 { return rangeLower }
            if p == 1 { return rangeUpper }
            if isDiscreteDistribution {
                let probs = Self.normalisedFrequencies(
                    counts: counts,
                    total: totalCount,
                    alpha: smoothingAlpha
                )
                var cumulative: T = 0
                for idx in uniqueValues.indices {
                    cumulative += probs[idx]
                    if cumulative >= p {
                        return uniqueValues[idx]
                    }
                }
                return rangeUpper
            } else {
                let position = T(totalCount - 1) * p
                let lowerIndex = Int(position.rounded(.down))
                let upperIndex = min(totalCount - 1, lowerIndex + 1)
                if lowerIndex == upperIndex {
                    return samplesSorted[lowerIndex]
                }
                let fraction = position - T(lowerIndex)
                return samplesSorted[lowerIndex] * (1 - fraction) + samplesSorted[upperIndex] * fraction
            }
        }

        /// Upper-tail quantile, equivalent to `quantile(1 - q)`.
        public func quantileComplement(_ q: T) throws -> T {
            guard q >= 0, q <= 1 else {
                throw DistributionError<T>.parameterOutOfRange(name: "q", min: 0, max: 1)
            }
            return try quantile(1 - q)
        }

        // MARK: - Summary statistics

        /// Sample mean calculated from the smoothed probability masses.
        public var mean: T? { meanCache }
        /// Sample variance based on smoothed masses; zero when all observations coincide.
        public var variance: T? { varianceCache }
        /// Mode of the distribution. Discrete samples report the most frequent value, continuous samples approximate via KDE.
        public var mode: T? { modeCache }
        /// Empirical median (0.5 quantile).
        public var median: T { medianCache }
        /// Smoothed sample skewness or `nil` when the variance is zero.
        public var skewness: T? { skewnessCache }
        /// Smoothed sample kurtosis or `nil` when the variance is zero.
        public var kurtosis: T? { kurtosisCache }
        public var kurtosisExcess: T? {
            guard let kurtosis = kurtosisCache else { return nil }
            return kurtosis - 3
        }

        /// Lattice spacing detected for discrete samples, otherwise `nil`.
        public var latticeStep: T? { latticeStepValue }
        /// Lattice origin for discrete samples, otherwise `nil`.
        public var latticeOrigin: T? { latticeOriginValue }

        /// Shannon entropy estimate (Miller–Madow corrected for discrete samples, KNN/KDE otherwise).
        public var entropy: T? { entropyCache }
        /// Flag indicating whether the sample was classified as discrete.
        public var isDiscrete: Bool { isDiscreteDistribution }

        /// Indicates whether the empirical sample exhibits multimodality according to KDE/local-max scans.
        public var isLikelyMultimodal: Bool { likelyMultimodal }

        // MARK: - KL divergence

        /// KL divergence between two empirical distributions using an automatic estimator.
        public func klDivergence(
            relativeTo other: Self,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            try klDivergence(relativeTo: other, estimator: .automatic)
        }

        /// KL divergence estimate using the requested density estimator.
        ///
        /// - Parameters:
        ///   - other: Reference empirical distribution.
        ///   - estimator: Density estimator selection (.knn or .kdeGaussian). `.automatic` chooses KNN when both
        ///     samples appear continuous, otherwise smoothed discrete estimator.
        /// - Returns: Estimated divergence (nats) or `nil` when undefined.
        public func klDivergence(
            relativeTo other: Self,
            estimator: DensityEstimator
        ) throws -> T? {
            if isDiscreteDistribution || other.isDiscreteDistribution {
                let smoothing = smoothingAlpha
                let merged = Self.unionSupport(left: self, right: other)
                var sum: Double = 0
                for entry in merged {
                    let p = Double(Self.smoothedProbability(
                        value: entry,
                        distribution: self,
                        alpha: smoothing
                    ))
                    let q = Double(Self.smoothedProbability(
                        value: entry,
                        distribution: other,
                        alpha: smoothing
                    ))
                    if p == 0 { continue }
                    sum += p * (log(p) - log(q))
                }
                return T(sum)
            }

            let chosenEstimator: DensityEstimator
            switch estimator {
            case .automatic:
                chosenEstimator = .knn(k: min(knnDefaultK, max(1, samplesSorted.count / 4)))
            default:
                chosenEstimator = estimator
            }

            switch chosenEstimator {
            case .knn(let k):
                return Self.klDivergenceKNN(
                    samplesP: samplesDouble,
                    samplesQ: other.samplesDouble,
                    k: max(1, k)
                ).map { T($0) }
            case .kdeGaussian(let bandwidth):
                return Self.klDivergenceKDE(
                    samplesP: samplesDouble,
                    samplesQ: other.samplesDouble,
                    bandwidth: bandwidth.map { Double($0) }
                ).map { T($0) }
            case .automatic:
                return nil // unreachable
            }
        }

        // MARK: - Bootstrap estimates

        /// Bootstrap estimate (point + confidence interval) for the entropy.
        ///
        /// - Parameters:
        ///   - estimator: Density estimator used for each replicate. `.automatic` mirrors ``entropy`` defaults.
        ///   - bootstrapSamples: Number of resampled replicates. Higher values tighten intervals at the cost of runtime.
        ///   - confidenceLevel: Central mass retained by the interval. Values near 1 request wider intervals.
        ///   - method: Resampling interval strategy (percentile or BCa).
        /// - Returns: Point estimate with optional confidence bounds. Bounds are `nil` when the interval could not be formed.
        public func entropyEstimate(
            estimator: DensityEstimator = .automatic,
            bootstrapSamples: Int = 200,
            confidenceLevel: T = 0.95,
            method: BootstrapMethod = .percentile
        ) throws -> Estimate<T> {
            let point: T
            switch estimator {
            case .automatic:
                point = entropy ?? 0
            case .knn(let k):
                if let value = Self.continuousEntropyKNN(samples: samplesDouble, k: max(1, k)) {
                    point = value
                } else {
                    point = Self.continuousEntropyKDE(samples: samplesDouble, bandwidth: nil)
                }
            case .kdeGaussian(let bandwidth):
                point = Self.continuousEntropyKDE(samples: samplesDouble, bandwidth: bandwidth.map { Double($0) })
            }

            let interval = try Self.bootstrapOneSample(
                data: samplesSorted,
                estimator: { sample -> T in
                    if isDiscreteDistribution {
                        return T(Self.discreteEntropy(
                            probabilities: Self.normalisedFrequencies(
                                counts: Self.multiplicities(for: sample),
                                total: sample.count,
                                alpha: smoothingAlpha
                            ),
                            uniqueCount: Self.uniqueValues(from: sample).count,
                            sampleSize: sample.count
                        ))
                    } else {
                        let doubles = sample.map { Double($0) }
                        switch estimator {
                        case .automatic, .knn:
                            if let value = Self.continuousEntropyKNN(
                                samples: doubles,
                                k: knnDefaultK
                            ) {
                                return value
                            }
                            return Self.continuousEntropyKDE(
                                samples: doubles,
                                bandwidth: nil
                            )
                        case .kdeGaussian(let bandwidth):
                            return Self.continuousEntropyKDE(
                                samples: doubles,
                                bandwidth: bandwidth.map { Double($0) }
                            )
                        }
                    }
                },
                bootstrapSamples: bootstrapSamples,
                confidenceLevel: confidenceLevel,
                method: method
            )

            return Estimate(value: point, confidenceInterval: interval, replicates: bootstrapSamples, method: method)
        }

        /// Bootstrap estimate (point + confidence interval) for the KL divergence.
        ///
        /// - Parameters mirror ``entropyEstimate(estimator:bootstrapSamples:confidenceLevel:method:)`` but operate on two samples.
        ///   Resampling is performed independently for `self` and `other`.
        /// - Returns: Point estimate with optional confidence bounds, or `nil` when the divergence could not be computed.
        public func klDivergenceEstimate(
            relativeTo other: Self,
            estimator: DensityEstimator = .automatic,
            bootstrapSamples: Int = 200,
            confidenceLevel: T = 0.95,
            method: BootstrapMethod = .percentile
        ) throws -> Estimate<T>? {
            guard let point = try klDivergence(relativeTo: other, estimator: estimator) else {
                return nil
            }

            let interval = try Self.bootstrapTwoSample(
                dataP: samplesSorted,
                dataQ: other.samplesSorted,
                estimator: { sampleP, sampleQ -> T in
                    let tempP = try Distribution.Empirical(samples: sampleP)
                    let tempQ = try Distribution.Empirical(samples: sampleQ)
                    return try tempP.klDivergence(relativeTo: tempQ, estimator: estimator) ?? 0
                },
                bootstrapSamples: bootstrapSamples,
                confidenceLevel: confidenceLevel,
                method: method
            )

            return Estimate(value: point, confidenceInterval: interval, replicates: bootstrapSamples, method: method)
        }
    }
}

// MARK: - Helpers (private)

private extension Distribution.Empirical {
    static func normalisedFrequencies(
        counts: [Int],
        total: Int,
        alpha: T
    ) -> [T] {
        let unique = counts.count
        let totalT = T(total)
        let alphaSum = alpha * T(unique)
        let denominator = totalT + alphaSum
        return counts.map { count in
            (T(count) + alpha) / denominator
        }
    }

    static func multiplicities(for sample: [T]) -> [Int] {
        if sample.isEmpty { return [] }
        var counts: [Int] = []
        var iterator = sample.sorted().makeIterator()
        if let first = iterator.next() {
            var current = first
            var currentCount = 1
            while let next = iterator.next() {
                if next == current {
                    currentCount += 1
                } else {
                    counts.append(currentCount)
                    current = next
                    currentCount = 1
                }
            }
            counts.append(currentCount)
        }
        return counts
    }

    static func uniqueValues(from sample: [T]) -> [T] {
        if sample.isEmpty { return [] }
        var values: [T] = []
        var iterator = sample.sorted().makeIterator()
        if let first = iterator.next() {
            values.append(first)
            var current = first
            while let next = iterator.next() {
                if next != current {
                    values.append(next)
                    current = next
                }
            }
        }
        return values
    }

    static func detectDiscreteSupport(
        values: [T],
        totalCount: Int,
        toleranceScale: T
    ) -> (isDiscrete: Bool, latticeStep: T?, latticeOrigin: T?) {
        let uniqueCount = values.count
        if uniqueCount == 0 {
            return (true, nil, nil)
        }
        if uniqueCount == 1 {
            return (true, nil, values.first)
        }
        if uniqueCount == totalCount {
            return (false, nil, nil)
        }

        var diffs: [T] = []
        diffs.reserveCapacity(values.count - 1)
        for idx in 1..<values.count {
            let diff = values[idx] - values[idx - 1]
            if diff > 0 {
                diffs.append(diff)
            }
        }

        guard let minDiff = diffs.min(), minDiff > 0 else {
            return (true, nil, values.first)
        }

        let tolerance = max(toleranceScale * minDiff, T.leastNonzeroMagnitude * 4)
        var isLattice = true
        for diff in diffs {
            let ratio = diff / minDiff
            let rounded = ratio.rounded()
            let delta = abs(ratio - rounded)
            if delta > tolerance {
                isLattice = false
                break
            }
        }

        if isLattice {
            return (true, minDiff, values.first)
        }

        let duplicateFraction = 1 - T(values.count) / T(values.count + diffs.count)
        if duplicateFraction > 0.25 && values.count <= Int(Double(values.count * 2).squareRoot()) {
            return (true, nil, values.first)
        }

        return (false, nil, nil)
    }

    static func discreteEntropy(
        probabilities: [T],
        uniqueCount: Int,
        sampleSize: Int
    ) -> T {
        var entropy: T = 0
        for p in probabilities where p > 0 {
            entropy -= p * T.log(p)
        }
        if sampleSize > 0 {
            entropy += (T(uniqueCount) - 1) / (2 * T(sampleSize))
        }
        return entropy
    }

    static func continuousEntropyKNN(samples: [Double], k: Int) -> T? {
        guard samples.count > max(k, 1) else { return nil }
        let n = samples.count
        var sumLogDistances: Double = 0
        let sorted = samples.sorted()
        let epsilon = max(1e-12 * (sorted.last! - sorted.first!), Double.leastNonzeroMagnitude)
        for idx in 0..<n {
            let distance = Self.distanceToKthNeighbour(in: sorted, index: idx, k: k)
            sumLogDistances += log(max(distance * 2, epsilon))
        }
        let psiN = Double(try! SpecialFunctions.digamma(T(n)))
        let psiK = Double(try! SpecialFunctions.digamma(T(k)))
        let entropy = psiN - psiK + log(2.0) + sumLogDistances / Double(n)
        return T(entropy)
    }

    static func continuousEntropyKDE(samples: [Double], bandwidth: Double?) -> T {
        let bw = bandwidth ?? Self.selectKDEBandwidth(samples: samples)
        let n = samples.count
        let normalization = 1.0 / (Double(n - 1) * bw * sqrt(2 * Double.pi))
        var sum: Double = 0
        for i in 0..<n {
            var density = 0.0
            for j in 0..<n where j != i {
                let z = (samples[i] - samples[j]) / bw
                density += exp(-0.5 * z * z)
            }
            density *= normalization
            density = max(density, Double.leastNonzeroMagnitude)
            sum -= log(density)
        }
        return T(sum / Double(n))
    }

    static func selectKDEBandwidth(samples: [Double]) -> Double {
        let n = samples.count
        if n < 2 { return 1 }
        let mean = samples.reduce(0, +) / Double(n)
        let variance = samples.reduce(0) { $0 + pow($1 - mean, 2) } / Double(max(n - 1, 1))
        let stdev = sqrt(max(variance, Double.leastNonzeroMagnitude))
        let silverman = 1.06 * stdev * pow(Double(n), -0.2)
        let candidates = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0].map { $0 * silverman }
        return candidates.max { lhs, rhs in
            Self.kdeLogLikelihood(samples: samples, bandwidth: lhs) <
                Self.kdeLogLikelihood(samples: samples, bandwidth: rhs)
        } ?? silverman
    }

    static func kdeLogLikelihood(samples: [Double], bandwidth: Double) -> Double {
        let n = samples.count
        if n < 2 { return .leastNormalMagnitude }
        let normalization = 1.0 / (Double(n - 1) * bandwidth * sqrt(2 * Double.pi))
        var logLikelihood = 0.0
        for i in 0..<n {
            var density = 0.0
            for j in 0..<n where j != i {
                let z = (samples[i] - samples[j]) / bandwidth
                density += exp(-0.5 * z * z)
            }
            density = max(density * normalization, Double.leastNonzeroMagnitude)
            logLikelihood += log(density)
        }
        return logLikelihood
    }

    static func distanceToKthNeighbour(in sorted: [Double], index: Int, k: Int) -> Double {
        if sorted.count <= 1 { return .leastNonzeroMagnitude }
        var count = 0
        var left = index - 1
        var right = index + 1
        var lastDistance: Double = .leastNonzeroMagnitude
        while count < k {
            let leftDistance = left >= 0 ? abs(sorted[index] - sorted[left]) : Double.infinity
            let rightDistance = right < sorted.count ? abs(sorted[right] - sorted[index]) : Double.infinity
            if leftDistance < rightDistance {
                lastDistance = leftDistance
                left -= 1
            } else {
                lastDistance = rightDistance
                right += 1
            }
            count += 1
            if lastDistance == 0 && count >= k {
                break
            }
            if leftDistance == .infinity && rightDistance == .infinity {
                break
            }
        }
        return max(lastDistance, Double.leastNonzeroMagnitude)
    }

    static func upperBound(in array: [T], value: T) -> Int {
        var low = 0
        var high = array.count
        while low < high {
            let mid = (low + high) >> 1
            if array[mid] <= value {
                low = mid + 1
            } else {
                high = mid
            }
        }
        return low
    }

    static func binarySearch(values: [T], target: T) -> Int? {
        var low = 0
        var high = values.count - 1
        while low <= high {
            let mid = (low + high) >> 1
            let value = values[mid]
            if value == target {
                return mid
            } else if value < target {
                low = mid + 1
            } else {
                if mid == 0 { break }
                high = mid - 1
            }
        }
        return nil
    }

    static func quantile(
        samples: [T],
        probabilities: [T],
        p: T
    ) -> T {
        if probabilities.isEmpty { return samples.first! }
        var cumulative: T = 0
        for idx in probabilities.indices {
            cumulative += probabilities[idx]
            if cumulative >= p {
                return samples[idx]
            }
        }
        return samples.last!
    }

    static func smoothedProbability(value: T, distribution: Distribution.Empirical<T>, alpha: T) -> T {
        if let index = Self.binarySearch(values: distribution.uniqueValues, target: value) {
            let counts = distribution.counts
            let total = distribution.totalCount
            let unique = distribution.uniqueValues.count
            return (T(counts[index]) + alpha) / (T(total) + alpha * T(unique))
        }
        let unique = distribution.uniqueValues.count
        return alpha / (T(distribution.totalCount) + alpha * T(unique))
    }

    static func unionSupport(
        left: Distribution.Empirical<T>,
        right: Distribution.Empirical<T>
    ) -> [T] {
        var merged = left.uniqueValues
        for value in right.uniqueValues {
            if !merged.contains(where: { $0 == value }) {
                merged.append(value)
            }
        }
        merged.sort()
        return merged
    }

    static func klDivergenceKNN(samplesP: [Double], samplesQ: [Double], k: Int) -> Double? {
        let n = samplesP.count
        let m = samplesQ.count
        guard n > k, m >= k else { return nil }

        let sortedP = samplesP.sorted()
        let sortedQ = samplesQ.sorted()
        var sum: Double = 0
        let epsilon = max(1e-12 * (sortedP.last! - sortedP.first!), Double.leastNonzeroMagnitude)
        for i in 0..<n {
            let rho = distanceToKthNeighbour(in: sortedP, index: i, k: k)
            let nu = distanceToKthNeighbour(in: sortedQ, index: Self.closestIndex(in: sortedQ, to: sortedP[i]), k: k)
            sum += log(max(nu, epsilon) / max(rho, epsilon))
        }
        return sum / Double(n) + log(Double(m) / Double(n - 1))
    }

    static func closestIndex(in array: [Double], to value: Double) -> Int {
        var low = 0
        var high = array.count
        while low < high {
            let mid = (low + high) >> 1
            if array[mid] < value {
                low = mid + 1
            } else {
                high = mid
            }
        }
        if low == 0 { return 0 }
        if low >= array.count { return array.count - 1 }
        return abs(array[low] - value) < abs(value - array[low - 1]) ? low : low - 1
    }

    static func klDivergenceKDE(samplesP: [Double], samplesQ: [Double], bandwidth: Double?) -> Double? {
        let bwP = bandwidth ?? selectKDEBandwidth(samples: samplesP)
        let bwQ = bandwidth ?? selectKDEBandwidth(samples: samplesQ)
        let n = samplesP.count
        if n == 0 { return nil }
        var sum: Double = 0
        for value in samplesP {
            let p = kdeDensityGaussian(samples: samplesP, x: value, bandwidth: bwP, omitSelf: true)
            let q = kdeDensityGaussian(samples: samplesQ, x: value, bandwidth: bwQ, omitSelf: false)
            let densityP = max(p, Double.leastNonzeroMagnitude)
            let densityQ = max(q, Double.leastNonzeroMagnitude)
            sum += log(densityP / densityQ)
        }
        return sum / Double(n)
    }

    static func kdeDensityGaussian(
        samples: [Double],
        x: Double,
        bandwidth: Double?,
        omitSelf: Bool = false
    ) -> Double {
        let bw = bandwidth ?? selectKDEBandwidth(samples: samples)
        let n = samples.count
        if n == 0 { return 0 }
        let normalization = 1.0 / (Double(n - (omitSelf ? 1 : 0)) * bw * sqrt(2 * Double.pi))
        var density = 0.0
        for sample in samples {
            if omitSelf && sample == x { continue }
            let z = (x - sample) / bw
            density += exp(-0.5 * z * z)
        }
        density = max(density * normalization, Double.leastNonzeroMagnitude)
        return density
    }

    static func estimateContinuousMode(samples: [Double]) -> T? {
        let bandwidth = selectKDEBandwidth(samples: samples)
        var bestValue = samples.first ?? 0
        var bestDensity = -Double.infinity
        for point in samples {
            let density = kdeDensityGaussian(samples: samples, x: point, bandwidth: bandwidth)
            if density > bestDensity {
                bestDensity = density
                bestValue = point
            }
        }
        return T(bestValue)
    }

    static func detectMultimodality(samples: [Double], isDiscrete: Bool) -> Bool {
        if samples.count < 5 { return false }
        if isDiscrete {
            var localMaxima = 0
            let counts = Self.multiplicities(for: samples.map { T($0) })
            let probabilities = counts.map { Double($0) / Double(samples.count) }
            for idx in probabilities.indices where idx > 0 && idx < probabilities.endIndex - 1 {
                if probabilities[idx] >= probabilities[idx - 1] && probabilities[idx] > probabilities[idx + 1] {
                    localMaxima += 1
                }
            }
            return localMaxima > 1
        }
        guard let minSample = samples.min(), let maxSample = samples.max(), maxSample > minSample else {
            return false
        }
        let bandwidth = selectKDEBandwidth(samples: samples)
        let gridCount = min(128, max(16, samples.count * 6))
        let step = (maxSample - minSample) / Double(gridCount - 1)
        var peaks = 0
        var prevDensity = 0.0
        var prevPrevDensity = 0.0
        for index in 0..<gridCount {
            let x = minSample + Double(index) * step
            let density = kdeDensityGaussian(samples: samples, x: x, bandwidth: bandwidth)
            if index >= 2 {
                if prevDensity > prevPrevDensity && prevDensity > density {
                    peaks += 1
                }
            }
            prevPrevDensity = prevDensity
            prevDensity = density
        }
        return peaks > 1
    }

    static func bootstrapOneSample(
        data: [T],
        estimator: ([T]) throws -> T,
        bootstrapSamples: Int,
        confidenceLevel: T,
        method: BootstrapMethod
    ) throws -> (lower: T, upper: T)? {
        guard bootstrapSamples > 1 else { return nil }
        let original = try estimator(data)
        var replicates: [T] = []
        replicates.reserveCapacity(bootstrapSamples)
        var rng = SystemRandomNumberGenerator()
        for _ in 0..<bootstrapSamples {
            var resampled: [T] = []
            resampled.reserveCapacity(data.count)
            for _ in 0..<data.count {
                let index = Int.random(in: 0..<data.count, using: &rng)
                resampled.append(data[index])
            }
            replicates.append(try estimator(resampled))
        }
        return computeBootstrapInterval(
            replicates: replicates,
            original: original,
            confidenceLevel: confidenceLevel,
            method: method,
            data: data,
            estimator: estimator
        )
    }

    static func bootstrapTwoSample(
        dataP: [T],
        dataQ: [T],
        estimator: ([T], [T]) throws -> T,
        bootstrapSamples: Int,
        confidenceLevel: T,
        method: BootstrapMethod
    ) throws -> (lower: T, upper: T)? {
        guard bootstrapSamples > 1 else { return nil }
        let original = try estimator(dataP, dataQ)
        var replicates: [T] = []
        replicates.reserveCapacity(bootstrapSamples)
        var rng = SystemRandomNumberGenerator()
        for _ in 0..<bootstrapSamples {
            var resampledP: [T] = []
            var resampledQ: [T] = []
            resampledP.reserveCapacity(dataP.count)
            resampledQ.reserveCapacity(dataQ.count)
            for _ in 0..<dataP.count {
                let index = Int.random(in: 0..<dataP.count, using: &rng)
                resampledP.append(dataP[index])
            }
            for _ in 0..<dataQ.count {
                let index = Int.random(in: 0..<dataQ.count, using: &rng)
                resampledQ.append(dataQ[index])
            }
            replicates.append(try estimator(resampledP, resampledQ))
        }
        return computeBootstrapInterval(
            replicates: replicates,
            original: original,
            confidenceLevel: confidenceLevel,
            method: method,
            data: dataP,
            estimator: { _ in original } // Jackknife not available for two-sample; fall back to percentile
        )
    }

    static func computeBootstrapInterval(
        replicates: [T],
        original: T,
        confidenceLevel: T,
        method: BootstrapMethod,
        data: [T],
        estimator: ([T]) throws -> T
    ) -> (lower: T, upper: T)? {
        guard !replicates.isEmpty else { return nil }
        let alpha = (1 - confidenceLevel) / 2
        let sorted = replicates.sorted()
        switch method {
        case .percentile:
            let lower = quantile(sortedValues: sorted, probability: alpha)
            let upper = quantile(sortedValues: sorted, probability: 1 - alpha)
            return (lower, upper)
        case .bca:
            let lower = bcaQuantile(
                sortedReplicates: sorted,
                original: original,
                alpha: Double(alpha),
                data: data,
                estimator: estimator
            )
            let upper = bcaQuantile(
                sortedReplicates: sorted,
                original: original,
                alpha: Double(1 - alpha),
                data: data,
                estimator: estimator
            )
            return (lower, upper)
        }
    }

    static func quantile(sortedValues: [T], probability: T) -> T {
        let clampedP = min(max(probability, 0), 1)
        let position = clampedP * T(sortedValues.count - 1)
        let index = Int(position.rounded(.down))
        let nextIndex = min(sortedValues.count - 1, index + 1)
        let fraction = position - T(index)
        return sortedValues[index] * (1 - fraction) + sortedValues[nextIndex] * fraction
    }

    static func bcaQuantile(
        sortedReplicates: [T],
        original: T,
        alpha: Double,
        data: [T],
        estimator: ([T]) throws -> T
    ) -> T {
        let b = sortedReplicates.filter { $0 < original }.count
        let z0 = normalInverseCDF(Double(b) / Double(sortedReplicates.count))
        let jackknife = try? jackknifeReplicates(data: data, estimator: estimator)
        let acceleration = jackknife.map { accelerationCoefficient(from: $0) } ?? 0
        let zAlpha = normalInverseCDF(alpha)
        let adjusted = normalCDF(
            z0 + (z0 + zAlpha) / max(1 - acceleration * (z0 + zAlpha), Double.leastNonzeroMagnitude)
        )
        return quantile(sortedValues: sortedReplicates, probability: T(adjusted))
    }

    static func jackknifeReplicates(
        data: [T],
        estimator: ([T]) throws -> T
    ) throws -> [T] {
        var replicates: [T] = []
        replicates.reserveCapacity(data.count)
        for idx in data.indices {
            var subset = data
            subset.remove(at: idx)
            replicates.append(try estimator(subset))
        }
        return replicates
    }

    static func accelerationCoefficient(from jackknife: [T]) -> Double {
        guard !jackknife.isEmpty else { return 0 }
        let mean = jackknife.reduce(0) { $0 + Double($1) } / Double(jackknife.count)
        var numerator = 0.0
        var denominator = 0.0
        for value in jackknife {
            let diff = mean - Double(value)
            numerator += diff * diff * diff
            denominator += diff * diff
        }
        if denominator == 0 { return 0 }
        return numerator / (6 * pow(denominator, 1.5))
    }

    static func normalCDF(_ x: Double) -> Double {
        0.5 * (1 + erf(x / sqrt(2)))
    }

    static func normalInverseCDF(_ p: Double) -> Double {
        let clamped = min(max(p, Double.leastNonzeroMagnitude), 1 - Double.leastNonzeroMagnitude)
        if let inv = try? SpecialFunctions.inverseErrorFunction(T(2 * clamped - 1)) {
            return Double(inv) * sqrt(2)
        }
        // Fallback approximation (Acklam)
        let a: [Double] = [
            -3.969683028665376e+01,
            2.209460984245205e+02,
            -2.759285104469687e+02,
            1.383577518672690e+02,
            -3.066479806614716e+01,
            2.506628277459239e+00
        ]
        let b: [Double] = [
            -5.447609879822406e+01,
            1.615858368580409e+02,
            -1.556989798598866e+02,
            6.680131188771972e+01,
            -1.328068155288572e+01
        ]
        let c: [Double] = [
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838,
            -2.549732539343734,
            4.374664141464968,
            2.938163982698783
        ]
        let d: [Double] = [
            7.784695709041462e-03,
            3.224671290700398e-01,
            2.445134137142996,
            3.754408661907416
        ]
        let plow = 0.02425
        let phigh = 1 - plow
        if clamped < plow {
            let q = sqrt(-2 * log(clamped))
            return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        } else if clamped > phigh {
            let q = sqrt(-2 * log(1 - clamped))
            return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        } else {
            let q = clamped - 0.5
            let r = q * q
            return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
                (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
        }
    }
}
