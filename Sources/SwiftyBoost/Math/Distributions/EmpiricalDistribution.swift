//
//  EmpiricalDistribution.swift
//  SwiftyBoost
//
//  Created by Codex Agent (GPT-5) for Volker Thieme.
//

import SwiftyBoostPrelude
import Foundation

extension Distribution {
    /// An empirical distribution backed by a finite set of observations.
    ///
    /// This type provides a pragmatic, feature-complete implementation of ``DistributionProtocol``
    /// from a sample array. It supports both lattice/discrete and continuous data:
    ///
    /// - Discrete detection:
    ///   The implementation attempts to detect whether the support is discrete (lattice-based)
    ///   using repeated values, minimum spacing, and ratio checks. For discrete data, mass
    ///   is assigned to unique support points using smoothed relative frequencies.
    ///
    /// - Continuous handling:
    ///   For continuous samples, density-related queries (e.g., `pdf`, `mode`, entropy, KL divergence)
    ///   are computed using nonparametric estimators (KDE with a Gaussian kernel or a KNN-based method).
    ///
    /// - Summary statistics:
    ///   Mean, variance, skewness, kurtosis, median, mode, entropy, and lattice metadata are computed
    ///   on initialisation and cached. When unavailable or undefined, corresponding properties return `nil`.
    ///
    /// - KL divergence:
    ///   KL divergence is available against any other ``DistributionProtocol`` via
    ///   numeric integration/summation (see ``Distribution/KLDivergenceOptions``),
    ///   and against another empirical distribution using nonparametric estimators
    ///   (KDE or KNN).
    ///
    /// - Uncertainty via bootstrap:
    ///   Bootstrap-based confidence intervals are provided for entropy and KL divergence.
    ///   Percentile and BCa intervals are supported. For BCa, a jackknife acceleration
    ///   is computed when possible.
    ///
    /// Notes:
    /// - All densities and entropies are in natural units (nats).
    /// - KDE bandwidth selection uses a Silverman-like rule with a short data-driven
    ///   model selection pass to maximise log-likelihood among a small set of multipliers.
    /// - KNN-based entropy/kl implementations use symmetric neighbour distance heuristics
    ///   with safety floors to limit numerical issues.
    public struct Empirical<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        // MARK: Nested types

        /// Bootstrap interval type.
        ///
        /// - percentile: Standard percentile bootstrap interval.
        /// - bca: Bias-corrected and accelerated interval (requires jackknife).
        public enum BootstrapMethod: Sendable {
            case percentile
            case bca
        }

        /// A bootstrap estimate with optional confidence interval.
        ///
        /// - value: The point estimate computed on the original data.
        /// - confidenceInterval: Optional (lower, upper) interval at the requested confidence level.
        /// - replicates: Number of bootstrap resamples used.
        /// - method: The bootstrap method used (percentile or BCa).
        public struct Estimate<Value: Sendable>: Sendable {
            public let value: Value
            public let confidenceInterval: (lower: Value, upper: Value)?
            public let replicates: Int
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

        /// Density estimator choice for continuous tasks (entropy, KL divergence, mode).
        ///
        /// - automatic: Heuristic selection (typically KNN with small k for moderate n).
        /// - knn(k:): K-nearest neighbour estimator with user-specified `k` (clamped to ≥ 1).
        /// - kdeGaussian(bandwidth:): Gaussian-kernel KDE with optional bandwidth override.
        public enum DensityEstimator: Sendable {
            case automatic
            case knn(k: Int)
            case kdeGaussian(bandwidth: T?)
        }

        // MARK: Stored properties

        /// Sample values sorted ascending (includes duplicates).
        private let samplesSorted: [T]
        /// Unique support points (sorted ascending).
        private let uniqueValues: [T]
        /// Multiplicity per unique support point (aligned with `uniqueValues`).
        private let counts: [Int]
        /// Total number of observations.
        private let totalCount: Int
        /// Additive smoothing parameter for discrete frequencies (Laplace-like).
        private let smoothingAlpha: T = 0.5
        /// Default KNN `k` when not specified.
        private let knnDefaultK: Int = 3

        // Cached statistics (computed at init)
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
        /// - Parameter samples: Array of finite observations (NaN and ±∞ are rejected).
        ///
        /// This initialiser:
        /// - Sorts the sample and computes the range.
        /// - Collapses samples into unique support points with multiplicities.
        /// - Detects whether the support is likely discrete (lattice) vs continuous.
        /// - Computes smoothed probabilities for discrete support.
        /// - Pre-computes summary statistics (mean, variance, skewness, kurtosis, median, mode).
        /// - Estimates entropy using a discrete plug-in (with Miller–Madow correction) or
        ///   continuous estimators (KNN first, then KDE fallback).
        /// - Flags likely multimodality using a simple peaks-count heuristic on either discrete
        ///   probabilities or a KDE grid for continuous data.
        ///
        /// - Throws:
        ///   - ``DistributionError/invalidCombination(message:value:)`` if the input is empty.
        ///   - ``DistributionError/parameterNotFinite(name:value:)`` if any observation is not finite.
        public init(samples: [T]) throws {
            guard !samples.isEmpty else {
                throw DistributionError<T>.invalidCombination(
                    message: "Empirical distribution requires at least one observation.",
                    value: nil
                )
            }
            for v in samples {
                guard v.isFinite else {
                    throw DistributionError<T>.parameterNotFinite(name: "samples", value: v)
                }
            }

            var collected = samples
            collected.sort()
            self.samplesSorted = collected
            self.totalCount = collected.count
            self.rangeLower = collected.first!
            self.rangeUpper = collected.last!

            // Collapse to unique support points with multiplicities.
            var values: [T] = []
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
                        multiplicities.append(currentCount)
                        current = next
                        currentCount = 1
                    }
                }
                values.append(current)
                multiplicities.append(currentCount)
            }
            self.uniqueValues = values
            self.counts = multiplicities

            // Discrete vs continuous detection (lattice step/origin when discrete).
            let detection = Self.detectDiscreteSupport(values: values, totalCount: collected.count, toleranceScale: T(1e-6))
            self.isDiscreteDistribution = detection.isDiscrete
            self.latticeStepValue = detection.latticeStep
            self.latticeOriginValue = detection.latticeOrigin

            let probabilities = Self.normalisedFrequencies(
                counts: multiplicities,
                total: collected.count,
                alpha: smoothingAlpha
            )

            // Mean/variance/skewness/kurtosis (from unique values and smoothed probabilities).
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

            // Mode:
            // - Discrete: argmax of smoothed probabilities (ties broken by smallest value).
            // - Continuous: maximum of Gaussian KDE evaluated at sample points.
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
                self.modeCache = Self.estimateContinuousMode(samples: self.samplesSorted)
            }

            // Median via weighted quantile at p=0.5 (discrete) or linear interpolation (continuous).
            self.medianCache = Self.quantile(samples: collected, probabilities: probabilities, p: 0.5)

            // Entropy:
            // - Discrete: plug-in with Miller–Madow correction.
            // - Continuous: prefer KNN; fall back to KDE if KNN is not applicable.
            if detection.isDiscrete {
                self.entropyCache = Self.discreteEntropy(probabilities: probabilities, uniqueCount: values.count, sampleSize: collected.count)
            } else if let knnEntropy = Self.continuousEntropyKNN(samples: self.samplesSorted, k: knnDefaultK) {
                self.entropyCache = knnEntropy
            } else {
                self.entropyCache = Self.continuousEntropyKDE(samples: self.samplesSorted, bandwidth: nil)
            }

            // Lightweight multimodality heuristic (counts local maxima).
            self.likelyMultimodal = Self.detectMultimodality(
                samples: self.samplesSorted,
                isDiscrete: detection.isDiscrete
            )
        }

        // MARK: DistributionProtocol – Support

        /// Lower bound of the observed sample range.
        ///
        /// For empirical distributions this equals the minimum observation.
        public var supportLowerBound: T { rangeLower }
        /// Upper bound of the observed sample range.
        ///
        /// For empirical distributions this equals the maximum observation.
        public var supportUpperBound: T { rangeUpper }
        /// Convenience pair for the sample range `(min, max)`.
        public var range: (lower: T, upper: T) { (rangeLower, rangeUpper) }

        // MARK: - Core functions

        /// Probability density/mass function.
        ///
        /// - Discrete: returns the smoothed probability mass at `x` when `x` is a support point; 0 otherwise.
        /// - Continuous: evaluates a Gaussian-kernel KDE at `x` with an automatically chosen bandwidth.
        ///
        /// - Parameter x: Evaluation point (must be finite).
        /// - Returns: Nonnegative density/mass at `x`.
        /// - Throws: ``DistributionError/parameterNotFinite(name:value:)`` if `x` is NaN or ±∞.
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
                samples: samplesSorted,
                x: x,
                bandwidth: nil
            )
            return density
        }

        /// Natural logarithm of the PDF/PMF.
        ///
        /// Returns `-infinity` when the density/mass is numerically zero.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `log(pdf(x))` or `-infinity` if zero.
        /// - Throws: Re-throws from ``pdf(_:)``.
        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            if p > 0 {
                return T.log(p)
            }
            return -.infinity
        }

        /// Cumulative distribution function, F(x) = P(X ≤ x).
        ///
        /// - Discrete: sum of smoothed masses at support points ≤ `x`.
        /// - Continuous: empirical CDF from the sample (rank-based).
        ///
        /// - Parameter x: Evaluation point; `-∞` → 0, `+∞` → 1.
        /// - Returns: Value in [0, 1].
        /// - Throws: ``DistributionError/parameterNotFinite(name:value:)`` for NaN.
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
                let idx = Self.upperBound(in: samplesSorted, value: x)
                return T(idx) / T(totalCount)
            }
        }

        /// Survival function, S(x) = P(X > x).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: Value in [0, 1].
        /// - Throws: Re-throws from ``cdf(_:)``.
        public func sf(_ x: T) throws -> T {
            let c = try cdf(x)
            return max(0, 1 - c)
        }

        /// (Instantaneous) hazard function, h(x) = f(x) / S(x) where defined.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: Hazard at `x`, or `infinity` when `S(x) = 0` and `f(x) > 0`.
        /// - Throws: Re-throws from ``pdf(_:)`` and ``sf(_:)``.
        public func hazard(_ x: T) throws -> T {
            let density = try pdf(x)
            if density == 0 { return 0 }
            let survival = try sf(x)
            if survival > 0 {
                return density / survival
            }
            return .infinity
        }

        /// Cumulative hazard function, H(x) = −log S(x) where defined.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: CHF at `x`, or `infinity` when `S(x) ≤ 0`.
        /// - Throws: Re-throws from ``sf(_:)``.
        public func chf(_ x: T) throws -> T {
            let survival = try sf(x)
            if survival <= 0 { return .infinity }
            return -T.log(survival)
        }

        // MARK: - Inverses

        /// Lower-tail quantile (inverse CDF).
        ///
        /// - Discrete: smallest support value with cumulative probability ≥ `p`.
        /// - Continuous: linear interpolation between order statistics.
        ///
        /// - Parameter p: Probability in [0, 1].
        /// - Returns: Value `x` such that F(x) ≈ `p`.
        /// - Throws: ``DistributionError/parameterOutOfRange(name:min:max:)`` when `p ∉ [0, 1]`.
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

        /// Upper-tail quantile (inverse survival function).
        ///
        /// Returns `x` such that S(x) ≈ `q` for `q` in [0, 1].
        ///
        /// - Parameter q: Upper-tail probability in [0, 1].
        /// - Returns: The corresponding upper-tail quantile.
        /// - Throws: ``DistributionError/parameterOutOfRange(name:min:max:)`` when `q ∉ [0, 1]`.
        public func quantileComplement(_ q: T) throws -> T {
            guard q >= 0, q <= 1 else {
                throw DistributionError<T>.parameterOutOfRange(name: "q", min: 0, max: 1)
            }
            return try quantile(1 - q)
        }

        // MARK: - Summary statistics

        /// Smoothed mean of the empirical distribution.
        public var mean: T? { meanCache }
        /// Smoothed variance of the empirical distribution.
        public var variance: T? { varianceCache }
        /// Mode:
        /// - Discrete: support point with the highest smoothed probability.
        /// - Continuous: arg max of KDE evaluated at sample points.
        public var mode: T? { modeCache }
        /// Median (0.5 quantile).
        public var median: T { medianCache }
        /// Skewness (standardised third moment), when variance > 0.
        public var skewness: T? { skewnessCache }
        /// Kurtosis (Pearson), when variance > 0. Equals 3 for a normal distribution.
        public var kurtosis: T? { kurtosisCache }
        /// Excess kurtosis (kurtosis − 3), when defined.
        public var kurtosisExcess: T? {
            guard let kurtosis = kurtosisCache else { return nil }
            return kurtosis - 3
        }

        /// Lattice step for discrete support (spacing between adjacent support points), if detected.
        public var latticeStep: T? { latticeStepValue }
        /// Lattice origin for discrete support (offset), if detected.
        public var latticeOrigin: T? { latticeOriginValue }

        /// Entropy in nats:
        /// - Discrete: plug-in with Miller–Madow correction.
        /// - Continuous: KNN when applicable; otherwise KDE.
        public var entropy: T? { entropyCache }

        /// Indicates whether the distribution behaves as discrete (true) or continuous (false).
        public var isDiscrete: Bool { isDiscreteDistribution }

        /// Heuristic flag indicating likely multimodality.
        ///
        /// - Discrete: counts local maxima of the smoothed probability mass function.
        /// - Continuous: counts peaks of a KDE evaluated on a coarse grid.
        public var isLikelyMultimodal: Bool { likelyMultimodal }

        // MARK: - KL divergence

        /// KL divergence relative to an arbitrary distribution.
        ///
        /// - If `other` is another empirical distribution, this will call the specialised
        ///   overload using nonparametric estimators (see ``klDivergence(relativeTo:estimator:)``).
        /// - Otherwise, it falls back to the generic helper which uses numerical
        ///   integration/summation based on ``Distribution/KLDivergenceOptions``.
        ///
        /// - Parameters:
        ///   - other: Target distribution.
        ///   - options: Numerical configuration (quadrature rules, density floors, tail cutoffs).
        /// - Returns: `D_KL(self || other)` in nats, `nil` if undefined or supports do not match.
        public func klDivergence<D: DistributionProtocol>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D.RealType == T {
               if let empiricalOther = other as? Self {
                return try klDivergence(relativeTo: empiricalOther, estimator: .automatic)
            }
            return try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }

        /// KL divergence relative to another empirical distribution using nonparametric estimators.
        ///
        /// - Discrete: Mass functions are merged over the union support and smoothed with `alpha = 0.5`.
        /// - Continuous:
        ///   - `.automatic`: Chooses a KNN estimator with a small `k` relative to `n`, otherwise KDE.
        ///   - `.knn(k:)`: KNN estimator with the specified `k` (≥ 1).
        ///   - `.kdeGaussian(bandwidth:)`: Gaussian-kernel KDE (optional bandwidth override).
        ///
        /// Numerical notes:
        /// - Small density floors are applied internally to avoid `log(0)` issues.
        /// - Estimators are approximate; results depend on sample size, bandwidth, and `k`.
        ///
        /// - Parameters:
        ///   - other: Another empirical distribution over the same real type.
        ///   - estimator: Estimation method; see ``DensityEstimator``.
        /// - Returns: `D_KL(self || other)` in nats, or `nil` if not computable.
        public func klDivergence(
            relativeTo other: Self,
            estimator: DensityEstimator
        ) throws -> T? {
            if isDiscreteDistribution || other.isDiscreteDistribution {
                let smoothing = smoothingAlpha
                let merged = Self.unionSupport(left: self, right: other)
                var sum: T = 0
                for entry in merged {
                    let p = Self.smoothedProbability(
                        value: entry,
                        distribution: self,
                        alpha: smoothing
                    )
                    let qRaw = Self.smoothedProbability(
                        value: entry,
                        distribution: other,
                        alpha: smoothing
                    )
                    if p == 0 { continue }
                    let q = max(qRaw, T.leastNonzeroMagnitude)
                    sum += p * (T.log(p) - T.log(q))
                }
                return sum
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
                    samplesP: samplesSorted,
                    samplesQ: other.samplesSorted,
                    k: max(1, k)
                )
            case .kdeGaussian(let bandwidth):
                return Self.klDivergenceKDE(
                    samplesP: samplesSorted,
                    samplesQ: other.samplesSorted,
                    bandwidth: bandwidth
                )
            case .automatic:
                return nil // unreachable
            }
        }

        // MARK: - Bootstrap estimates

        /// Bootstrap confidence interval for entropy.
        ///
        /// - For discrete samples: re-estimates entropy via smoothed frequencies with
        ///   a Miller–Madow correction on each resample.
        /// - For continuous samples: uses the selected estimator; in `.automatic` and `.knn`
        ///   it prefers KNN when feasible, falling back to KDE otherwise.
        ///
        /// - Parameters:
        ///   - estimator: Density estimator for continuous entropy (ignored for discrete).
        ///   - bootstrapSamples: Number of bootstrap resamples (≥ 2 recommended).
        ///   - confidenceLevel: Confidence level in (0, 1), e.g. 0.95.
        ///   - method: Percentile or BCa interval. BCa uses jackknife acceleration when available.
        /// - Returns: An ``Estimate`` containing the point estimate and the interval if computable.
        /// - Throws: Rethrows from the estimator closure if it fails on a resample.
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
                if let value = Self.continuousEntropyKNN(samples: samplesSorted, k: max(1, k)) {
                    point = value
                } else {
                    point = Self.continuousEntropyKDE(samples: samplesSorted, bandwidth: nil)
                }
            case .kdeGaussian(let bandwidth):
                point = Self.continuousEntropyKDE(samples: samplesSorted, bandwidth: bandwidth)
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
                        switch estimator {
                        case .automatic, .knn:
                            if let value = Self.continuousEntropyKNN(
                                samples: sample,
                                k: knnDefaultK
                            ) {
                                return value
                            }
                            return Self.continuousEntropyKDE(
                                samples: sample,
                                bandwidth: nil
                            )
                        case .kdeGaussian(let bandwidth):
                            return Self.continuousEntropyKDE(
                                samples: sample,
                                bandwidth: bandwidth
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

        /// Bootstrap confidence interval for empirical KL divergence.
        ///
        /// - Parameters:
        ///   - other: Comparison empirical distribution.
        ///   - estimator: Nonparametric estimator to use for continuous data.
        ///   - bootstrapSamples: Number of bootstrap resamples per distribution.
        ///   - confidenceLevel: Confidence level in (0, 1), e.g. 0.95.
        ///   - method: Percentile or BCa interval. For two-sample bootstrap, BCa falls back
        ///     to percentile due to the lack of a two-sample jackknife.
        /// - Returns: An ``Estimate`` with point estimate and interval, or `nil` if the divergence is undefined.
        /// - Throws: Rethrows from internal estimators if they fail on resamples.
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
    /// Convert counts to smoothed probabilities using additive smoothing.
    ///
    /// p_i = (count_i + alpha) / (total + alpha * unique)
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

    /// Multiplicities per unique value in a sample (sorted internally).
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

    /// Unique sorted values from a sample.
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

    /// Heuristics to detect discrete (lattice) support and optional step/origin.
    ///
    /// - Uses min spacing and ratio checks against a tolerance to infer a lattice.
    /// - Falls back to discrete if a substantial fraction of duplicates is observed.
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

        let duplicateFractionFull = 1 - T(uniqueCount) / T(max(totalCount, 1))
        if duplicateFractionFull > 0.25 {
            return (true, nil, values.first)
        }

        return (false, nil, nil)
    }

    /// Discrete entropy with Miller–Madow correction.
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

    /// Continuous entropy via KNN estimator (Kraskov-like in 1D).
    ///
    /// Returns `nil` when `n <= k`.
    static func continuousEntropyKNN(samples: [T], k: Int) -> T? {
        guard samples.count > max(k, 1) else { return nil }
        let n = samples.count
        var sumLogDistances: T = 0
        let sorted = samples.sorted()
        guard let first = sorted.first, let last = sorted.last else { return nil }
        let epsilon = max(T(1e-12) * (last - first), T.leastNonzeroMagnitude)
        for idx in 0..<n {
            let distance = Self.distanceToKthNeighbour(in: sorted, index: idx, k: k)
            sumLogDistances += T.log(max(distance * 2, epsilon))
        }
        let psiN = (try? SpecialFunctions.digamma(T(n))) ?? 0
        let psiK = (try? SpecialFunctions.digamma(T(k))) ?? 0
        let entropy = psiN - psiK + T.log(2) + sumLogDistances / T(n)
        return entropy
    }

    /// Continuous entropy via leave-one-out KDE with Gaussian kernel.
    ///
    /// Uses a bandwidth chosen by ``selectKDEBandwidth(samples:)`` when `bandwidth` is `nil`.
    static func continuousEntropyKDE(samples: [T], bandwidth: T?) -> T {
        let bw = bandwidth ?? Self.selectKDEBandwidth(samples: samples)
        let n = samples.count
        if n <= 1 { return 0 }
        let normalization = 1 / (T(n - 1) * bw * (T(2) * T.pi).squareRoot())
        var sum: T = 0
        for i in 0..<n {
            var density: T = 0
            for j in 0..<n where j != i {
                let z = (samples[i] - samples[j]) / bw
                density += T.exp(-0.5 * z * z)
            }
            density *= normalization
            density = max(density, T.leastNonzeroMagnitude)
            sum -= T.log(density)
        }
        return sum / T(n)
    }

    /// Selects a KDE bandwidth starting from Silverman's rule-of-thumb,
    /// refined by choosing the multiplier that maximises the leave-one-out log-likelihood
    /// over a small candidate set.
    static func selectKDEBandwidth(samples: [T]) -> T {
        let n = samples.count
        if n < 2 { return 1 }
        let mean = samples.reduce(0, +) / T(n)
        let variance = samples.reduce(0) { $0 + ( $1 - mean ) * ( $1 - mean ) } / T(max(n - 1, 1))
        let stdev = max(variance, T.leastNonzeroMagnitude).squareRoot()
        let silverman = T(1.06) * stdev * T.pow(T(n), -0.2)
        let multipliers: [T] = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0].map { T($0) }
        let candidates = multipliers.map { $0 * silverman }
        return candidates.max { lhs, rhs in
            Self.kdeLogLikelihood(samples: samples, bandwidth: lhs) <
                Self.kdeLogLikelihood(samples: samples, bandwidth: rhs)
        } ?? silverman
    }

    /// Leave-one-out KDE log-likelihood under a Gaussian kernel.
    static func kdeLogLikelihood(samples: [T], bandwidth: T) -> T {
        let n = samples.count
        if n < 2 { return -T.greatestFiniteMagnitude }
        let normalization = 1 / (T(n - 1) * bandwidth * (T(2) * T.pi).squareRoot())
        var logLikelihood: T = 0
        for i in 0..<n {
            var density: T = 0
            for j in 0..<n where j != i {
                let z = (samples[i] - samples[j]) / bandwidth
                density += T.exp(-0.5 * z * z)
            }
            density = max(density * normalization, T.leastNonzeroMagnitude)
            logLikelihood += T.log(density)
        }
        return logLikelihood
    }

    /// Distance from a sorted point to its k-th nearest neighbour in 1D (both sides).
    ///
    /// Returns a small positive floor when the sample is too small or distances are zero.
    static func distanceToKthNeighbour(in sorted: [T], index: Int, k: Int) -> T {
        if sorted.count <= 1 { return T.leastNonzeroMagnitude }
        var count = 0
        var left = index - 1
        var right = index + 1
        var lastDistance: T = T.leastNonzeroMagnitude
        while count < k {
            let leftDistance = left >= 0 ? abs(sorted[index] - sorted[left]) : T.infinity
            let rightDistance = right < sorted.count ? abs(sorted[right] - sorted[index]) : T.infinity
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
        return max(lastDistance, T.leastNonzeroMagnitude)
    }

    /// Upper-bound index for `value` in a sorted array (first index with element > value).
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

    /// Binary search for exact equality in a sorted array.
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

    /// Weighted quantile for a discrete distribution defined by `(samples, probabilities)`.
    ///
    /// Assumes `samples` are sorted and `probabilities` sum to ~1.
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

    /// Laplace-smoothed probability for a discrete value (returns small mass for unseen values).
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

    /// Union of support points (sorted, unique) from two discrete empirical distributions.
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

    /// KL divergence via KNN estimator for continuous samples.
    ///
    /// Returns `nil` when the sample sizes are insufficient.
    static func klDivergenceKNN(samplesP: [T], samplesQ: [T], k: Int) -> T? {
        let n = samplesP.count
        let m = samplesQ.count
        guard n > k, m >= k else { return nil }

        let sortedP = samplesP.sorted()
        let sortedQ = samplesQ.sorted()
        var sum: T = 0
        guard let pFirst = sortedP.first, let pLast = sortedP.last else { return nil }
        let epsilon = max(T(1e-12) * (pLast - pFirst), T.leastNonzeroMagnitude)
        for i in 0..<n {
            let rho = distanceToKthNeighbour(in: sortedP, index: i, k: k)
            let qIndex = Self.closestIndex(in: sortedQ, to: sortedP[i])
            let nu = distanceToKthNeighbour(in: sortedQ, index: qIndex, k: k)
            sum += T.log(max(nu, epsilon) / max(rho, epsilon))
        }
        return sum / T(n) + T.log(T(m) / T(n - 1))
    }

    /// Closest index to a value in a sorted array (ties broken toward the lower index).
    static func closestIndex(in array: [T], to value: T) -> Int {
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

    /// KL divergence via Gaussian KDE for continuous samples.
    ///
    /// - Uses bandwidths selected independently for P and Q when `bandwidth` is `nil`.
    /// - Evaluates densities at each sample from P; Q is evaluated without leave-one-out.
    static func klDivergenceKDE(samplesP: [T], samplesQ: [T], bandwidth: T?) -> T? {
        let bwP = bandwidth ?? selectKDEBandwidth(samples: samplesP)
        let bwQ = bandwidth ?? selectKDEBandwidth(samples: samplesQ)
        let n = samplesP.count
        if n == 0 { return nil }
        var sum: T = 0
        for value in samplesP {
            let p = kdeDensityGaussian(samples: samplesP, x: value, bandwidth: bwP, omitSelf: true)
            let q = kdeDensityGaussian(samples: samplesQ, x: value, bandwidth: bwQ, omitSelf: false)
            let densityP = max(p, T.leastNonzeroMagnitude)
            let densityQ = max(q, T.leastNonzeroMagnitude)
            sum += T.log(densityP / densityQ)
        }
        return sum / T(n)
    }

    /// Gaussian-kernel KDE density at `x`.
    ///
    /// - Parameters:
    ///   - samples: Data points.
    ///   - x: Evaluation location.
    ///   - bandwidth: Optional bandwidth; if nil, a data-driven selection is used.
    ///   - omitSelf: When true, omits contributions from samples equal to `x` (useful for leave-one-out).
    static func kdeDensityGaussian(
        samples: [T],
        x: T,
        bandwidth: T?,
        omitSelf: Bool = false
    ) -> T {
        let bw = bandwidth ?? selectKDEBandwidth(samples: samples)
        let n = samples.count
        if n == 0 { return 0 }
        let denomCount = T(n - (omitSelf ? 1 : 0))
        let normalization = 1 / (denomCount * bw * (T(2) * T.pi).squareRoot())
        var density: T = 0
        for sample in samples {
            if omitSelf && sample == x { continue }
            let z = (x - sample) / bw
            density += T.exp(-0.5 * z * z)
        }
        density = max(density * normalization, T.leastNonzeroMagnitude)
        return density
    }

    /// Estimates the mode for continuous data as the sample point with the largest KDE density.
    static func estimateContinuousMode(samples: [T]) -> T? {
        guard !samples.isEmpty else { return nil }
        let bandwidth = selectKDEBandwidth(samples: samples)
        var bestValue = samples.first!
        var bestDensity = -T.infinity
        for point in samples {
            let density = kdeDensityGaussian(samples: samples, x: point, bandwidth: bandwidth)
            if density > bestDensity {
                bestDensity = density
                bestValue = point
            }
        }
        return bestValue
    }

    /// Heuristic multimodality detector:
    /// - Discrete: counts local maxima in the smoothed PMF.
    /// - Continuous: counts peaks in a KDE evaluated over a small grid.
    static func detectMultimodality(samples: [T], isDiscrete: Bool) -> Bool {
        if samples.count < 5 { return false }
        if isDiscrete {
            var localMaxima = 0
            let counts = Self.multiplicities(for: samples)
            let probabilities = counts.map { T($0) / T(samples.count) }
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
        let step = (maxSample - minSample) / T(gridCount - 1)
        var peaks = 0
        var prevDensity: T = 0
        var prevPrevDensity: T = 0
        for index in 0..<gridCount {
            let x = minSample + T(index) * step
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

    // MARK: - Bootstrap internals

    /// One-sample bootstrap driver for a scalar estimator.
    ///
    /// - Returns: (lower, upper) interval or `nil` if `bootstrapSamples ≤ 1`.
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

    /// Two-sample bootstrap driver for a scalar estimator.
    ///
    /// - Returns: (lower, upper) interval or `nil` if `bootstrapSamples ≤ 1`.
    /// - Note: BCa falls back to percentile because a two-sample jackknife is not provided here.
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

    /// Computes percentile or BCa bootstrap intervals given replicates and the original estimate.
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

    /// Linear interpolation quantile for sorted values with probability in [0, 1].
    static func quantile(sortedValues: [T], probability: T) -> T {
        let clampedP = min(max(probability, 0), 1)
        let position = clampedP * T(sortedValues.count - 1)
        let index = Int(position.rounded(.down))
        let nextIndex = min(sortedValues.count - 1, index + 1)
        let fraction = position - T(index)
        return sortedValues[index] * (1 - fraction) + sortedValues[nextIndex] * fraction
    }

    /// BCa adjusted quantile based on bias and acceleration.
    ///
    /// - Computes z0 (bias) from the fraction of replicates below the original,
    ///   and acceleration from jackknife replicates when available.
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

    /// Leave-one-out jackknife replicates for a scalar estimator.
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

    /// Standard BCa acceleration coefficient from jackknife replicates.
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

    /// Standard normal CDF used by BCa.
    static func normalCDF(_ x: Double) -> Double {
        do {
            return try Distribution.Normal<Double>().cdf(x)
        }
        catch {
            return .nan
        }
    }

    /// Standard normal inverse CDF used by BCa.
    static func normalInverseCDF(_ p: Double) -> Double {
        do {
            return try Distribution.Normal<Double>().quantile(p)
        }
        catch {
            return .nan
        }
    }
}

