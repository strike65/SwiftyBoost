//
//  DynamicDistribution.swift
//  Unified runtime distribution over CBoostBridge generic vtable
//
//  Example:
//  let d = try Distribution.Dynamic<Double>(distributionName: "gamma", parameters: ["shape": 4.5, "scale": 1.0])
//  let p = try d.pdf(2.0)
//

import SwiftyBoostPrelude

extension Distribution {
    /// Internal runtime-backed distribution wrapper that binds Swift generic entry points
    /// to the Boost.Math vtable exposed by `CBoostBridge`.
    ///
    /// Each instance owns the underlying C handle and forwards calls through function
    /// pointers supplied by the factory. Public distribution types (e.g. ``Distribution/Gamma``)
    /// delegate to this implementation for consistency.
    internal struct Dynamic<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Type-erased storage that holds whichever Boost-backed handle matches `T`.
        private final class Box: @unchecked Sendable {
            /// Double-precision vtable, present when `T == Double` or used as a fallback.
            var d: bs_dist_d?
            /// Single-precision vtable, present when `T == Float`.
            var f: bs_dist_f?
            /// Extended-precision vtable, present when `T == Float80` on supported architectures.
            var l: bs_dist_l?

            init(d: bs_dist_d) { self.d = d }
            init(f: bs_dist_f) { self.f = f }
            init(l: bs_dist_l) { self.l = l }
            deinit {
                if let h = d { h.free?(h.ctx) }
                if let h = f { h.free?(h.ctx) }
                if let h = l { h.free?(h.ctx) }
            }
        }

        /// Mutable flag used to propagate divergence state across integration callbacks.
        private final class DivergenceFlag: @unchecked Sendable {
            /// Indicates whether the KL divergence evaluation exceeded finite bounds.
            var isInfinite: Bool = false
        }

        /// Precision-specific distribution handle retained for the lifetime of `self`.
        private let box: Box
        /// Normalized distribution identifier (lowercase, underscores) used for runtime fallbacks.
        private let nameNormalized: String
        /// Lowercased parameter dictionary for convenience lookups.
        private let paramsLower: [String: T]

        public init(distributionName name: String, parameters: [String: T]) throws {
            // Build C parameter array with strdup'd keys (freed immediately after factory returns)
            var cStrings: [UnsafeMutablePointer<CChar>] = []
            defer {
                cStrings.forEach { free($0) }
            }

            // Normalize name and parameter keys for Swift-side fallbacks
            self.nameNormalized = name
                .trimmingCharacters(in: .whitespacesAndNewlines)
                .lowercased()
                .replacingOccurrences(of: "-", with: "_")

            var pl: [String: T] = [:]
            pl.reserveCapacity(parameters.count)
            for (k, v) in parameters {
                pl[k.lowercased()] = v
            }
            self.paramsLower = pl

            if T.self == Double.self {
                var cParams: [bs_param_d] = []
                cParams.reserveCapacity(parameters.count)
                for (k, v) in parameters {
                    guard let ks = strdup(k) else { continue }
                    cStrings.append(ks)
                    cParams.append(
                        bs_param_d(
                            key: UnsafePointer(ks),
                            value: Double(v)
                        )
                    )
                }
                let ptr = UnsafeMutablePointer<bs_dist_d>.allocate(capacity: 1)
                defer { ptr.deallocate() }
                let ok = bs_dist_make_d(name, cParams, cParams.count, ptr)
                guard ok else {
                    throw DistributionError<Double>.invalidCombination(message: "Unknown distribution or bad parameters: \(name)", value: nil)
                }
                self.box = Box(d: ptr.pointee)
            } else if T.self == Float.self {
                var cParams: [bs_param_f] = []
                cParams.reserveCapacity(parameters.count)
                for (k, v) in parameters {
                    guard let ks = strdup(k) else { continue }
                    cStrings.append(ks)
                    cParams.append(
                        bs_param_f(
                            key: UnsafePointer(ks),
                            value: Float(v)
                        )
                    )
                }
                let ptr = UnsafeMutablePointer<bs_dist_f>.allocate(capacity: 1)
                defer { ptr.deallocate() }
                let ok = bs_dist_make_f(name, cParams, cParams.count, ptr)
                guard ok else {
                    throw DistributionError<Double>.invalidCombination(message: "Unknown distribution or bad parameters: \(name)", value: nil)
                }
                self.box = Box(f: ptr.pointee)
            } else {
                #if arch(x86_64) || arch(i386)
                var cParams: [bs_param_l] = []
                cParams.reserveCapacity(parameters.count)
                for (k, v) in parameters {
                    guard let ks = strdup(k) else { continue }
                    cStrings.append(ks)
                    cParams.append(
                        bs_param_l(
                            key: UnsafePointer(ks),
                            value: Float80(v)
                        )
                    )
                }
                let ptr = UnsafeMutablePointer<bs_dist_l>.allocate(capacity: 1)
                defer { ptr.deallocate() }
                let ok = bs_dist_make_l(name, cParams, cParams.count, ptr)
                guard ok else {
                    throw DistributionError<Double>.invalidCombination(message: "Unknown distribution or bad parameters: \(name)", value: nil)
                }
                self.box = Box(l: ptr.pointee)
                #else
                // Fallback to double factory
                var cParams: [bs_param_d] = []
                cParams.reserveCapacity(parameters.count)
                for (k, v) in parameters {
                    guard let ks = strdup(k) else { continue }
                    cStrings.append(ks)
                    cParams.append(
                        bs_param_d(
                            key: UnsafePointer(ks),
                            value: Double(v)
                        )
                    )
                }
                let ptr = UnsafeMutablePointer<bs_dist_d>.allocate(capacity: 1)
                defer { ptr.deallocate() }
                let ok = bs_dist_make_d(name, cParams, cParams.count, ptr)
                guard ok else {
                    throw DistributionError<Double>.invalidCombination(message: "Unknown distribution or bad parameters: \(name)", value: nil)
                }
                self.box = Box(d: ptr.pointee)
                #endif
            }
        }

        // MARK: - DistributionProtocol conformance

        /// Lower bound of the distribution support as reported by the underlying Boost handle.
        public var supportLowerBound: T {
            if let h = box.d {
                let r = h.range?(h.ctx) ?? bs_range_d(lower: .nan, upper: .nan)
                return T(r.lower)
            }
            if let h = box.f {
                let r = h.range?(h.ctx) ?? bs_range_f(lower: .nan, upper: .nan)
                return T(r.lower)
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l {
                let r = h.range?(h.ctx) ?? bs_range_l(lower: .nan, upper: .nan)
                return T(r.lower)
            }
            #endif
            return .nan
        }

        /// Upper bound of the distribution support as reported by the underlying Boost handle.
        public var supportUpperBound: T {
            if let h = box.d {
                let r = h.range?(h.ctx) ?? bs_range_d(lower: .nan, upper: .nan)
                return T(r.upper)
            }
            if let h = box.f {
                let r = h.range?(h.ctx) ?? bs_range_f(lower: .nan, upper: .nan)
                return T(r.upper)
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l {
                let r = h.range?(h.ctx) ?? bs_range_l(lower: .nan, upper: .nan)
                return T(r.upper)
            }
            #endif
            return .nan
        }

        /// Convenience tuple containing both support endpoints.
        public var range: (lower: T, upper: T) {
            (supportLowerBound, supportUpperBound)
        }

        /// Returns the density/mass evaluated at `x` using the active vtable entry.
        public func pdf(_ x: T) throws -> T {
            if let h = box.d, let f = h.pdf {
                return T(f(h.ctx, Double(x)))
            }
            if let h = box.f, let f = h.pdf {
                return T(f(h.ctx, Float(x)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.pdf {
                return T(f(h.ctx, Float80(Double(x))))
            }
            #endif
            return .nan
        }

        /// Returns the natural logarithm of the density/mass at `x`.
        public func logPdf(_ x: T) throws -> T {
            if let h = box.d, let f = h.logpdf {
                return T(f(h.ctx, Double(x)))
            }
            if let h = box.f, let f = h.logpdf {
                return T(f(h.ctx, Float(x)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.logpdf {
                return T(f(h.ctx, Float80(Double(x))))
            }
            #endif
            return .nan
        }

        /// Returns the lower-tail cumulative probability at `x`.
        public func cdf(_ x: T) throws -> T {
            if let h = box.d, let f = h.cdf {
                return T(f(h.ctx, Double(x)))
            }
            if let h = box.f, let f = h.cdf {
                return T(f(h.ctx, Float(x)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.cdf {
                return T(f(h.ctx, Float80(Double(x))))
            }
            #endif
            return .nan
        }

        /// Returns the survival (upper-tail) probability at `x`.
        public func sf(_ x: T) throws -> T {
            if let h = box.d, let f = h.sf {
                return T(f(h.ctx, Double(x)))
            }
            if let h = box.f, let f = h.sf {
                return T(f(h.ctx, Float(x)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.sf {
                return T(f(h.ctx, Float80(Double(x))))
            }
            #endif
            return .nan
        }

        /// Returns the hazard rate `f(x) / S(x)` when provided by the backend.
        public func hazard(_ x: T) throws -> T {
            if let h = box.d, let f = h.hazard {
                return T(f(h.ctx, Double(x)))
            }
            if let h = box.f, let f = h.hazard {
                return T(f(h.ctx, Float(x)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.hazard {
                return T(f(h.ctx, Float80(Double(x))))
            }
            #endif
            return .nan
        }

        /// Returns the cumulative hazard `-log S(x)` when provided by the backend.
        public func chf(_ x: T) throws -> T {
            if let h = box.d, let f = h.chf {
                return T(f(h.ctx, Double(x)))
            }
            if let h = box.f, let f = h.chf {
                return T(f(h.ctx, Float(x)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.chf {
                return T(f(h.ctx, Float80(Double(x))))
            }
            #endif
            return .nan
        }

        /// Returns the lower-tail quantile associated with probability `p`.
        public func quantile(_ p: T) throws -> T {
            if let h = box.d, let f = h.quantile {
                return T(f(h.ctx, Double(p)))
            }
            if let h = box.f, let f = h.quantile {
                return T(f(h.ctx, Float(p)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.quantile {
                return T(f(h.ctx, Float80(Double(p))))
            }
            #endif
            return .nan
        }

        /// Returns the upper-tail quantile associated with survival probability `q`.
        public func quantileComplement(_ q: T) throws -> T {
            if let h = box.d, let f = h.quantile_complement {
                return T(f(h.ctx, Double(q)))
            }
            if let h = box.f, let f = h.quantile_complement {
                return T(f(h.ctx, Float(q)))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.quantile_complement {
                return T(f(h.ctx, Float80(Double(q))))
            }
            #endif
            return .nan
        }

        /// Finite mean when available; otherwise `nil`.
        public var mean: T? {
            if let h = box.d, let f = h.mean {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.mean {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.mean {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            return nil
        }

        /// Representative mode of the distribution when provided.
        public var mode: T? {
            if let h = box.d, let f = h.mode {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.mode {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.mode {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            return nil
        }

        /// Finite variance when available; otherwise `nil`.
        public var variance: T? {
            if let h = box.d, let f = h.variance {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.variance {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.variance {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            return nil
        }

        /// Median of the distribution, falling back to the 0.5 quantile if necessary.
        public var median: T {
            if let h = box.d, let f = h.median {
                return T(f(h.ctx))
            }
            if let h = box.f, let f = h.median {
                return T(f(h.ctx))
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.median {
                return T(f(h.ctx))
            }
            #endif
            return .nan
        }

        /// Finite skewness measure when supplied by the backend.
        public var skewness: T? {
            if let h = box.d, let f = h.skewness {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.skewness {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.skewness {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            return nil
        }

        /// Pearson kurtosis when supplied by the backend.
        public var kurtosis: T? {
            if let h = box.d, let f = h.kurtosis {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.kurtosis {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.kurtosis {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            return nil
        }

        /// Excess kurtosis (`kurtosis - 3`) when supplied by the backend.
        public var kurtosisExcess: T? {
            if let h = box.d, let f = h.kurtosis_excess {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.kurtosis_excess {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.kurtosis_excess {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            return nil
        }

        /// Lattice spacing extracted from the discrete backend, when applicable.
        public var latticeStep: T? { discreteLattice?.step }

        /// Lattice origin extracted from the discrete backend, when applicable.
        public var latticeOrigin: T? { discreteLattice?.origin }

        /// Shannon or differential entropy in nats, when available.
        public var entropy: T? {
            if let h = box.d, let f = h.entropy {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            if let h = box.f, let f = h.entropy {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #if arch(x86_64) || arch(i386)
            if let h = box.l, let f = h.entropy {
                let v = f(h.ctx)
                return v.isFinite ? T(v) : nil
            }
            #endif
            // Swift-side fallback for distributions whose entropy is not provided by the vtable
            return fallbackEntropy()
        }
        
        public func klDivergence(
            relativeTo other: Distribution.Dynamic<T>,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            guard isDiscrete == other.isDiscrete else { return nil }
            if isDiscrete {
                return try discreteKLDivergence(relativeTo: other, options: options)
            } else {
                return try continuousKLDivergence(relativeTo: other, options: options)
            }
        }
        
        /// Boolean indicating whether the distribution operates on a discrete lattice.
        public var isDiscrete: Bool { discreteLattice != nil }

        // MARK: - KL divergence helpers

        private func continuousKLDivergence(
            relativeTo other: Distribution.Dynamic<T>,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            if let analytic = analyticKLDivergence(relativeTo: other) {
                return analytic
            }
            guard let bounds = sharedBounds(with: other) else { return nil }
            if !(bounds.lower < bounds.upper) {
                return nil
            }
            let densityFloor = max(options.densityFloor, T.leastNonzeroMagnitude)
            let flag = DivergenceFlag()
            let integrand: @Sendable (T) -> T = { x in
                let p = self.positiveDensity(self, at: x)
                if !p.isFinite || p <= densityFloor { return .zero }
                var q = self.positiveDensity(other, at: x)
                if !q.isFinite || q <= .zero {
                    flag.isInfinite = true
                    return .zero
                }
                if q < densityFloor {
                    q = densityFloor
                }
                let ratio = p / q
                if !ratio.isFinite || ratio <= .zero {
                    return .zero
                }
                return p * T.log(ratio)
            }

            let value: T
            switch (bounds.lower.isFinite, bounds.upper.isFinite) {
            case (true, true):
                let lowerD = Double(bounds.lower)
                let upperD = Double(bounds.upper)
                guard lowerD.isFinite, upperD.isFinite, lowerD < upperD else { return nil }
                let integrator = try Quadrature.Integrator<T>(rule: options.finiteRule)
                value = try integrator
                    .integrate(over: .finite(lower: lowerD, upper: upperD), integrand: integrand)
                    .value
            case (true, false):
                let lower = bounds.lower
                guard lower.isFinite else { return nil }
                let integrator = try Quadrature.Integrator<T>(rule: options.semiInfiniteRule)
                value = try integrator.integrate { shift in
                    integrand(lower + shift)
                }.value
            case (false, true):
                let upper = bounds.upper
                guard upper.isFinite else { return nil }
                let integrator = try Quadrature.Integrator<T>(rule: options.semiInfiniteRule)
                value = try integrator.integrate { shift in
                    integrand(upper - shift)
                }.value
            default:
                let integrator = try Quadrature.Integrator<T>(rule: options.infiniteRule)
                value = try integrator.integrate(over: .automatic, integrand: integrand).value
            }

            if flag.isInfinite {
                return T.infinity
            }
            return value.isFinite ? value : nil
        }

        private func discreteKLDivergence(
            relativeTo other: Distribution.Dynamic<T>,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            guard let bounds = sharedBounds(with: other) else { return nil }
            let densityFloor = max(options.densityFloor, T.leastNonzeroMagnitude)
            let tailTolerance = max(options.discreteTailCutoff, T.leastNonzeroMagnitude)

            let startIndex: Int
            if bounds.lower.isInfinite {
                return nil
            }
            else if let s = discreteCeil(bounds.lower) {
                startIndex = s
            }
            else {
                return nil
            }

            var endIndex: Int? = nil
            if bounds.upper.isFinite {
                endIndex = discreteFloor(bounds.upper)
                if let end = endIndex, end < startIndex {
                    return nil
                }
            }

            var idx = startIndex
            var divergence = T.zero
            var iterations = 0
            var infiniteFlag = false
            while true {
                if let end = endIndex, idx > end { break }
                let point = T(idx)
                let p = positiveDensity(self, at: point)
                if p > densityFloor {
                    var q = positiveDensity(other, at: point)
                    if !q.isFinite || q <= .zero {
                        infiniteFlag = true
                    } else {
                        if q < densityFloor { q = densityFloor }
                        let ratio = p / q
                        if ratio.isFinite && ratio > .zero {
                            divergence += p * T.log(ratio)
                        }
                    }
                }

                iterations += 1
                if iterations >= options.maxDiscreteEvaluations {
                    return nil
                }

                if endIndex == nil {
                    let tailSelf = positiveSurvival(self, above: point)
                    let tailOther = positiveSurvival(other, above: point)
                    if tailSelf <= tailTolerance && tailOther <= tailTolerance {
                        break
                    }
                }

                idx += 1
            }

            if infiniteFlag {
                return T.infinity
            }
            return divergence
        }

        private func sharedBounds(
            with other: Distribution.Dynamic<T>
        ) -> (lower: T, upper: T)? {
            let lower = Swift.max(supportLowerBound, other.supportLowerBound)
            let upper = Swift.min(supportUpperBound, other.supportUpperBound)
            if lower.isNaN || upper.isNaN || lower > upper {
                return nil
            }
            return (lower, upper)
        }

        @inline(__always)
        private func positiveDensity(
            _ distribution: Distribution.Dynamic<T>,
            at x: T
        ) -> T {
            guard let value = try? distribution.pdf(x), value.isFinite, value > 0 else {
                return .zero
            }
            return value
        }

        @inline(__always)
        private func positiveSurvival(
            _ distribution: Distribution.Dynamic<T>,
            above x: T
        ) -> T {
            guard let value = try? distribution.sf(x), value.isFinite else {
                return T.infinity
            }
            return value < 0 ? T.infinity : value
        }

        private func discreteCeil(_ value: T) -> Int? {
            let dv = Double(value)
            guard dv.isFinite else { return nil }
            if dv >= Double(Int.max) { return nil }
            if dv <= Double(Int.min) { return nil }
            return Int(Foundation.ceil(dv))
        }

        private func discreteFloor(_ value: T) -> Int? {
            let dv = Double(value)
            guard dv.isFinite else { return nil }
            if dv >= Double(Int.max) { return Int(Int.max) }
            if dv <= Double(Int.min) { return nil }
            return Int(Foundation.floor(dv))
        }

        private func analyticKLDivergence(
            relativeTo other: Distribution.Dynamic<T>
        ) -> T? {
            if Self.isNormalAlias(nameNormalized),
               Self.isNormalAlias(other.nameNormalized),
               let paramsSelf = normalParamValues(),
               let paramsOther = other.normalParamValues()
            {
                let (meanSelf, sdSelf) = paramsSelf
                let (meanOther, sdOther) = paramsOther
                guard sdSelf > 0, sdOther > 0 else { return nil }
                let logTerm = T.log(sdOther / sdSelf)
                let meanDiff = meanSelf - meanOther
                let varianceTerm = (sdSelf * sdSelf + meanDiff * meanDiff) / (2 * sdOther * sdOther)
                return logTerm + varianceTerm - 0.5
            }
            return nil
        }

        private func normalParamValues() -> (mean: T, sd: T)? {
            let meanKeys = ["mean", "mu", "location"]
            let sdKeys = ["sd", "sigma", "stddev", "standard_deviation"]
            guard let sd = firstParam(sdKeys), sd > 0 else { return nil }
            let mean = firstParam(meanKeys) ?? 0
            return (mean, sd)
        }

        /// Derived lattice metadata (origin and step) for discrete distributions.
        private var discreteLattice: (origin: T, step: T)? {
            switch nameNormalized {
            case "bernoulli", "bernoulli_distribution":
                return (0, 1)
            case "binomial", "binomial_distribution":
                return (0, 1)
            case "negative_binomial", "negative_binomial_distribution", "negativebinomial", "neg_binomial", "nbinom":
                return (0, 1)
            case "geometric", "geometric_distribution":
                return (0, 1)
            case "hypergeometric", "hypergeometric_distribution":
                return (0, 1)
            case "poisson", "poisson_distribution":
                return (0, 1)
            default:
                return nil
            }
        }

        // MARK: - Swift fallbacks for missing metrics

        private func fallbackEntropy() -> T? {
            switch nameNormalized {
            default:
                return nil
            }
        }

        private func firstParam(_ keys: [String]) -> T? {
            for k in keys {
                if let v = paramsLower[k] {
                    return v
                }
            }
            return nil
        }

        private static func isNormalAlias(_ name: String) -> Bool {
            switch name {
            case "normal", "normal_distribution", "gauss", "gaussian", "gaussian_distribution", "gauss_distribution":
                return true
            default:
                return false
            }
        }
    }
}
