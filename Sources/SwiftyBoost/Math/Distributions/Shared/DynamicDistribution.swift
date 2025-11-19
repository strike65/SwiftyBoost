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
        
        public func klDivergence<D: DistributionProtocol>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D.RealType == T {
            if options.integrationLowerBound == nil,
               options.integrationUpperBound == nil,
               let otherDynamic = other as? Distribution.Dynamic<T>,
               let analytic = analyticKLDivergence(relativeTo: otherDynamic)
            {
                return analytic
            }
            return try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
        
        /// Boolean indicating whether the distribution operates on a discrete lattice.
        public var isDiscrete: Bool { discreteLattice != nil }

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
