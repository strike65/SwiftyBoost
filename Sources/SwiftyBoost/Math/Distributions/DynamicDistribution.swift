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

        private final class Box: @unchecked Sendable {
            var d: bs_dist_d?
            var f: bs_dist_f?
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

        private let box: Box
        private let nameNormalized: String
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
                    throw DistributionError.invalidCombination(message: "Unknown distribution or bad parameters: \(name)")
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
                    throw DistributionError.invalidCombination(message: "Unknown distribution or bad parameters: \(name)")
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
                    throw DistributionError.invalidCombination(message: "Unknown distribution or bad parameters: \(name)")
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
                    throw DistributionError.invalidCombination(message: "Unknown distribution or bad parameters: \(name)")
                }
                self.box = Box(d: ptr.pointee)
                #endif
            }
        }

        // MARK: - DistributionProtocol conformance

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

        public var range: (lower: T, upper: T) {
            (supportLowerBound, supportUpperBound)
        }

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

        public var latticeStep: T? { nil }

        public var latticeOrigin: T? { nil }

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
    }
}
