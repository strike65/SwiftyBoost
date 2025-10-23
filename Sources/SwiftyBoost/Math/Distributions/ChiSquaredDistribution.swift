//
//  Created by Volker Thieme
//  Copyright © 2025 Volker Thieme. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
import SwiftyBoostPrelude

extension Distribution {
    /// A Chi-squared probability distribution.
    ///
    /// The Chi-squared distribution with ν degrees of freedom is the distribution of a sum of the squares of ν independent standard normal random variables.
    /// It is a special case of the Gamma distribution with shape k = ν/2 and scale θ = 2.
    ///
    /// - Generic parameter `T`: The floating-point type used by this distribution (`Float`, `Double`, or platform-available `Float80`).
    ///
    /// - Important: All methods are value-type and thread-safe. This type is a thin Swift wrapper around a dynamic backend
    ///   provided by Boost via `CBoostBridge`, while also exposing some analytically derived properties implemented in Swift.
    ///
    /// - SeeAlso: Wikipedia: Chi-squared distribution, Boost.Math Chi-squared distribution
    public struct ChiSquared<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The floating-point type used by this distribution.
        typealias Real = T

        /// The degrees of freedom ν (> 0).
        ///
        /// - Note: Must be finite and strictly greater than zero.
        public let degreesOfFreedom: T

        /// The dynamic backend that provides core distribution functionality (PDF/CDF/quantiles, etc.).
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Chi-squared distribution with the given degrees of freedom.
        ///
        /// - Parameters:
        ///   - df: The degrees of freedom ν. Must be finite and greater than `0`.
        ///
        /// - Throws: `DistributionError.invalidCombination` if `df <= 0`. `DistributionError.parameterNotFinite` if `df` is not finite.
        ///
        /// - Example:
        ///   ```swift
        ///   let chi2 = try Distribution.ChiSquared<Double>(degreesOfFreedom: 10)
        ///   let p = try chi2.cdf(15.0)
        ///   ```
        public init(degreesOfFreedom df: T) throws {
            guard df > 0 else {
                throw DistributionError.invalidCombination(message: "deegreesOfFreedom must be > 0", value: df)
            }
            guard df.isFinite else {
                throw DistributionError.parameterNotFinite(name: "df", value: df) }
            self.degreesOfFreedom = df
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "chisquared",
                parameters: [
                    "df": df,
                ]
            )
        }

        /// The lower bound of the support.
        ///
        /// For the Chi-squared distribution this is typically `0`.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support.
        ///
        /// For the Chi-squared distribution this is `+∞`.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// The closed-open range of the support as a tuple `(lower, upper)`.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: - Core functions

        /// The probability density function (PDF) evaluated at `x`.
        ///
        /// - Parameter x: Point at which to evaluate the PDF. For `x < 0` the PDF is `0`.
        /// - Returns: `f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// The natural logarithm of the PDF evaluated at `x`.
        ///
        /// - Parameter x: Point at which to evaluate the log-PDF.
        /// - Returns: `log(f(x))`.
        /// - Throws: Rethrows any underlying backend error.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// The cumulative distribution function (CDF) evaluated at `x`.
        ///
        /// - Parameter x: Point at which to evaluate the CDF.
        /// - Returns: `P(X ≤ x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// The survival function (SF) evaluated at `x`.
        ///
        /// - Parameter x: Point at which to evaluate the SF.
        /// - Returns: `P(X > x) = 1 - CDF(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// The quantile function (inverse CDF) evaluated at probability `p`.
        ///
        /// - Parameter p: A probability in `[0, 1]`.
        /// - Returns: The value `x` such that `P(X ≤ x) = p`.
        /// - Throws: `DistributionError.domain` if `p` is outside `[0, 1]`, or underlying backend errors.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// The quantile of the complement distribution (inverse SF) evaluated at probability `q`.
        ///
        /// - Parameter q: A probability in `[0, 1]` representing the upper tail probability.
        /// - Returns: The value `x` such that `P(X > x) = q`.
        /// - Throws: `DistributionError.domain` if `q` is outside `[0, 1]`, or underlying backend errors.
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        /// The mean of the distribution, if defined.
        ///
        /// For Chi-squared, `E[X] = ν`.
        public var mean: T? { dyn.mean }

        /// The variance of the distribution, if defined.
        ///
        /// For Chi-squared, `Var[X] = 2ν`.
        public var variance: T? { dyn.variance }

        /// The mode of the distribution, if defined.
        ///
        /// For Chi-squared, `mode = max(ν - 2, 0)`.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// There is no simple closed form; this value is computed numerically by the backend.
        public var median: T { dyn.median }

        /// The skewness of the distribution, if defined.
        ///
        /// For Chi-squared, `skewness = sqrt(8/ν)`.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis of the distribution, if defined.
        ///
        /// For Chi-squared, `kurtosis = 12/ν + 3`.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution, if defined.
        ///
        /// For Chi-squared, `excess = 12/ν`.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        /// The hazard (failure) rate function `h(x) = f(x) / S(x)`.
        ///
        /// - Parameter x: Point at which to evaluate the hazard function.
        /// - Returns: The hazard value at `x`.
        /// - Throws: Rethrows any underlying backend error.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// The cumulative hazard function `H(x) = -ln(S(x))`.
        ///
        /// - Parameter x: Point at which to evaluate the cumulative hazard function.
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: Rethrows any underlying backend error.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        /// The lattice step of the distribution if it is discrete.
        ///
        /// Chi-squared is continuous, so this is `nil`.
        public var latticeStep: T? { nil }

        /// The lattice origin of the distribution if it is discrete.
        ///
        /// Chi-squared is continuous, so this is `nil`.
        public var latticeOrigin: T? { nil }

        /// The (differential) entropy of the Chi-squared distribution, if it exists.
        ///
        /// For ν > 0 the entropy is:
        ///
        ///     H = (ν/2) + ln(2) + ln Γ(ν/2) + (1 - ν/2) ψ(ν/2)
        ///
        /// where `Γ` is the Gamma function and `ψ` is the digamma function.
        ///
        /// - Returns: The entropy as `T`,⁄ or `nil` if it cannot be computed (e.g., special functions failed).
        ///
        /// - Note: This implementation computes the formula directly using `SpecialFunctions.logGamma` and `SpecialFunctions.digamma`.
        ///   It will choose an appropriate `log` overload for the generic `T`.
        public var entropy: T? {
            guard degreesOfFreedom > 0 else { return nil }
            do {
                let lgA: T = try SpecialFunctions.logGamma(self.degreesOfFreedom / 2)
                let dgA: T = try SpecialFunctions.digamma(self.degreesOfFreedom / 2)
                let term1: T
                if T.self == Double.self {
                    term1 = self.degreesOfFreedom / 2 + T(log(Double(2)))
                } else if T.self == Float.self {
                    term1 = self.degreesOfFreedom / 2 + T(log(Float(2)))
                } else {
                    #if arch(x86_64) || arch(i386)
                    term1 = self.degreesOfFreedom / 2 + T(log(Float80(2)))
                    #else
                    term1 = self.degreesOfFreedom / 2 + T(log(Double(2)))
                    #endif
                }
                let term2: T = lgA + (1 - (self.degreesOfFreedom / 2)) * dgA
                return term1 + term2
            } catch {
                return nil
            }
        }

        /// Find a degrees-of-freedom value that achieves a desired difference-from-variance criterion.
        ///
        /// This is a helper that delegates to architecture-appropriate implementations in `CBoostBridge`.
        /// It can be used to solve for ν given a target relationship between variance and other parameters,
        /// for example in certain power or sample-size style calculations.
        ///
        /// - Parameters:
        ///   - difference_from_variance: The target difference from the variance (definition depends on the underlying routine).
        ///   - alpha: A tuning parameter (routine-specific).
        ///   - beta: A tuning parameter (routine-specific).
        ///   - variance: The variance value to compare against (routine-specific).
        ///   - hint: An optional initial guess for ν, default is `100`.
        ///
        /// - Returns: The estimated degrees of freedom ν as `T`.
        ///
        /// - Important: The precise interpretation of `difference_from_variance`, `alpha`, `beta`, and `variance`
        ///   follows the semantics of the Boost-based routine. Consult the corresponding C/Boost documentation
        ///   for details about how these inputs are used.
        ///
        /// - Example:
        ///   ```swift
        ///   // Solve for ν (Double) using an initial hint of 50
        ///   let nu = Distribution.ChiSquared<Double>.find_degreesOfFreedom(
        ///       difference_from_variance: 0.1,
        ///       alpha: 0.05,
        ///       beta: 0.2,
        ///       variance: 4.0,
        ///       hint: 50
        ///   )
        ///   ```
        public static func find_degreesOfFreedom(
            difference_from_variance: T,
            alpha: T,
            beta: T,
            variance: T,
            hint: T = 100
        ) -> Int {
            let res: T
            if T.self == Double.self {
                res = T(
                    bs_chisquare_find_degreesOfFreedom_d(
                        Double(difference_from_variance),
                        Double(alpha),
                        Double(beta),
                        Double(variance),
                        Double(hint)
                    )
                )
            }
            else if T.self == Float.self {
                res = T(
                    bs_chisquare_find_degreesOfFreedom_f(
                        Float(difference_from_variance),
                        Float(alpha),
                        Float(beta),
                        Float(variance),
                        Float(hint)
                    )
                )
            }
            else {
#if arch(x86_64) || arch(i386)
                res = T(
                    bs_chisquare_find_degreesOfFreedom_l(
                        Float80(difference_from_variance),
                        Float80(alpha),
                        Float80(beta),
                        Float80(variance),
                        Float80(hint)
                    )
                )
#else
                res = T(
                    bs_chisquare_find_degreesOfFreedom_d(
                        Double(difference_from_variance),
                        Double(alpha),
                        Double(beta),
                        Double(variance),
                        Double(hint)
                    )
                )
#endif
            }
            return Int(res) + 1
        }
    }
}
