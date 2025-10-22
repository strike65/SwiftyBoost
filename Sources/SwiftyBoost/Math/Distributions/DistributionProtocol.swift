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
//

/// A common interface for real-valued probability distributions.
///
/// Conforming types provide:
/// - The distribution's support (domain where the density/mass is nonzero).
/// - Core functions such as the PDF, log-PDF, CDF, survival function, hazard, and cumulative hazard.
/// - Inverse functions (quantiles) for both lower- and upper-tail probabilities.
/// - Summary statistics such as mean, variance, median, skewness, and kurtosis.
///
/// Unless otherwise noted, methods are expected to operate on values within the
/// distribution's support. Calling them with values outside the support may throw
/// an error specific to the conforming type (for example, a domain error).

import SwiftyBoostPrelude
internal protocol DistributionProtocol {
    /// The real number type used by this distribution (for example, `Float`, `Double`, `Float80`).
    ///
    /// Conforming types should choose a floating-point type that suits their precision
    /// and performance needs. The associated type must be `Sendable` so values can be
    /// safely passed across concurrency domains.
    associatedtype RealType: Real & BinaryFloatingPoint & Sendable

    // MARK: - Support

    /// The lower bound of the distribution's support (domain where the distribution is defined).
    ///
    /// Depending on the distribution, the lower bound may be finite or infinite (represented
    /// by `-infinity` if supported by `Real`). Whether the bound is inclusive or exclusive
    /// is distribution-specific and should be documented by the conforming type.
    var supportLowerBound: RealType { get }

    /// The upper bound of the distribution's support (domain where the distribution is defined).
    ///
    /// Depending on the distribution, the upper bound may be finite or infinite (represented
    /// by `+infinity` if supported by `Real`). Whether the bound is inclusive or exclusive
    /// is distribution-specific and should be documented by the conforming type.
    var supportUpperBound: RealType { get }

    // MARK: - Core functions

    /// The cumulative distribution function (CDF).
    ///
    /// For continuous distributions, this is F(x) = P(X ≤ x).
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The CDF value in [0, 1].
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func cdf(_ x: RealType) throws -> RealType

    /// The natural logarithm of the probability density function (log-PDF).
    ///
    /// Implementations should prefer numerically stable computations and avoid computing
    /// `log(pdf(x))` directly where possible.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The log-PDF at `x`.
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func logPdf(_ x: RealType) throws -> RealType

    /// The probability density function (PDF).
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The PDF at `x` (nonnegative, integrating to 1 over the support).
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func pdf(_ x: RealType) throws -> RealType

    /// The survival function (SF), also known as the complementary CDF.
    ///
    /// For continuous distributions, this is S(x) = P(X > x). Where numerically advantageous,
    /// implementations should compute this directly rather than using `1 - cdf(x)`.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The survival probability in [0, 1].
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func sf(_ x: RealType) throws -> RealType

    /// The (instantaneous) hazard function.
    ///
    /// For continuous distributions, h(x) = f(x) / S(x) where defined.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The hazard at `x`.
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func hazard(_ x: RealType) throws -> RealType

    /// The cumulative hazard function (CHF).
    ///
    /// For continuous distributions, H(x) = −log S(x) where defined.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The cumulative hazard at `x`.
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func chf(_ x: RealType) throws -> RealType

    // MARK: - Inverses

    /// The lower-tail quantile function, i.e. the inverse CDF.
    ///
    /// Returns the value `x` such that F(x) = p for `p` in [0, 1].
    /// Implementations should be monotonic in `p` and numerically stable near 0 and 1.
    ///
    /// - Parameter p: A probability in [0, 1].
    /// - Returns: The quantile corresponding to probability `p`.
    /// - Throws: If `p` is outside [0, 1] or if the computation fails to converge.
    func quantile(_ p: RealType) throws -> RealType

    /// The upper-tail quantile function, i.e. the inverse survival function.
    ///
    /// Returns the value `x` such that S(x) = p for `p` in [0, 1].
    /// Implementations should be monotonic in `p` and numerically stable near 0 and 1.
    ///
    /// - Parameter p: An upper-tail probability in [0, 1].
    /// - Returns: The upper-tail quantile corresponding to probability `p`.
    /// - Throws: If `p` is outside [0, 1] or if the computation fails to converge.
    func quantileComplement(_ p: RealType) throws -> RealType

    // MARK: - Summary statistics

    /// The expected value (mean) of the distribution.
    ///
    /// If the mean does not exist (diverges), conforming types should document and choose
    /// an appropriate behavior (for example, `±infinity`, `nan`, or `nil`).
    var mean: RealType? { get }

    /// The mode (value at which the PDF attains its maximum), if it exists and is unique.
    ///
    /// Some distributions do not have a unique mode; in such cases, this may be `nil`.
    var mode: RealType? { get }

    /// The variance of the distribution.
    ///
    /// If the variance does not exist (diverges), conforming types should document and choose
    /// an appropriate behavior (for example, `±infinity`, `nan`, or `nil`).
    var variance: RealType? { get }

    /// The median of the distribution.
    ///
    /// For continuous unimodal distributions this is typically unique; for others,
    /// a conventional choice should be documented by the conforming type.
    var median: RealType { get }

    /// The skewness of the distribution.
    ///
    /// Defined as E[((X − μ)/σ)^3] where it exists.
    var skewness: RealType? { get }

    /// The kurtosis of the distribution (Pearson's definition).
    ///
    /// Defined as E[((X − μ)/σ)^4] where it exists. This value equals 3 for the normal distribution.
    var kurtosis: RealType? { get }

    /// The excess kurtosis of the distribution.
    ///
    /// Defined as `kurtosis − 3`. Equals 0 for the normal distribution.
    var kurtosisExcess: RealType? { get }

    /// A convenient representation of the distribution's overall range.
    ///
    /// Unless documented otherwise by the conforming type, this is typically equivalent
    /// to `(supportLowerBound, supportUpperBound)`, which may include infinite values.
    var range: (lower: RealType, upper: RealType) { get }

    /// For lattice (discrete) distributions, the spacing between adjacent support points.
    ///
    /// Continuous distributions should return `nil`.
    var latticeStep: RealType? { get }

    /// For lattice (discrete) distributions, the origin (offset) of the lattice.
    ///
    /// Continuous distributions should return `nil`.
    var latticeOrigin: RealType? { get }

    /// The entropy H[X] (Shannon) where defined.
    ///
    /// Implementations should return `nil` when the entropy is undefined or does not
    /// exist (for example, diverges). Units are in nats (natural logarithm base).
    var entropy: RealType? { get }
}
