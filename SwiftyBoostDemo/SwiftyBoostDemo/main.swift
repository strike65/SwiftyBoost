//
//  Created by VT on 11.10.25.
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

import Foundation
import SwiftyBoost
// import CBoostBridge

// Example pattern (instantiation + property access):
// let b = try Distribution.Bernoulli(p: 0.5)
// print("\(String(describing: b.kurtosis))")

do {
    // Bernoulli(p)
    let b = try Distribution.Bernoulli(p: 0.5)
    print("Bernoulli:")
    print("  stored: p(success_fraction) = \(b.success_fraction)")
    print("  support: [\(b.supportLowerBound), \(b.supportUpperBound)]  range: \(b.range)")
    print("  mean = \(String(describing: b.mean))")
    print("  variance = \(String(describing: b.variance))")
    print("  mode = \(String(describing: b.mode))")
    print("  median = \(b.median)")
    print("  skewness = \(String(describing: b.skewness))")
    print("  kurtosis = \(String(describing: b.kurtosis))")
    print("  kurtosisExcess = \(String(describing: b.kurtosisExcess))")
    print("  entropy = \(String(describing: b.entropy))")
    print("  latticeStep = \(String(describing: b.latticeStep))")
    print("  latticeOrigin = \(String(describing: b.latticeOrigin))")
    print("  pdf(1) = \(try b.pdf(1.0))  cdf(1) = \(try b.cdf(1.0))  sf(0) = \(try b.sf(0.0))")
    print("  quantile(0.5) = \(try b.quantile(0.5))  q^c(0.5) = \(try b.quantileComplement(0.5))")
    print("  hazard(1) = \(try b.hazard(1.0))  chf(1) = \(try b.chf(1.0))")
} catch {
    print("Bernoulli error: \(error)")
}

do {
    // Arcsine(minX, maxX)
    let a = try Distribution.Arcsine(minX: 0.0, maxX: 10.0)
    print("Arcsine:")
    print("  stored: minX = \(a.minX), maxX = \(a.maxX)")
    print("  support: [\(a.supportLowerBound), \(a.supportUpperBound)]  range: \(a.range)")
    print("  mean = \(String(describing: a.mean))")
    print("  variance = \(String(describing: a.variance))")
    print("  mode = \(String(describing: a.mode))")
    print("  median = \(a.median)")
    print("  skewness = \(String(describing: a.skewness))")
    print("  kurtosis = \(String(describing: a.kurtosis))")
    print("  kurtosisExcess = \(String(describing: a.kurtosisExcess))")
    print("  entropy = \(String(describing: a.entropy))")
    print("  latticeStep = \(String(describing: a.latticeStep))")
    print("  latticeOrigin = \(String(describing: a.latticeOrigin))")
    print("  pdf(0.5) = \(try a.pdf(0.5))  cdf(0.5) = \(try a.cdf(0.5))  sf(0.5) = \(try a.sf(0.5))")
    print("  quantile(0.25) = \(try a.quantile(0.25))  q^c(0.25) = \(try a.quantileComplement(0.25))")
    print("  hazard(0.5) = \(try a.hazard(0.5))  chf(0.5) = \(try a.chf(0.5))")
} catch {
    print("Arcsine error: \(error)")
}

do {
    // Gamma(shape, scale)
    let g = try Distribution.Gamma(shape: 2.0, scale: 3.0)
    print("Gamma:")
    print("  stored: shape = \(g.shape), scale = \(g.scale)")
    print("  support: [\(g.supportLowerBound), \(g.supportUpperBound)]  range: \(g.range)")
    print("  mean = \(String(describing: g.mean))")
    print("  variance = \(String(describing: g.variance))")
    print("  mode = \(String(describing: g.mode))")
    print("  median = \(g.median)")
    print("  skewness = \(String(describing: g.skewness))")
    print("  kurtosis = \(String(describing: g.kurtosis))")
    print("  kurtosisExcess = \(String(describing: g.kurtosisExcess))")
    print("  entropy = \(String(describing: g.entropy))")
    print("  latticeStep = \(String(describing: g.latticeStep))")
    print("  latticeOrigin = \(String(describing: g.latticeOrigin))")
    print("  pdf(1) = \(try g.pdf(1.0))  cdf(1) = \(try g.cdf(1.0))  sf(1) = \(try g.sf(1.0))")
    print("  quantile(0.9) = \(try g.quantile(0.9))  q^c(0.1) = \(try g.quantileComplement(0.1))")
    print("  hazard(1) = \(try g.hazard(1.0))  chf(1) = \(try g.chf(1.0))")
} catch {
    print("Gamma error: \(error)")
}

do {
    // StudentT(df)
    let t = try Distribution.StudentT(degreesOfFreedom: 5.0)
    print("StudentT:")
    print("  stored: degreesOfFreedom = \(t.degreesOfFreedom)")
    print("  support: [\(t.supportLowerBound), \(t.supportUpperBound)]  range: \(t.range)")
    print("  mean = \(String(describing: t.mean))")
    print("  variance = \(String(describing: t.variance))")
    print("  mode = \(String(describing: t.mode))")
    print("  median = \(t.median)")
    print("  skewness = \(String(describing: t.skewness))")
    print("  kurtosis = \(String(describing: t.kurtosis))")
    print("  kurtosisExcess = \(String(describing: t.kurtosisExcess))")
    print("  entropy = \(String(describing: t.entropy))")
    print("  latticeStep = \(String(describing: t.latticeStep))")
    print("  latticeOrigin = \(String(describing: t.latticeOrigin))")
    print("  pdf(0) = \(try t.pdf(0.0))  cdf(0) = \(try t.cdf(0.0))  sf(0) = \(try t.sf(0.0))")
    print("  quantile(0.95) = \(try t.quantile(0.95))  q^c(0.05) = \(try t.quantileComplement(0.05))")
    print("  hazard(0) = \(try t.hazard(0.0))  chf(0) = \(try t.chf(0.0))")
} catch {
    print("StudentT error: \(error)")
}
