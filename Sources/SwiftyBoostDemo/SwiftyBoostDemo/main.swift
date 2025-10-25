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
internal import Numerics
// import CBoostBridge

// Example pattern (instantiation + property access):
// let b = try Distribution.Bernoulli(p: 0.5)
// print("\(String(describing: b.kurtosis))")

let ja: Double = try SpecialFunctions.jacobiEllipticCd(0.4, theta: 4.5)

print("\(try SpecialFunctions.legendreStieltjes(12, 0.3))")
do {
    let Zeros: [Double] = try SpecialFunctions.legendreStieltjesZeros(order: 12)
    print("\(Zeros)")
}
catch let error {
    print("\(error)")
}
do {
    let Zeros: [Float] = try SpecialFunctions.legendreStieltjesZeros(order: 12)
    print("\(Zeros)")
}
catch let error {
    print("\(error)")
}

do {
    print("Quadrature examples:")
    
    let legendre = try SpecialFunctions.Quadrature.Integrator<Double>(rule: .gaussLegendre(points: 30))
    
    let reciprocalLinear = try legendre.integrate(over: .finite(lower: 0.0, upper: 1.0)) { x in 1.0 / (x + 1.0) }
    print("  Gauss-Legendre(30) integral[0,1] 1/(x+1) dx ≈ \(reciprocalLinear.value) (+/-\(reciprocalLinear.estimatedError))")
    
    let sineIntegral = try legendre.integrate(
        over: .finite(lower: 0.0, upper: Double.pi / 2.0 )
    ) {
        x in Double.sin(x)
    }
    print("  Gauss-Legendre(30) integral[0,pi / 2] sin(x) dx ≈ \(sineIntegral.value)")
    
    let inverseSine = try SpecialFunctions.Quadrature.integrate(
        using: SpecialFunctions.Quadrature.Rule.gaussKronrod(points: 21),
        over: .finite(lower: 0.0, upper: Double.pi)
    ) { (x: Double) in 1.0 / (1.0 + Double.sin(x)) }
    print("  Gauss-Kronrod(21) integral[0,pi] 1/(1+sin(x)) dx ≈ \(inverseSine.value)")
    
    let cauchy = try Distribution.Cauchy<Double>(location: -2.0, scale: 0.4)
    let loc = cauchy.location
    let scale = cauchy.scale
    let cauchyArea = try SpecialFunctions.Quadrature.integrate(
        using: SpecialFunctions.Quadrature.Rule.sinhSinh(),
        integrand: { (x: Double) in
            let z = (x - loc) / scale
            return 1.0 / (Double.pi * scale * (1.0 + z * z))
        }
    )
    print("  Sinh-Sinh integral(-inf,+inf) Cauchy1(loc:-2, scale:0.4) dx ≈ \(cauchyArea.value)")
    let cauchyArea1 = try SpecialFunctions.Quadrature.integrate(
        using: SpecialFunctions.Quadrature.Rule.sinhSinh(),
        integrand: { (x: Double) in try! cauchy.pdf(x)
        }
    )
    print("  Sinh-Sinh integral(-inf,+inf) Cauchy2(loc:-2, scale:0.4) dx ≈ \(cauchyArea.value)")
    let hm = try Distribution.Holtsmark(loc : 33, scale: 1.0)
    let hmArea = try SpecialFunctions.Quadrature.integrate(
        using: SpecialFunctions.Quadrature.Rule.gaussLegendre(points: 100),
        integrand: { (x: Double) in try! hm.pdf(x)
        }
    )
    print("  gaussLegendre(-inf,+inf) Holtsmark(scale:1 shape: 0.5) dx ≈ \(cauchyArea.value)")
    
    
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
    do {
        // ExtremeValue
        let e = try Distribution.ExtremeValueGumpel(location: 0.0, scale: 1.0)
        print("StudentT:")
        print("  stored: location = \(e.location)")
        print("  stored: scale = \(e.scale)")
        print("  support: [\(e.supportLowerBound), \(e.supportUpperBound)]  range: \(e.range)")
        print("  mean = \(String(describing: e.mean))")
        print("  variance = \(String(describing: e.variance))")
        print("  mode = \(String(describing: e.mode))")
        print("  median = \(e.median)")
        print("  skewness = \(String(describing: e.skewness))")
        print("  kurtosis = \(String(describing: e.kurtosis))")
        print("  kurtosisExcess = \(String(describing: e.kurtosisExcess))")
        print("  entropy = \(String(describing: e.entropy))")
        print("  latticeStep = \(String(describing: e.latticeStep))")
        print("  latticeOrigin = \(String(describing: e.latticeOrigin))")
        print("  pdf(0) = \(try e.pdf(0.0))  cdf(0) = \(try e.cdf(0.0))  sf(0) = \(try e.sf(0.0))")
        print("  quantile(0.95) = \(try e.quantile(0.95))  q^c(0.05) = \(try e.quantileComplement(0.05))")
        print("  hazard(0) = \(try e.hazard(0.0))  chf(0) = \(try e.chf(0.0))")
    } catch {
        print("Gumpel error: \(error)")
    }
    
    do {
        // Geometric
        let t = try Distribution.Geometric(probabibilityOfSuccess: 0.4)
        print("StudentT:")
        print(
            "  stored: probabilityOfSuccess = \(t.pSuccess)"
        )
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
        print("ExtremeValue error: \(error)")
    }
    
    do {
        // Binomial(n, p)
        let bin = try Distribution.Binomial(numberOfTrials: 10.0, probabibilityOfSuccess: 0.3)
        print("Binomial:")
        print("  stored: nTrials = \(bin.nTrials), pSuccess = \(bin.pSuccess)")
        print("  support: [\(bin.supportLowerBound), \(bin.supportUpperBound)]  range: \(bin.range)")
        print("  mean = \(String(describing: bin.mean))")
        print("  variance = \(String(describing: bin.variance))")
        print("  mode = \(String(describing: bin.mode))")
        print("  median = \(bin.median)")
        print("  skewness = \(String(describing: bin.skewness))")
        print("  kurtosis = \(String(describing: bin.kurtosis))")
        print("  kurtosisExcess = \(String(describing: bin.kurtosisExcess))")
        print("  entropy = \(String(describing: bin.entropy))")
        print("  latticeStep = \(String(describing: bin.latticeStep))")
        print("  latticeOrigin = \(String(describing: bin.latticeOrigin))")
        print("  pdf(3) = \(try bin.pdf(3.0))  logpdf(3) = \(try bin.logPdf(3.0))")
        print("  cdf(3) = \(try bin.cdf(3.0))  sf(3) = \(try bin.sf(3.0))")
        print("  quantile(0.5) = \(try bin.quantile(0.5))  q^c(0.1) = \(try bin.quantileComplement(0.1))")
        print("  hazard(3) = \(try bin.hazard(3.0))  chf(3) = \(try bin.chf(3.0))")
    } catch {
        print("Binomial error: \(error)")
    }
    
    do {
        // Cauchy)
        let ca = try Distribution.Cauchy(location: 0.0, scale: 2)
        print("Cauchy:")
        print("  stored: location = \(ca.location)")
        print("  stored: scale = \(ca.scale)")
        print("  support: [\(ca.supportLowerBound), \(ca.supportUpperBound)]  range: \(ca.range)")
        print("  mean = \(String(describing: ca.mean))")
        print("  variance = \(String(describing: ca.variance))")
        print("  mode = \(String(describing: ca.mode))")
        print("  median = \(ca.median)")
        print("  skewness = \(String(describing: ca.skewness))")
        print("  kurtosis = \(String(describing: ca.kurtosis))")
        print("  kurtosisExcess = \(String(describing: ca.kurtosisExcess))")
        print("  entropy = \(String(describing: ca.entropy))")
        print("  latticeStep = \(String(describing: ca.latticeStep))")
        print("  latticeOrigin = \(String(describing: ca.latticeOrigin))")
        print("  pdf(3) = \(try ca.pdf(3.0))  logpdf(3) = \(try ca.logPdf(3.0))")
        print("  cdf(3) = \(try ca.cdf(3.0))  sf(3) = \(try ca.sf(3.0))")
        print("  quantile(0.5) = \(try ca.quantile(0.5))  q^c(0.1) = \(try ca.quantileComplement(0.1))")
        print("  hazard(3) = \(try ca.hazard(3.0))  chf(3) = \(try ca.chf(3.0))")
    } catch {
        print("Cauchy error: \(error)")
    }
    do {
        // Holtsmark
        let ca = try Distribution.Holtsmark(loc: 1.0, scale: 0.2)
        print("Holtsmark:")
        print("  stored: location = \(ca.location)")
        print("  stored: scale = \(ca.scale)")
        print("  support: [\(ca.supportLowerBound), \(ca.supportUpperBound)]  range: \(ca.range)")
        print("  mean = \(String(describing: ca.mean))")
        print("  variance = \(String(describing: ca.variance))")
        print("  mode = \(String(describing: ca.mode))")
        print("  median = \(ca.median)")
        print("  skewness = \(String(describing: ca.skewness))")
        print("  kurtosis = \(String(describing: ca.kurtosis))")
        print("  kurtosisExcess = \(String(describing: ca.kurtosisExcess))")
        print("  entropy = \(String(describing: ca.entropy))")
        print("  latticeStep = \(String(describing: ca.latticeStep))")
        print("  latticeOrigin = \(String(describing: ca.latticeOrigin))")
        print("  pdf(3) = \(try ca.pdf(3.0))  logpdf(3) = \(try ca.logPdf(3.0))")
        print("  cdf(3) = \(try ca.cdf(3.0))  sf(3) = \(try ca.sf(3.0))")
        print("  quantile(0.5) = \(try ca.quantile(0.5))  q^c(0.1) = \(try ca.quantileComplement(0.1))")
        print("  hazard(3) = \(try ca.hazard(3.0))  chf(3) = \(try ca.chf(3.0))")
    } catch {
        print("Holtsmark error: \(error)")
    }
    
    
    // MARK: - Distribution static helpers
    
    do {
        print("Distribution static functions:")
        // StudentT.findDegreesOfFreedom
        do {
            print("  StudentT.findDegreesOfFreedom:")
            let diffD: Double = 1.0
            let alphaD: Double = 0.05
            let betaD: Double = 0.2
            let sdD: Double = 1.5
            let hintD: Double = 2.0
            let nuD = Distribution.StudentT<Double>.findDegreesOfFreedom(
                differenceFromMean: diffD,
                alpha: alphaD,
                beta: betaD,
                sd: sdD,
                hint: hintD
            )
            print("    Double -> ν = \(nuD)")
            
            let diffF: Float = 0.8
            let alphaF: Float = 0.05
            let betaF: Float = 0.2
            let sdF: Float = 1.25
            let hintF: Float = 2.0
            let nuF = Distribution.StudentT<Float>.findDegreesOfFreedom(
                differenceFromMean: diffF,
                alpha: alphaF,
                beta: betaF,
                sd: sdF,
                hint: hintF
            )
            print("    Float -> ν = \(nuF)")
            
#if arch(x86_64)
            let diffL: Float80 = 0.9
            let alphaL: Float80 = 0.05
            let betaL: Float80 = 0.2
            let sdL: Float80 = 1.1
            let hintL: Float80 = 2.0
            let nuL = Distribution.StudentT<Float80>.findDegreesOfFreedom(
                differenceFromMean: diffL,
                alpha: alphaL,
                beta: betaL,
                sd: sdL,
                hint: hintL
            )
            print("    Float80 -> ν = \(nuL)")
#endif
        }
        
        // ChiSquared.find_degreesOfFreedom
        do {
            print("  ChiSquared.find_degreesOfFreedom:")
            let dfvD = Distribution.ChiSquared<Double>.find_degreesOfFreedom(
                difference_from_variance: 0.1,
                alpha: 0.05,
                beta: 0.2,
                variance: 4.0,
                hint: 50.0
            )
            print("    Double -> ν = \(dfvD)")
            
            let dfvF = Distribution.ChiSquared<Float>.find_degreesOfFreedom(
                difference_from_variance: 0.1,
                alpha: 0.05,
                beta: 0.2,
                variance: 3.0,
                hint: 25.0
            )
            print("    Float -> ν = \(dfvF)")
            
#if arch(x86_64)
            let dfvL = Distribution.ChiSquared<Float80>.find_degreesOfFreedom(
                difference_from_variance: 0.05,
                alpha: 0.05,
                beta: 0.2,
                variance: 5.0,
                hint: 75.0
            )
            print("    Float80 -> ν = \(dfvL)")
#endif
        }
        
        // Binomial planning helpers
        do {
            print("  Binomial planning helpers:")
            // One-sided bounds on p (alpha in [0,1])
            let n = 50, k = 20
            let alphaD: Double = 0.05
            let lbD = Distribution.Binomial<Double>.findLowerBoundOnP(nTrials: n, nSuccesses: k, proposedSuccessFraction: alphaD, useJeffreys: false)
            let ubD = Distribution.Binomial<Double>.findUpperBoundOnP(nTrials: n, nSuccesses: k, proposedSuccessFraction: alphaD, useJeffreys: false)
            print("    Double -> lowerBoundOnP = \(lbD), upperBoundOnP = \(ubD)")
            
            let alphaF: Float = 0.05
            let lbF = Distribution.Binomial<Float>.findLowerBoundOnP(nTrials: n, nSuccesses: k, proposedSuccessFraction: alphaF, useJeffreys: true)
            let ubF = Distribution.Binomial<Float>.findUpperBoundOnP(nTrials: n, nSuccesses: k, proposedSuccessFraction: alphaF, useJeffreys: true)
            print("    Float -> lowerBoundOnP (Jeffreys) = \(lbF), upperBoundOnP (Jeffreys) = \(ubF)")
            
            // Trial count solvers (alpha in [0,1], p0 in [0,1])
            do {
                let s = 5
                let p0D: Double = 0.3
                let nMinD = try Distribution.Binomial<Double>.findMinimumNumberOfTrials(successes: s, proposedSuccessFraction: p0D, alpha: 0.05)
                let nMaxD = try Distribution.Binomial<Double>.findMaximumNumberOfTrials(successes: s, proposedSuccessFraction: p0D, alpha: 0.05)
                print("    Double -> nMin = \(nMinD), nMax = \(nMaxD)")
            } catch {
                print("    Double trial-count solvers error: \(error)")
            }
            do {
                let s = 7
                let p0F: Float = 0.4
                let nMinF = try Distribution.Binomial<Float>.findMinimumNumberOfTrials(successes: s, proposedSuccessFraction: p0F, alpha: 0.1)
                let nMaxF = try Distribution.Binomial<Float>.findMaximumNumberOfTrials(successes: s, proposedSuccessFraction: p0F, alpha: 0.1)
                print("    Float -> nMin = \(nMinF), nMax = \(nMaxF)")
            } catch {
                print("    Float trial-count solvers error: \(error)")
            }
#if arch(x86_64)
            do {
                let s = 6
                let alphaL: Float80 = 0.05
                let lbL = Distribution.Binomial<Float80>.findLowerBoundOnP(nTrials: n, nSuccesses: k, proposedSuccessFraction: alphaL, useJeffreys: false)
                let ubL = Distribution.Binomial<Float80>.findUpperBoundOnP(nTrials: n, nSuccesses: k, proposedSuccessFraction: alphaL, useJeffreys: false)
                print("    Float80 -> lowerBoundOnP = \(lbL), upperBoundOnP = \(ubL)")
                let p0L: Float80 = 0.35
                let nMinL = try Distribution.Binomial<Float80>.findMinimumNumberOfTrials(successes: s, proposedSuccessFraction: p0L, alpha: alphaL)
                let nMaxL = try Distribution.Binomial<Float80>.findMaximumNumberOfTrials(successes: s, proposedSuccessFraction: p0L, alpha: alphaL)
                print("    Float80 -> nMin = \(nMinL), nMax = \(nMaxL)")
            } catch {
                print("    Float80 trial-count solvers error: \(error)")
            }
#endif
        }
    }
    
    // MARK: - Complex.swift API coverage
    
    // Double-precision complex
    do {
        print("Complex<Double> (ComplexD) API:")
        var z = ComplexD(3.0, 4.0)
        print("  z = \(z)")
        // imag getter/setter
        print("  real = \(z.real), imag = \(z.imag)")
        z.imag = 5.0
        print("  after set imag -> \(z)")
        // zeros/ones/i
        let z0 = ComplexD.zero
        let z1 = ComplexD.one
        let zi = ComplexD.i
        print("  zero = \(z0), one = \(z1), i = \(zi)")
        // basic props
        print("  conjugate = \(z.conjugate)")
        print("  norm = \(z.norm)")
        print("  magnitude = \(z.magnitude)")
        print("  isFinite = \(z.isFinite), isInfinite = \(z.isInfinite), isNaN = \(z.isNaN)")
        // reciprocal, normalized, pow, squareRoot, description
        print("  reciprocal = \(z.reciprocal())")
        print("  normalized = \(z.normalized())")
        print("  z.pow(3) = \(z.pow(3))")
        print("  squareRoot = \(z.squareRoot)")
        print("  description = \(z.description)")
        // mixed operators with scalar
        let s: Double = 2.5
        print("  z + s = \(z + s)")
        print("  s + z = \(s + z)")
        print("  z - s = \(z - s)")
        print("  s - z = \(s - z)")
        print("  z * s = \(z * s)")
        print("  s * z = \(s * z)")
        print("  z / s = \(z / s)")
        // complex division
        let w = ComplexD(1.5, -0.75)
        print("  z / w = \(z / w)")
        // Double specialization: fromPolar, phase, exp/log/trig/hyperbolic/atan
        let zp = ComplexD.fromPolar(radius: 2.0, phase: .pi / 3.0)
        print("  fromPolar(2, pi/3) = \(zp)")
        print("  phase(z) = \(z.phase)")
        print("  exp(z) = \(z.exp)")
        print("  log(z) = \(z.log)")
        print("  sin(z) = \(z.sin)")
        print("  cos(z) = \(z.cos)")
        print("  tan(z) = \(z.tan)")
        print("  sinh(z) = \(z.sinh)")
        print("  cosh(z) = \(z.cosh)")
        print("  tanh(z) = \(z.tanh)")
        print("  atan(z) = \(z.atan)")
    }
    
    // Float-precision complex
    do {
        print("Complex<Float> (ComplexF) API:")
        var zf = ComplexF(3.0, 4.0)
        print("  zf = \(zf)")
        // imag getter/setter
        print("  real = \(zf.real), imag = \(zf.imag)")
        zf.imag = 5.0
        print("  after set imag -> \(zf)")
        // zeros/ones/i
        let zf0 = ComplexF.zero
        let zf1 = ComplexF.one
        let zfi = ComplexF.i
        print("  zero = \(zf0), one = \(zf1), i = \(zfi)")
        // basic props
        print("  conjugate = \(zf.conjugate)")
        print("  norm = \(zf.norm)")
        print("  magnitude = \(zf.magnitude)")
        print("  isFinite = \(zf.isFinite), isInfinite = \(zf.isInfinite), isNaN = \(zf.isNaN)")
        // reciprocal, normalized, pow, squareRoot, description
        print("  reciprocal = \(zf.reciprocal())")
        print("  normalized = \(zf.normalized())")
        print("  zf.pow(3) = \(zf.pow(3))")
        print("  squareRoot = \(zf.squareRoot)")
        print("  description = \(zf.description)")
        // mixed operators with scalar
        let sf: Float = 2.5
        print("  zf + sf = \(zf + sf)")
        print("  sf + zf = \(sf + zf)")
        print("  zf - sf = \(zf - sf)")
        print("  sf - zf = \(sf - zf)")
        print("  zf * sf = \(zf * sf)")
        print("  sf * zf = \(sf * zf)")
        print("  zf / sf = \(zf / sf)")
        // complex division
        let wf = ComplexF(1.5, -0.75)
        print("  zf / wf = \(zf / wf)")
        // Float specialization: fromPolar, phase, exp/log/trig/hyperbolic/atan
        let zfp = ComplexF.fromPolar(radius: 2.0, phase: .pi / 3.0)
        print("  fromPolar(2, pi/3) = \(zfp)")
        print("  phase(zf) = \(zf.phase)")
        print("  exp(zf) = \(zf.exp)")
        print("  log(zf) = \(zf.log)")
        print("  sin(zf) = \(zf.sin)")
        print("  cos(zf) = \(zf.cos)")
        print("  tan(zf) = \(zf.tan)")
        print("  sinh(zf) = \(zf.sinh)")
        print("  cosh(zf) = \(zf.cosh)")
        print("  tanh(zf) = \(zf.tanh)")
        print("  atan(zf) = \(zf.atan)")
    }
    
#if arch(x86_64)
    // Float80-precision complex
    do {
        print("Complex<Float80> (ComplexL) API:")
        var zl = ComplexL(3.0, 4.0)
        print("  zl = \(zl)")
        // imag getter/setter
        print("  real = \(zl.real), imag = \(zl.imag)")
        zl.imag = 5.0
        print("  after set imag -> \(zl)")
        // zeros/ones/i
        let zl0 = ComplexL.zero
        let zl1 = ComplexL.one
        let zli = ComplexL.i
        print("  zero = \(zl0), one = \(zl1), i = \(zli)")
        // basic props
        print("  conjugate = \(zl.conjugate)")
        print("  norm = \(zl.norm)")
        print("  magnitude = \(zl.magnitude)")
        print("  isFinite = \(zl.isFinite), isInfinite = \(zl.isInfinite), isNaN = \(zl.isNaN)")
        // reciprocal, normalized, pow, squareRoot, description
        print("  reciprocal = \(zl.reciprocal())")
        print("  normalized = \(zl.normalized())")
        print("  zl.pow(3) = \(zl.pow(3))")
        print("  squareRoot = \(zl.squareRoot)")
        print("  description = \(zl.description)")
        // mixed operators with scalar
        let sl: Float80 = 2.5
        print("  zl + sl = \(zl + sl)")
        print("  sl + zl = \(sl + zl)")
        print("  zl - sl = \(zl - sl)")
        print("  sl - zl = \(sl - zl)")
        print("  zl * sl = \(zl * sl)")
        print("  sl * zl = \(sl * zl)")
        print("  zl / sl = \(zl / sl)")
        // complex division
        let wl = ComplexL(1.5, -0.75)
        print("  zl / wl = \(zl / wl)")
        // Float80 specialization: fromPolar, phase, exp/log/trig/hyperbolic/atan
        let zlp = ComplexL.fromPolar(radius: 2.0, phase: .pi / 3.0)
        print("  fromPolar(2, pi/3) = \(zlp)")
        print("  phase(zl) = \(zl.phase)")
        print("  exp(zl) = \(zl.exp)")
        print("  log(zl) = \(zl.log)")
        print("  sin(zl) = \(zl.sin)")
        print("  cos(zl) = \(zl.cos)")
        print("  tan(zl) = \(zl.tan)")
        print("  sinh(zl) = \(zl.sinh)")
        print("  cosh(zl) = \(zl.cosh)")
        print("  tanh(zl) = \(zl.tanh)")
        print("  atan(zl) = \(zl.atan)")
    }
#endif
    
    // MARK: - SpecialFunctions API coverage
    
    // Digamma / Trigamma / Polygamma / Riemann Zeta
    do {
        print("SpecialFunctions: Digamma/Trigamma/Polygamma/Zeta")
        let xd: Double = 3.5
        let xf: Float = 3.5
        print("  digamma(Double): \(try SpecialFunctions.digamma(xd))")
        print("  trigamma(Double): \(try SpecialFunctions.trigamma(xd))")
        print("  polygamma(n:2, Double): \(try SpecialFunctions.polygamma(order: 2, xd))")
        print("  zeta(Double, 2.0): \(try SpecialFunctions.riemannZeta(2.0 as Double))")
        print("  digamma(Float): \(try SpecialFunctions.digamma(xf))")
        print("  trigamma(Float): \(try SpecialFunctions.trigamma(xf))")
        print("  polygamma(n:3, Float): \(try SpecialFunctions.polygamma(order: 3, xf))")
        print("  zeta(Float, 2.0): \(try SpecialFunctions.riemannZeta(2.0 as Float))")
        // Mixed generic (T: Double) path is already covered by Double above.
#if arch(x86_64)
        let xl: Float80 = 3.5
        print("  digamma(Float80): \(try SpecialFunctions.digamma(xl))")
        print("  trigamma(Float80): \(try SpecialFunctions.trigamma(xl))")
        print("  polygamma(n:1, Float80): \(try SpecialFunctions.polygamma(order: 1, xl))")
        print("  zeta(Float80, 2.0): \(try SpecialFunctions.riemannZeta(2.0 as Float80))")
#endif
    } catch {
        print("SpecialFunctions digamma/trigamma/polygamma/zeta error: \(error)")
    }
    
    // Error functions: erf, erfc, inverse erf
    do {
        print("SpecialFunctions: Error functions (erf/erfc/inverse)")
        let xd: Double = 0.5
        let zD: Double = 0.3
        print("  erf(Double): \(try SpecialFunctions.errorFunction(xd))")
        print("  erfc(Double): \(try SpecialFunctions.complementaryErrorFunction(xd))")
        print("  erfInv(Double): \(try SpecialFunctions.inverseErrorFunction(zD))")
        
        let xf: Float = 0.5
        let zF: Float = 0.3
        print("  erf(Float): \(try SpecialFunctions.errorFunction(xf))")
        print("  erfc(Float): \(try SpecialFunctions.complementaryErrorFunction(xf))")
        print("  erfInv(Float): \(try SpecialFunctions.inverseErrorFunction(zF))")
#if arch(x86_64)
        let xl: Float80 = 0.5
        let zL: Float80 = 0.3
        print("  erf(Float80): \(try SpecialFunctions.errorFunction(xl))")
        print("  erfc(Float80): \(try SpecialFunctions.complementaryErrorFunction(xl))")
        print("  erfInv(Float80): \(try SpecialFunctions.inverseErrorFunction(zL))")
#endif
    } catch {
        print("SpecialFunctions error functions error: \(error)")
    }
    
    // Hypergeometric functions: 1F0, 0F1, 2F0, 1F1, and pFq
    do {
        print("SpecialFunctions: Hypergeometric 1F0 / 0F1 / 2F0 / 1F1 / pFq")
        // 1F0
        print("  1F0<Double>: \(try SpecialFunctions.hypergeometric1F0(a: 0.75 as Double, z: 0.2 as Double))")
        print("  1F0<Float>: \(try SpecialFunctions.hypergeometric1F0(a: 0.75 as Float, z: 0.2 as Float))")
        print("  1F0<Float,Double mixed -> Double>: \(try SpecialFunctions.hypergeometric1F0(a: 0.5 as Float, z: 0.25 as Double))")
        // 0F1
        print("  0F1<Double>: \(try SpecialFunctions.hypergeometric0F1(b: 1.25 as Double, z: 0.3 as Double))")
        print("  0F1<Float>: \(try SpecialFunctions.hypergeometric0F1(b: 1.25 as Float, z: 0.3 as Float))")
        print("  0F1<Double,Float mixed -> Double>: \(try SpecialFunctions.hypergeometric0F1(b: 1.5 as Double, z: 0.2 as Float))")
        // 2F0
        print("  2F0<Double>: \(try SpecialFunctions.hypergeometric2F0(a: 0.5 as Double, b: 1.25 as Double, z: 0.1 as Double))")
        print("  2F0<Float>: \(try SpecialFunctions.hypergeometric2F0(a: 0.5 as Float, b: 1.25 as Float, z: 0.1 as Float))")
        print("  2F0<Float,Double mixed -> Double>: \(try SpecialFunctions.hypergeometric2F0(a: 0.5 as Float, b: 1.25 as Double, z: 0.1 as Double))")
        // 1F1 (Kummer’s M)
        print("  1F1<Double>: \(try SpecialFunctions.hypergeometric1F1(a: 0.75 as Double, b: 1.5 as Double, z: 0.2 as Double))")
        print("  1F1<Float>: \(try SpecialFunctions.hypergeometric1F1(a: 0.75 as Float, b: 1.5 as Float, z: 0.2 as Float))")
        print("  1F1<Double,Float mixed -> Double>: \(try SpecialFunctions.hypergeometric1F1(a: 0.75 as Double, b: 1.5 as Float, z: 0.2 as Double))")
        // pFq arrays
        let aD: [Double] = [0.5, 1.25]
        let bD: [Double] = [1.5, 2.25, 3.0]
        let aF: [Float] = [0.5, 1.25]
        let bF: [Float] = [1.5, 2.25, 3.0]
        print("  pFq<Double>: \(SpecialFunctions.hypergeometricPFQ(a: aD, b: bD, z: 0.2))")
        print("  pFq<Float>: \(SpecialFunctions.hypergeometricPFQ(a: aF, b: bF, z: 0.2 as Float))")
        print("  pFq<[Float], Double z -> Double>: \(SpecialFunctions.hypergeometricPFQ(a: aF, b: bF, z: 0.2 as Double))")
        print("  pFq<[Double], Float z -> Double>: \(SpecialFunctions.hypergeometricPFQ(a: aD, b: bD, z: 0.2 as Float))")
#if arch(x86_64)
        let aL: [Float80] = [0.5, 1.25]
        let bL: [Float80] = [1.5, 2.25, 3.0]
        let pFqL1 = SpecialFunctions.hypergeometricPFQ(a: aL, b: bL, z: 0.2 as Float80)
        print("  pFq<Float80>: \(pFqL1)")
        let pFqL2 = SpecialFunctions.hypergeometricPFQ(a: aL, b: bL, z: 0.2 as Double)
        print("  pFq<[Float80], Double z -> Float80>: \(pFqL2)")
        let pFqL3 = SpecialFunctions.hypergeometricPFQ(a: aD, b: bD, z: 0.2 as Float80)
        print("  pFq<[Double], Float80 z -> Float80>: \(pFqL3)")
        let pFqL4 = SpecialFunctions.hypergeometricPFQ(a: aF, b: bF, z: 0.2 as Float80)
        print("  pFq<[Float], Float80 z -> Float80>: \(pFqL4)")
#endif
    } catch {
        print("SpecialFunctions hypergeometric error: \(error)")
    }
    
    // Beta and incomplete/regularized/inverse/derivative + parameter solvers
    do {
        print("SpecialFunctions: Beta family")
        let aD = 2.5 as Double, bD = 3.5 as Double, xD = 0.4 as Double, pD = 0.6 as Double
        print("  beta<Double>: \(try SpecialFunctions.beta(aD, bD))")
        print("  incompleteBetaUnnormalized<Double>: \(try SpecialFunctions.incompleteBetaUnnormalized(aD, bD, x: xD))")
        print("  regularizedIncompleteBeta<Double>: \(try SpecialFunctions.regularizedIncompleteBeta(aD, bD, x: xD))")
        print("  complementaryRegularizedIncompleteBeta<Double>: \(try SpecialFunctions.complementaryRegularizedIncompleteBeta(aD, bD, x: xD))")
        print("  inverseRegularizedIncompleteBeta<Double>: \(try SpecialFunctions.inverseRegularizedIncompleteBeta(aD, bD, p: pD))")
        print("  inverseComplementaryRegularizedIncompleteBeta<Double>: \(try SpecialFunctions.inverseComplementaryRegularizedIncompleteBeta(aD, bD, p: pD))")
        print("  solveAForRegularizedIncompleteBeta<Double>: \(SpecialFunctions.solveAForRegularizedIncompleteBeta(b: bD, x: xD, p: pD))")
        print("  solveBForRegularizedIncompleteBeta<Double>: \(SpecialFunctions.solveBForRegularizedIncompleteBeta(a: aD, x: xD, p: pD))")
        print("  regularizedIncompleteBetaDerivative<Double>: \(try SpecialFunctions.regularizedIncompleteBetaDerivative(aD, bD, x: xD))")
        // Float overloads
        let aF = 1.75 as Float, bF = 2.25 as Float, xF = 0.3 as Float, pF = 0.4 as Float
        print("  beta<Float>: \(try SpecialFunctions.beta(aF, bF))")
        print("  incompleteBetaUnnormalized<Float>: \(try SpecialFunctions.incompleteBetaUnnormalized(aF, bF, x: xF))")
        print("  regularizedIncompleteBeta<Float>: \(try SpecialFunctions.regularizedIncompleteBeta(aF, bF, x: xF))")
        print("  complementaryRegularizedIncompleteBeta<Float>: \(try SpecialFunctions.complementaryRegularizedIncompleteBeta(aF, bF, x: xF))")
        print("  inverseRegularizedIncompleteBeta<Float>: \(try SpecialFunctions.inverseRegularizedIncompleteBeta(aF, bF, p: pF))")
        print("  inverseComplementaryRegularizedIncompleteBeta<Float>: \(try SpecialFunctions.inverseComplementaryRegularizedIncompleteBeta(aF, bF, p: pF))")
        print("  solveAForRegularizedIncompleteBeta<Float>: \(SpecialFunctions.solveAForRegularizedIncompleteBeta(b: bF, x: xF, p: pF))")
        print("  solveBForRegularizedIncompleteBeta<Float>: \(SpecialFunctions.solveBForRegularizedIncompleteBeta(a: aF, x: xF, p: pF))")
        print("  regularizedIncompleteBetaDerivative<Float>: \(try SpecialFunctions.regularizedIncompleteBetaDerivative(aF, bF, x: xF))")
        // Mixed promotions examples
        print("  beta<Float,Double mixed -> Double>: \(try SpecialFunctions.beta(aF, bD))")
        print("  regularizedIncompleteBeta<Double,Float mixed -> Double>: \(try SpecialFunctions.regularizedIncompleteBeta(aD, bF, x: xD))")
        print("  inverseRegularizedIncompleteBeta<Double,Float mixed -> Double>: \(try SpecialFunctions.inverseRegularizedIncompleteBeta(aD, bF, p: pD))")
#if arch(x86_64)
        // Float80 overloads
        let aL = 2.25 as Float80, bL = 1.75 as Float80, xL = 0.35 as Float80, pL = 0.55 as Float80
        print("  beta<Float80>: \(try SpecialFunctions.beta(aL, bL))")
        print("  incompleteBetaUnnormalized<Float80>: \(try SpecialFunctions.incompleteBetaUnnormalized(aL, bL, x: xL))")
        print("  regularizedIncompleteBeta<Float80>: \(try SpecialFunctions.regularizedIncompleteBeta(aL, bL, x: xL))")
        print("  complementaryRegularizedIncompleteBeta<Float80>: \(try SpecialFunctions.complementaryRegularizedIncompleteBeta(aL, bL, x: xL))")
        print("  inverseRegularizedIncompleteBeta<Float80>: \(try SpecialFunctions.inverseRegularizedIncompleteBeta(aL, bL, p: pL))")
        print("  inverseComplementaryRegularizedIncompleteBeta<Float80>: \(try SpecialFunctions.inverseComplementaryRegularizedIncompleteBeta(aL, bL, p: pL))")
        print("  solveAForRegularizedIncompleteBeta<Float80>: \(SpecialFunctions.solveAForRegularizedIncompleteBeta(b: bL, x: xL, p: pL))")
        print("  solveBForRegularizedIncompleteBeta<Float80>: \(SpecialFunctions.solveBForRegularizedIncompleteBeta(a: aL, x: xL, p: pL))")
        print("  regularizedIncompleteBetaDerivative<Float80>: \(try SpecialFunctions.regularizedIncompleteBetaDerivative(aL, bL, x: xL))")
        // Mixed with Float80
        print("  beta<Float80,Double mixed -> Float80>: \(try SpecialFunctions.beta(aL, bD))")
        print("  regularizedIncompleteBeta<Float80, Double x mixed -> Float80>: \(try SpecialFunctions.regularizedIncompleteBeta(aL, bL, x: Double(xL)))")
#endif
    } catch {
        print("SpecialFunctions beta family error: \(error)")
    }
    
    // Gegenbauer polynomials and derivatives
    do {
        print("SpecialFunctions: Gegenbauer")
        let n = 3
        // Double
        print("  C_n^λ<Double>: \(try SpecialFunctions.gegenbauer(n: n, lambda: 0.5 as Double, x: 0.3 as Double))")
        print("  d/dx C_n^λ<Double>: \(try SpecialFunctions.gegenbauerPrime(n: n, lambda: 0.5 as Double, x: 0.3 as Double))")
        print("  d^k/dx^k C_n^λ<Double> (k=2): \(try SpecialFunctions.gegenbauerDerivative(n: n, lambda: 0.5 as Double, x: 0.3 as Double, k: 2))")
        // Float
        print("  C_n^λ<Float>: \(try SpecialFunctions.gegenbauer(n: n, lambda: 0.5 as Float, x: 0.3 as Float))")
        print("  d/dx C_n^λ<Float>: \(try SpecialFunctions.gegenbauerPrime(n: n, lambda: 0.5 as Float, x: 0.3 as Float))")
        print("  d^k/dx^k C_n^λ<Float> (k=2): \(try SpecialFunctions.gegenbauerDerivative(n: n, lambda: 0.5 as Float, x: 0.3 as Float, k: 2))")
        // Mixed promotions examples
        print("  C_n^λ<Float,Double mixed -> Double>: \(try SpecialFunctions.gegenbauer(n: n, lambda: 0.5 as Float, x: 0.3 as Double))")
        print("  d/dx C_n^λ<Double,Float mixed -> Double>: \(try SpecialFunctions.gegenbauerPrime(n: n, lambda: 0.5 as Double, x: 0.3 as Float))")
        print("  d^k/dx^k C_n^λ<Double,Float mixed -> Double>: \(try SpecialFunctions.gegenbauerDerivative(n: n, lambda: 0.5 as Double, x: 0.3 as Float, k: 2))")
#if arch(x86_64)
        // Float80
        print("  C_n^λ<Float80>: \(try SpecialFunctions.gegenbauer(n: n, lambda: 0.5 as Float80, x: 0.3 as Float80))")
        print("  d/dx C_n^λ<Float80>: \(try SpecialFunctions.gegenbauerPrime(n: n, lambda: 0.5 as Float80, x: 0.3 as Float80))")
        print("  d^k/dx^k C_n^λ<Float80> (k=2): \(try SpecialFunctions.gegenbauerDerivative(n: n, lambda: 0.5 as Float80, x: 0.3 as Float80, k: 2))")
        // Mixed with Float80
        print("  C_n^λ<Float80,Double mixed -> Float80>: \(try SpecialFunctions.gegenbauer(n: n, lambda: 0.5 as Float80, x: 0.3 as Double))")
        print("  d/dx C_n^λ<Double,Float80 mixed -> Float80>: \(try SpecialFunctions.gegenbauerPrime(n: n, lambda: 0.5 as Double, x: 0.3 as Float80))")
        print("  d^k/dx^k C_n^λ<Float80,Float mixed -> Float80>: \(try SpecialFunctions.gegenbauerDerivative(n: n, lambda: 0.5 as Float80, x: 0.3 as Float, k: 2))")
#endif
    } catch {
        print("SpecialFunctions Gegenbauer error: \(error)")
    }
    
    // Chebyshev polynomials (T, U), Clenshaw T-series, and recurrence helper
    do {
        print("SpecialFunctions: Chebyshev T/U, Clenshaw, recurrence")
        // T and U (Double)
        print("  T_n<Double>: \(try SpecialFunctions.chebyshevT(3, 0.5 as Double))")
        print("  U_n<Double>: \(try SpecialFunctions.chebyshevU(3, 0.5 as Double))")
        // T and U (Float)
        print("  T_n<Float>: \(try SpecialFunctions.chebyshevT(4, 0.3 as Float))")
        print("  U_n<Float>: \(try SpecialFunctions.chebyshevU(4, 0.3 as Float))")
        // Clenshaw T-series
        let cD: [Double] = [1.0, 0.5, -0.25, 0.125] // c0..c3
        let cF: [Float] = [1.0, 0.5, -0.25, 0.125]
        print("  Clenshaw<Double> (halfWeightC0=true): \(try SpecialFunctions.chebyshevClenshawRecurrence(cD, x: 0.4))")
        print("  Clenshaw<Double> (halfWeightC0=false): \(try SpecialFunctions.chebyshevClenshawRecurrence(cD, x: 0.4, halfWeightC0: false))")
        print("  Clenshaw<Float> (halfWeightC0=true): \(try SpecialFunctions.chebyshevClenshawRecurrence(cF, x: 0.4 as Float))")
        print("  Clenshaw<Float coeffs, Double x -> Double>: \(try SpecialFunctions.chebyshevClenshawRecurrence(cF, x: 0.4 as Double))")
        print("  Clenshaw<Double coeffs, Float x -> Double>: \(try SpecialFunctions.chebyshevClenshawRecurrence(cD, x: 0.4 as Float))")
        // Recurrence helper next-term (Double/Float)
        let xRecD = 0.5 as Double
        let T0D = 1.0 as Double, T1D = xRecD
        print("  chebyshev_next<Double> (T2): \(try SpecialFunctions.chebyshev_next(T1D, T0D, xRecD))")
        let xRecF = 0.5 as Float
        let T0F: Float = 1.0, T1F = xRecF
        print("  chebyshev_next<Float> (T2): \(try SpecialFunctions.chebyshev_next(T1F, T0F, xRecF))")
#if arch(x86_64)
        // Float80 variants
        print("  T_n<Float80>: \(try SpecialFunctions.chebyshevT(5, 0.2 as Float80))")
        print("  U_n<Float80>: \(try SpecialFunctions.chebyshevU(5, 0.2 as Float80))")
        let cL: [Float80] = [1.0, 0.25, -0.125, 0.0625]
        print("  Clenshaw<Float80> (halfWeightC0=true): \(try SpecialFunctions.chebyshevClenshawRecurrence(cL, x: 0.3 as Float80))")
        print("  Clenshaw<[Float80], Double x -> Float80>: \(try SpecialFunctions.chebyshevClenshawRecurrence(cL, x: 0.3 as Double))")
        print("  Clenshaw<[Double], Float80 x -> Float80>: \(try SpecialFunctions.chebyshevClenshawRecurrence(cD, x: 0.3 as Float80))")
        print("  Clenshaw<[Float], Float80 x -> Float80>: \(try SpecialFunctions.chebyshevClenshawRecurrence(cF, x: 0.3 as Float80))")
        let xRecL = 0.5 as Float80
        let T0L: Float80 = 1.0, T1L = xRecL
        print("  chebyshev_next<Float80> (T2): \(try SpecialFunctions.chebyshev_next(T1L, T0L, xRecL))")
#endif
    } catch {
        print("SpecialFunctions Chebyshev error: \(error)")
    }
}
