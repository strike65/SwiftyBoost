//
//  Created by VT on 15.10.25.
//

import Testing
@testable import SwiftyBoost

@Suite("Quadrature integration (Double)")
struct QuadratureDoubleIntegrationTests {

    @Test("Gauss-Legendre integrates affine functions on finite intervals")
    func gaussLegendreFiniteInterval() throws {
        let integrator = try Quadrature.Integrator<Double>(rule: .gaussLegendre(points: 30))

        let result = try integrator.integrate(over: .finite(lower: 0.0, upper: 1.0)) { x in x + 1.0 }

        #expect(abs(result.value - 1.5) < 1e-12)
        #expect(result.didConverge)
        #expect(result.metadata.rule == .gaussLegendre(points: 30))
        #expect(result.metadata.points == 30)
        #expect(!result.metadata.isAdaptive)

        var abscissa = Array(repeating: 0.0, count: result.metadata.points)
        var weights = Array(repeating: 0.0, count: result.metadata.points)
        let copied = abscissa.withUnsafeMutableBufferPointer { abPtr in
            weights.withUnsafeMutableBufferPointer { wtPtr in
                integrator.copyAbscissaWeights(into: abPtr, weights: wtPtr)
            }
        }
        #expect(copied)
        #expect(abscissa.allSatisfy { $0.isFinite && abs($0) <= 1.0 })
        #expect(weights.allSatisfy { $0.isFinite && $0 > 0.0 })
    }

    @Test("Convenience integrate for sin(x) on [0, π] via Gauss-Kronrod")
    func gaussKronrodSine() throws {
        let result = try Quadrature.integrate(
            using: Quadrature.Rule.gaussKronrod(points: 21),
            over: .finite(lower: 0.0, upper: Double.pi)
        ) { (x: Double) in Double.sin(x) }

        #expect(abs(result.value - 2.0) < 1e-12)
        #expect(result.didConverge)
        #expect(result.metadata.rule == Quadrature.Rule.gaussKronrod(points: 21))
        #expect(!result.metadata.isAdaptive)
        #expect(!result.metadata.supportsInfiniteBounds)
    }

    @Test("Tanh-Sinh handles infinite bounds and declines abscissa copy")
    func tanhSinhAutomaticInterval() throws {
        let integrator = try Quadrature.Integrator<Double>(rule: .tanhSinh())

        let result = try integrator.integrate { x in 1.0 / (1.0 + x * x) }

        #expect(abs(result.value - Double.pi / 2.0) < 1e-9)
        #expect(result.didConverge)
        #expect(result.metadata.isAdaptive)
        #expect(!result.metadata.supportsInfiniteBounds)
        #expect(result.metadata.points == -1)

        var abscissa = Array(repeating: 0.0, count: 4)
        var weights = Array(repeating: 0.0, count: 4)
        let copied = abscissa.withUnsafeMutableBufferPointer { abPtr in
            weights.withUnsafeMutableBufferPointer { wtPtr in
                integrator.copyAbscissaWeights(into: abPtr, weights: wtPtr)
            }
        }
        #expect(!copied)
    }

    @Test("Invalid interval rejects unordered bounds")
    func invalidInterval() throws {
        let integrator = try Quadrature.Integrator<Double>(rule: .gaussLegendre(points: 10))

        #expect(throws: Quadrature.Error.invalidInterval(lower: 1.0, upper: 1.0)) {
            _ = try integrator.integrate(over: .finite(lower: 1.0, upper: 1.0)) { $0 }
        }
    }
}
