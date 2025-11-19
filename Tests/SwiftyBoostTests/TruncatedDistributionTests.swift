import Testing
@testable import SwiftyBoost

@Suite("Truncated distribution coverage")
struct TruncatedDistributionTests {
    private static let probabilities: [Double] = [0.1, 0.3, 0.5, 0.7, 0.9]
    private static let tolerance: Double = 1e-10

    private static func expectClose(_ actual: Double, _ expected: Double, context: @autoclosure () -> String) {
        let scale = max(1.0, max(abs(actual), abs(expected)))
        #expect(abs(actual - expected) <= tolerance * scale, "\(context()) mismatch - got \(actual), expected \(expected)")
    }

    private static func exercise<Base: DistributionProtocol>(
        _ name: StaticString,
        build: () throws -> Base,
        lower: Double,
        upper: Double
    ) throws where Base.RealType == Double {
        let base = try build()
        let truncated = try Distribution.TruncatedDistribution(base: base, lower: lower, upper: upper)

        let lowerMass = try base.cdf(lower)
        let upperMass = try base.cdf(upper)
        let normalization = upperMass - lowerMass
        #expect(normalization > 0, "\(name) produced zero probability mass in truncation interval")

        let width = upper - lower
        #expect(width > 0, "\(name) truncation interval collapsed")
        let interiorPoints = [0.2, 0.5, 0.8].map { lower + width * $0 }

        for x in interiorPoints {
            let truncatedPdf = try truncated.pdf(x)
            let truncatedCdf = try truncated.cdf(x)
            let truncatedSf = try truncated.sf(x)

            let basePdf = try base.pdf(x)
            let baseCdf = try base.cdf(x)

            expectClose(truncatedPdf, basePdf / normalization, context: "\(name) pdf @ \(x)")
            expectClose(truncatedCdf, (baseCdf - lowerMass) / normalization, context: "\(name) cdf @ \(x)")
            expectClose(truncatedSf, (upperMass - baseCdf) / normalization, context: "\(name) sf @ \(x)")
        }

        let guardDelta = max(width * 0.25, 0.1)
        let below = lower - guardDelta
        let above = upper + guardDelta

        #expect(try truncated.pdf(below) == 0, "\(name) pdf should vanish below the truncated support")
        #expect(try truncated.cdf(below) == 0, "\(name) cdf should be zero below the truncated support")
        #expect(try truncated.sf(below) == 1, "\(name) sf should be one below the truncated support")
        #expect(try truncated.pdf(above) == 0, "\(name) pdf should vanish above the truncated support")
        #expect(try truncated.cdf(above) == 1, "\(name) cdf should be one above the truncated support")
        #expect(try truncated.sf(above) == 0, "\(name) sf should be zero above the truncated support")

        for p in probabilities {
            let target = lowerMass + p * normalization
            let expected = try base.quantile(target)
            let actual = try truncated.quantile(p)
            expectClose(actual, expected, context: "\(name) quantile for p=\(p)")

            let tailTarget = upperMass - p * normalization
            let expectedUpper = try base.quantile(tailTarget)
            let actualUpper = try truncated.quantileComplement(p)
            expectClose(actualUpper, expectedUpper, context: "\(name) quantileComplement for q=\(p)")
        }
    }

    // Bounded-support distributions
    @Test("Arcsine base")
    func arcsine() throws {
        try Self.exercise(
            "Arcsine",
            build: { try Distribution.Arcsine<Double>(minX: 0, maxX: 2) },
            lower: 0.3,
            upper: 1.7
        )
    }

    @Test("Beta base")
    func beta() throws {
        try Self.exercise(
            "Beta",
            build: { try Distribution.Beta<Double>(alpha: 3.5, beta: 2.2) },
            lower: 0.1,
            upper: 0.8
        )
    }

    @Test("Non-central beta base")
    func nonCentralBeta() throws {
        try Self.exercise(
            "NonCentralBeta",
            build: { try Distribution.NonCentralBeta<Double>(alpha: 3.0, beta: 2.5, lambda: 1.0) },
            lower: 0.1,
            upper: 0.8
        )
    }

    @Test("Kolmogorov-Smirnov base")
    func kolmogorovSmirnov() throws {
        try Self.exercise(
            "KolmogorovSmirnov",
            build: { try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: 25) },
            lower: 0.05,
            upper: 0.5
        )
    }

    @Test("Triangular base")
    func triangular() throws {
        try Self.exercise(
            "Triangular",
            build: { try Distribution.Triangular<Double>(lower: -1.0, mode: 0.2, upper: 2.0) },
            lower: -0.5,
            upper: 1.5
        )
    }

    @Test("Uniform base")
    func uniform() throws {
        try Self.exercise(
            "Uniform",
            build: { try Distribution.Uniform<Double>(lower: -2.0, upper: 3.0) },
            lower: -0.5,
            upper: 2.0
        )
    }

    // Positive-support distributions
    @Test("Chi-squared base")
    func chiSquared() throws {
        try Self.exercise(
            "ChiSquared",
            build: { try Distribution.ChiSquared<Double>(degreesOfFreedom: 6.0) },
            lower: 0.4,
            upper: 7.5
        )
    }

    @Test("Gamma base")
    func gamma() throws {
        try Self.exercise(
            "Gamma",
            build: { try Distribution.Gamma<Double>(shape: 3.5, scale: 0.8) },
            lower: 0.5,
            upper: 5.0
        )
    }

    @Test("Exponential base")
    func exponential() throws {
        try Self.exercise(
            "Exponential",
            build: { try Distribution.Exponential<Double>(lambda: 1.4) },
            lower: 0.2,
            upper: 2.0
        )
    }

    @Test("Weibull base")
    func weibull() throws {
        try Self.exercise(
            "Weibull",
            build: { try Distribution.Weibull<Double>(shape: 1.8, scale: 2.0) },
            lower: 0.5,
            upper: 4.0
        )
    }

    @Test("Rayleigh base")
    func rayleigh() throws {
        try Self.exercise(
            "Rayleigh",
            build: { try Distribution.Rayleigh<Double>(scale: 1.1) },
            lower: 0.2,
            upper: 3.0
        )
    }

    @Test("Pareto base")
    func pareto() throws {
        try Self.exercise(
            "Pareto",
            build: { try Distribution.Pareto<Double>(scale: 1.2, shape: 3.5) },
            lower: 1.4,
            upper: 4.0
        )
    }

    @Test("Hyperexponential base")
    func hyperexponential() throws {
        try Self.exercise(
            "Hyperexponential",
            build: {
                try Distribution.Hyperexponential<Double>(
                    probabilities: [0.2, 0.3, 0.5],
                    rates: [0.5, 1.0, 1.8]
                )
            },
            lower: 0.2,
            upper: 2.5
        )
    }

    @Test("Fisher F base")
    func fisherF() throws {
        try Self.exercise(
            "FisherF",
            build: { try Distribution.FisherF<Double>(degreesOfFreedom1: 5.0, degreesOfFreedom2: 12.0) },
            lower: 0.3,
            upper: 5.0
        )
    }

    @Test("Non-central F base")
    func nonCentralF() throws {
        try Self.exercise(
            "NonCentralF",
            build: { try Distribution.NonCentralF<Double>(degreesOfFreedom1: 6.0, degreesOfFreedom2: 9.0, nonCentrality: 1.2) },
            lower: 0.2,
            upper: 4.0
        )
    }

    @Test("Chi-squared non-central base")
    func nonCentralChiSquared() throws {
        try Self.exercise(
            "NonCentralChiSquared",
            build: { try Distribution.NonCentralChiSquared<Double>(degreesOfFreedom: 5.0, nonCentrality: 2.0) },
            lower: 0.3,
            upper: 6.0
        )
    }

    @Test("Inverse chi-squared base")
    func inverseChiSquared() throws {
        try Self.exercise(
            "InverseChiSquared",
            build: { try Distribution.InverseChiSquared<Double>(degreesOfFreedom: 7.0, scale: 0.6) },
            lower: 0.2,
            upper: 3.0
        )
    }

    @Test("Inverse gamma base")
    func inverseGamma() throws {
        try Self.exercise(
            "InverseGamma",
            build: { try Distribution.InverseGamma<Double>(shape: 4.0, scale: 1.25) },
            lower: 0.2,
            upper: 3.0
        )
    }

    @Test("Inverse normal base")
    func inverseNormal() throws {
        try Self.exercise(
            "InverseNormal",
            build: { try Distribution.InverseNormal<Double>(mean: 1.2, shape: 3.5) },
            lower: 0.4,
            upper: 2.5
        )
    }

    @Test("Log-normal base")
    func logNormal() throws {
        try Self.exercise(
            "LogNormal",
            build: { try Distribution.LogNormal<Double>(location: -0.3, scale: 0.7) },
            lower: 0.2,
            upper: 2.5
        )
    }

    // Real-line distributions
    @Test("Normal base")
    func normal() throws {
        try Self.exercise(
            "Normal",
            build: { try Distribution.Normal<Double>(mean: 0.3, sd: 1.2) },
            lower: -0.8,
            upper: 1.5
        )
    }

    @Test("Student T base")
    func studentT() throws {
        try Self.exercise(
            "StudentT",
            build: { try Distribution.StudentT<Double>(degreesOfFreedom: 8.0) },
            lower: -1.2,
            upper: 1.6
        )
    }

    @Test("Non-central Student T base")
    func nonCentralStudentT() throws {
        try Self.exercise(
            "NonCentralStudentT",
            build: { try Distribution.NonCentralStudentT<Double>(degreesOfFreedom: 7.0, nonCentrality: 1.1) },
            lower: -0.5,
            upper: 2.0
        )
    }

    @Test("Logistic base")
    func logistic() throws {
        try Self.exercise(
            "Logistic",
            build: { try Distribution.Logistic<Double>(location: 0.0, scale: 1.3) },
            lower: -1.0,
            upper: 1.5
        )
    }

    @Test("Laplace base")
    func laplace() throws {
        try Self.exercise(
            "Laplace",
            build: { try Distribution.Laplace<Double>(location: 0.1, scale: 0.9) },
            lower: -0.7,
            upper: 1.5
        )
    }

    @Test("Cauchy base")
    func cauchy() throws {
        try Self.exercise(
            "Cauchy",
            build: { try Distribution.Cauchy<Double>(location: 0.2, scale: 1.1) },
            lower: -1.5,
            upper: 1.7
        )
    }

    @Test("Skew-normal base")
    func skewNormal() throws {
        try Self.exercise(
            "SkewNormal",
            build: { try Distribution.SkewNormal<Double>(location: 0.4, scale: 1.2, shape: -2.5) },
            lower: -0.6,
            upper: 2.0
        )
    }

    @Test("Extreme value (Gumbel) base")
    func extremeValue() throws {
        try Self.exercise(
            "ExtremeValueGumbel",
            build: { try Distribution.ExtremeValueGumbel<Double>(location: -0.3, scale: 1.5) },
            lower: -1.2,
            upper: 2.0
        )
    }

    @Test("Landau base")
    func landau() throws {
        try Self.exercise(
            "Landau",
            build: { try Distribution.Landau<Double>(location: -0.1, scale: 1.0) },
            lower: -0.8,
            upper: 1.8
        )
    }

    @Test("Holtsmark base")
    func holtsmark() throws {
        try Self.exercise(
            "Holtsmark",
            build: { try Distribution.Holtsmark<Double>(loc: -0.2, scale: 1.1) },
            lower: -1.0,
            upper: 2.0
        )
    }

    @Test("Map-Airy base")
    func mapAiry() throws {
        try Self.exercise(
            "MapAiry",
            build: { try Distribution.MapAiry<Double>(location: 0.2, scale: 1.0) },
            lower: -0.5,
            upper: 1.5
        )
    }

    @Test("SAS point-five base")
    func sasPointFive() throws {
        try Self.exercise(
            "SASPoint5",
            build: { try Distribution.SASPoint5<Double>(location: -0.2, scale: 1.3) },
            lower: -1.0,
            upper: 1.8
        )
    }

}
