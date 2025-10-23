import Testing
@testable import SwiftyBoost

@Suite("Common helpers: rsqrt")
struct CommonHelpersRsqrtTests {

    // MARK: - Double
    @Test("rsqrt basic values (Double)")
    func rsqrtDouble() throws {
        let xs: [Double] = [1.0, 4.0, 0.25, 9.0, 0.0]
        for x in xs {
            let got = try SpecialFunctions.rsqrt(x)
            if x > 0 {
                let exp = 1.0 / x.squareRoot()
                #expect(abs(got - exp) <= 1e-15, "Mismatch for x=\(x)")
            } else {
                // x == 0 → 1 / sqrt(0) = +∞
                #expect(got.isInfinite && got > 0, "Expected +inf for x=0")
            }
        }
    }

    @Test("rsqrt negative domain throws (Double)")
    func rsqrtDoubleDomain() {
        #expect(throws: SpecialFunctionError<Double>.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity)) {
            _ = try SpecialFunctions.rsqrt(-1.0)
        }
    }

    // MARK: - Float
    @Test("rsqrt basic values (Float)")
    func rsqrtFloat() throws {
        let xs: [Float] = [1.0, 4.0, 0.25, 9.0, 0.0]
        for x in xs {
            let got = try SpecialFunctions.rsqrt(x)
            if x > 0 {
                let exp: Float = 1.0 / x.squareRoot()
                #expect(abs(got - exp) <= 1e-6, "Mismatch for x=\(x)")
            } else {
                #expect(got.isInfinite && got > 0, "Expected +inf for x=0")
            }
        }
    }

    @Test("rsqrt negative domain throws (Float)")
    func rsqrtFloatDomain() {
        #expect(throws: SpecialFunctionError<Float>.parameterOutOfRange(name: "x", min: 0.0, max: Float.infinity)) {
            _ = try SpecialFunctions.rsqrt(Float(-1.0))
        }
    }

    // MARK: - Float80 (x86_64)
    #if arch(x86_64)
    @Test("rsqrt basic values (Float80)")
    func rsqrtFloat80() throws {
        let xs: [Float80] = [1.0, 4.0, 0.25, 9.0, 0.0]
        for x in xs {
            let got = try SpecialFunctions.rsqrt(x)
            if x > 0 {
                let exp: Float80 = 1.0 / x.squareRoot()
                #expect(abs(got - exp) <= 1e-16, "Mismatch for x=\(x)")
            } else {
                #expect(got.isInfinite && got > 0, "Expected +inf for x=0")
            }
        }
    }

    @Test("rsqrt negative domain throws (Float80)")
    func rsqrtFloat80Domain() {
        #expect(throws: SpecialFunctionError<Float80>.parameterOutOfRange(name: "x", min: 0.0, max: Float80.infinity)) {
            _ = try SpecialFunctions.rsqrt(Float80(-1.0))
        }
    }
    #endif
}
