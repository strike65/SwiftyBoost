import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Airy Functions")
struct AiryTests {

    @Test("Double Ai/Bi match bridge")
    func doubleMatchesBridge() throws {
        let xs: [Double] = [-10, -2.5, -0.5, 0.0, 0.25, 1.5, 5.0]
        for x in xs {
            #expect(try SpecialFunctions.airyAi(x) == bs_airy_ai_d(x))
            #expect(try SpecialFunctions.airyBi(x) == bs_airy_bi_d(x))
            #expect(try SpecialFunctions.airyAiPrime(x) == bs_airy_ai_prime_d(x))
            #expect(try SpecialFunctions.airyBiPrime(x) == bs_airy_bi_prime_d(x))
        }
    }

    @Test("Float Ai/Bi match bridge")
    func floatMatchesBridge() throws {
        let xs: [Float] = [-5.0, -0.75, 0.0, 0.3, 1.0, 3.5]
        for x in xs {
            #expect(try SpecialFunctions.airyAi(x) == bs_airy_ai_f(x))
            #expect(try SpecialFunctions.airyBi(x) == bs_airy_bi_f(x))
            #expect(try SpecialFunctions.airyAiPrime(x) == bs_airy_ai_prime_f(x))
            #expect(try SpecialFunctions.airyBiPrime(x) == bs_airy_bi_prime_f(x))
        }
    }

    @Test("Zero finders (Double) match bridge")
    func doubleZerosMatchBridge() throws {
        for n in 0..<6 {
            #expect(try SpecialFunctions.airyAiZero(n) == bs_airy_ai_zero_d(Int32(n)))
            #expect(try SpecialFunctions.airyBiZero(n) == bs_airy_bi_zero_d(Int32(n)))
        }
        let aiZeros = try SpecialFunctions.airyAiZeros(startIndex: 0, count: 5) as [Double]
        var expectedAi = Array<Double>(repeating: 0, count: 5)
        expectedAi.withUnsafeMutableBufferPointer { buf in
            bs_airy_ai_zeros_d(0, 5, buf.baseAddress)
        }
        #expect(aiZeros == expectedAi)

        let biZeros = try SpecialFunctions.airyBiZeros(startIndex: 2, count: 3) as [Double]
        var expectedBi = Array<Double>(repeating: 0, count: 3)
        expectedBi.withUnsafeMutableBufferPointer { buf in
            bs_airy_bi_zeros_d(2, 3, buf.baseAddress)
        }
        #expect(biZeros == expectedBi)
    }

    @Test("Zero finders (Float) match bridge")
    func floatZerosMatchBridge() throws {
        for n in 0..<5 {
            #expect(try SpecialFunctions.airyAiZero(n) == bs_airy_ai_zero_f(Int32(n)))
            #expect(try SpecialFunctions.airyBiZero(n) == bs_airy_bi_zero_f(Int32(n)))
        }
        let aiZeros = try SpecialFunctions.airyAiZeros(startIndex: 1, count: 4)
        var expectedAi = Array<Float>(repeating: 0, count: 4)
        expectedAi.withUnsafeMutableBufferPointer { buf in
            bs_airy_ai_zeros_f(1, 4, buf.baseAddress)
        }
        #expect(aiZeros == expectedAi)

        let biZeros = try SpecialFunctions.airyBiZeros(startIndex: 0, count: 3)
        var expectedBi = Array<Float>(repeating: 0, count: 3)
        expectedBi.withUnsafeMutableBufferPointer { buf in
            bs_airy_bi_zeros_f(0, 3, buf.baseAddress)
        }
        #expect(biZeros == expectedBi)
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 Ai/Bi match bridge")
    func float80MatchesBridge() throws {
        let xs: [Float80] = [-6.0, -1.1, 0.0, 0.45, 2.0]
        for x in xs {
            #expect(try SpecialFunctions.airyAi(x) == bs_airy_ai_l(x))
            #expect(try SpecialFunctions.airyBi(x) == bs_airy_bi_l(x))
            #expect(try SpecialFunctions.airyAiPrime(x) == bs_airy_ai_prime_l(x))
            #expect(try SpecialFunctions.airyBiPrime(x) == bs_airy_bi_prime_l(x))
        }
    }

    @Test("Zero finders (Float80) match bridge")
    func float80ZerosMatchBridge() throws {
        for n in 0..<4 {
            #expect(try SpecialFunctions.airyAiZero(n) == bs_airy_ai_zero_l(Int32(n)))
            #expect(try SpecialFunctions.airyBiZero(n) == bs_airy_bi_zero_l(Int32(n)))
        }
        let aiZeros = try SpecialFunctions.airyAiZeros(startIndex: 0, count: 3) as [Float80]
        var expectedAi = Array<Float80>(repeating: 0, count: 3)
        expectedAi.withUnsafeMutableBufferPointer { buf in
            bs_airy_ai_zeros_l(0, 3, buf.baseAddress)
        }
        #expect(aiZeros == expectedAi)
    }
    #endif

    @Test("Rejects non-finite x")
    func rejectsNonFinite() {
        do {
            _ = try SpecialFunctions.airyAi(.infinity)
            #expect(Bool(false), "Expected parameterNotFinite for x=+∞")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: name, value: value) = error {
                #expect(name == "x")
                #expect(value.isInfinite)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch let error as SpecialFunctionError<Float> {
            if case let .parameterNotFinite(name: name, value: value) = error {
                #expect(name == "x")
                #expect(value.isInfinite)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }

        do {
            _ = try SpecialFunctions.airyBi(Double.nan)
            #expect(Bool(false), "Expected parameterNotFinite for x=NaN")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: name, value: value) = error {
                #expect(name == "x")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }

    @Test("Zero helpers reject invalid input")
    func zerosRejectInvalid() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "n", value: -1)) {
            _ = try SpecialFunctions.airyAiZero(-1) as Double
        }
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "count", value: -2)) {
            _ = try SpecialFunctions.airyAiZeros(startIndex: 0, count: -2) as [Double]
        }
    }
}
