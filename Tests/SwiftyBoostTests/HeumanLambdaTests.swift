import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Heuman Lambda")
struct HeumanLambdaTests {

    @Test("Double matches bridge")
    func doubleMatches() throws {
        let ks: [Double] = [0.0, 0.3, 0.9]
        let phis: [Double] = [-1.0, -0.25, 0.0, 0.8, 1.2]
        for k in ks {
            for phi in phis {
                let expected = bs_heuman_lambda_d(k, phi)
                let got = try SpecialFunctions.heumanLambda(k, phi: phi)
                #expect(got == expected)
            }
        }
    }

    @Test("Float matches bridge")
    func floatMatches() throws {
        let ks: [Float] = [0.1, 0.5, 0.95]
        let phis: [Float] = [-0.7, -0.2, 0.0, 0.6, 1.0]
        for k in ks {
            for phi in phis {
                let expected = bs_heuman_lambda_f(k, phi)
                let got = try SpecialFunctions.heumanLambda(k, phi: phi)
                #expect(got == expected)
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 matches bridge")
    func float80Matches() throws {
        let ks: [Float80] = [0.0, 0.6, 0.99]
        let phis: [Float80] = [-0.8, -0.3, 0.0, 0.75, 1.3]
        for k in ks {
            for phi in phis {
                let expected = bs_heuman_lambda_l(k, phi)
                let got = try SpecialFunctions.heumanLambda(k, phi: phi)
                #expect(got == expected)
            }
        }
    }
    #endif

    @Test("Rejects non-finite input")
    func rejectsInvalid() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotFinite(name: "k", value: .infinity)) {
            _ = try SpecialFunctions.heumanLambda(.infinity, phi: 0.2) as Double
        }
        do {
            _ = try SpecialFunctions.heumanLambda(0.3, phi: Double.nan)
            #expect(Bool(false), "Expected parameterNotFinite for phi")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: name, value: value) = error {
                #expect(name == "phi")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }
}
