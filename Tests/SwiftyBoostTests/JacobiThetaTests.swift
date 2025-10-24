import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Jacobi Theta Functions")
struct JacobiThetaTests {

    @Test("Double nome overloads match bridge")
    func doubleNomeMatchesBridge() throws {
        let xs: [Double] = [-1.25, -0.5, 0.0, 0.75]
        let qs: [Double] = [0.1, 0.25, 0.9]
        for x in xs {
            for q in qs {
                let theta1 = try SpecialFunctions.jacobiTheta1(x, q: q)
                let theta2 = try SpecialFunctions.jacobiTheta2(x, q: q)
                let theta3 = try SpecialFunctions.jacobiTheta3(x, q: q)
                let theta4 = try SpecialFunctions.jacobiTheta4(x, q: q)
                let expected1 = bs_jacobi_theta1_d(x, q)
                let expected2 = bs_jacobi_theta2_d(x, q)
                let expected3 = bs_jacobi_theta3_d(x, q)
                let expected4 = bs_jacobi_theta4_d(x, q)
                if theta1.isNaN || expected1.isNaN {
                    #expect(theta1.isNaN && expected1.isNaN, "θ1 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta1 == expected1, "θ1 mismatch at x=\(x), q=\(q)")
                }
                if theta2.isNaN || expected2.isNaN {
                    #expect(theta2.isNaN && expected2.isNaN, "θ2 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta2 == expected2, "θ2 mismatch at x=\(x), q=\(q)")
                }
                if theta3.isNaN || expected3.isNaN {
                    #expect(theta3.isNaN && expected3.isNaN, "θ3 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta3 == expected3, "θ3 mismatch at x=\(x), q=\(q)")
                }
                if theta4.isNaN || expected4.isNaN {
                    #expect(theta4.isNaN && expected4.isNaN, "θ4 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta4 == expected4, "θ4 mismatch at x=\(x), q=\(q)")
                }
            }
        }
    }

    @Test("Double τ overloads match bridge")
    func doubleTauMatchesBridge() throws {
        let xs: [Double] = [-1.25, -0.5, 0.0, 0.75]
        let taus: [Double] = [-0.8, 0.0, 0.5, 1.2]
        for x in xs {
            for tau in taus {
                let theta1 = try SpecialFunctions.jacobiTheta1Tau(x, tau: tau)
                let theta2 = try SpecialFunctions.jacobiTheta2Tau(x, tau: tau)
                let theta3 = try SpecialFunctions.jacobiTheta3Tau(x, tau: tau)
                let theta4 = try SpecialFunctions.jacobiTheta4Tau(x, tau: tau)
                let expected1 = bs_jacobi_theta1tau_d(x, tau)
                let expected2 = bs_jacobi_theta2tau_d(x, tau)
                let expected3 = bs_jacobi_theta3tau_d(x, tau)
                let expected4 = bs_jacobi_theta4tau_d(x, tau)
                if theta1.isNaN || expected1.isNaN {
                    #expect(theta1.isNaN && expected1.isNaN, "θ1τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta1 == expected1, "θ1τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta2.isNaN || expected2.isNaN {
                    #expect(theta2.isNaN && expected2.isNaN, "θ2τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta2 == expected2, "θ2τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta3.isNaN || expected3.isNaN {
                    #expect(theta3.isNaN && expected3.isNaN, "θ3τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta3 == expected3, "θ3τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta4.isNaN || expected4.isNaN {
                    #expect(theta4.isNaN && expected4.isNaN, "θ4τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta4 == expected4, "θ4τ mismatch at x=\(x), tau=\(tau)")
                }
            }
        }
    }

    @Test("Float nome overloads match bridge")
    func floatNomeMatchesBridge() throws {
        let xs: [Float] = [-1.25, -0.5, 0.0, 0.75]
        let qs: [Float] = [0.1, 0.25, 0.9]
        for x in xs {
            for q in qs {
                let theta1 = try SpecialFunctions.jacobiTheta1(x, q: q)
                let theta2 = try SpecialFunctions.jacobiTheta2(x, q: q)
                let theta3 = try SpecialFunctions.jacobiTheta3(x, q: q)
                let theta4 = try SpecialFunctions.jacobiTheta4(x, q: q)
                let expected1 = bs_jacobi_theta1_f(x, q)
                let expected2 = bs_jacobi_theta2_f(x, q)
                let expected3 = bs_jacobi_theta3_f(x, q)
                let expected4 = bs_jacobi_theta4_f(x, q)
                if theta1.isNaN || expected1.isNaN {
                    #expect(theta1.isNaN && expected1.isNaN, "θ1 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta1 == expected1, "θ1 mismatch at x=\(x), q=\(q)")
                }
                if theta2.isNaN || expected2.isNaN {
                    #expect(theta2.isNaN && expected2.isNaN, "θ2 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta2 == expected2, "θ2 mismatch at x=\(x), q=\(q)")
                }
                if theta3.isNaN || expected3.isNaN {
                    #expect(theta3.isNaN && expected3.isNaN, "θ3 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta3 == expected3, "θ3 mismatch at x=\(x), q=\(q)")
                }
                if theta4.isNaN || expected4.isNaN {
                    #expect(theta4.isNaN && expected4.isNaN, "θ4 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta4 == expected4, "θ4 mismatch at x=\(x), q=\(q)")
                }
            }
        }
    }

    @Test("Float τ overloads match bridge")
    func floatTauMatchesBridge() throws {
        let xs: [Float] = [-1.25, -0.5, 0.0, 0.75]
        let taus: [Float] = [-0.8, 0.0, 0.5, 1.2]
        for x in xs {
            for tau in taus {
                let theta1 = try SpecialFunctions.jacobiTheta1Tau(x, tau: tau)
                let theta2 = try SpecialFunctions.jacobiTheta2Tau(x, tau: tau)
                let theta3 = try SpecialFunctions.jacobiTheta3Tau(x, tau: tau)
                let theta4 = try SpecialFunctions.jacobiTheta4Tau(x, tau: tau)
                let expected1 = bs_jacobi_theta1tau_f(x, tau)
                let expected2 = bs_jacobi_theta2tau_f(x, tau)
                let expected3 = bs_jacobi_theta3tau_f(x, tau)
                let expected4 = bs_jacobi_theta4tau_f(x, tau)
                if theta1.isNaN || expected1.isNaN {
                    #expect(theta1.isNaN && expected1.isNaN, "θ1τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta1 == expected1, "θ1τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta2.isNaN || expected2.isNaN {
                    #expect(theta2.isNaN && expected2.isNaN, "θ2τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta2 == expected2, "θ2τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta3.isNaN || expected3.isNaN {
                    #expect(theta3.isNaN && expected3.isNaN, "θ3τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta3 == expected3, "θ3τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta4.isNaN || expected4.isNaN {
                    #expect(theta4.isNaN && expected4.isNaN, "θ4τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta4 == expected4, "θ4τ mismatch at x=\(x), tau=\(tau)")
                }
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 nome overloads match bridge")
    func float80NomeMatchesBridge() throws {
        let xs: [Float80] = [-1.25, -0.5, 0.0, 0.75]
        let qs: [Float80] = [0.1, 0.25, 0.9]
        for x in xs {
            for q in qs {
                let theta1 = try SpecialFunctions.jacobiTheta1(x, q: q)
                let theta2 = try SpecialFunctions.jacobiTheta2(x, q: q)
                let theta3 = try SpecialFunctions.jacobiTheta3(x, q: q)
                let theta4 = try SpecialFunctions.jacobiTheta4(x, q: q)
                let expected1 = bs_jacobi_theta1_l(x, q)
                let expected2 = bs_jacobi_theta2_l(x, q)
                let expected3 = bs_jacobi_theta3_l(x, q)
                let expected4 = bs_jacobi_theta4_l(x, q)
                if theta1.isNaN || expected1.isNaN {
                    #expect(theta1.isNaN && expected1.isNaN, "θ1 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta1 == expected1, "θ1 mismatch at x=\(x), q=\(q)")
                }
                if theta2.isNaN || expected2.isNaN {
                    #expect(theta2.isNaN && expected2.isNaN, "θ2 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta2 == expected2, "θ2 mismatch at x=\(x), q=\(q)")
                }
                if theta3.isNaN || expected3.isNaN {
                    #expect(theta3.isNaN && expected3.isNaN, "θ3 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta3 == expected3, "θ3 mismatch at x=\(x), q=\(q)")
                }
                if theta4.isNaN || expected4.isNaN {
                    #expect(theta4.isNaN && expected4.isNaN, "θ4 NaN mismatch at x=\(x), q=\(q)")
                } else {
                    #expect(theta4 == expected4, "θ4 mismatch at x=\(x), q=\(q)")
                }
            }
        }
    }

    @Test("Float80 τ overloads match bridge")
    func float80TauMatchesBridge() throws {
        let xs: [Float80] = [-1.25, -0.5, 0.0, 0.75]
        let taus: [Float80] = [-0.8, 0.0, 0.5, 1.2]
        for x in xs {
            for tau in taus {
                let theta1 = try SpecialFunctions.jacobiTheta1Tau(x, tau: tau)
                let theta2 = try SpecialFunctions.jacobiTheta2Tau(x, tau: tau)
                let theta3 = try SpecialFunctions.jacobiTheta3Tau(x, tau: tau)
                let theta4 = try SpecialFunctions.jacobiTheta4Tau(x, tau: tau)
                let expected1 = bs_jacobi_theta1tau_l(x, tau)
                let expected2 = bs_jacobi_theta2tau_l(x, tau)
                let expected3 = bs_jacobi_theta3tau_l(x, tau)
                let expected4 = bs_jacobi_theta4tau_l(x, tau)
                if theta1.isNaN || expected1.isNaN {
                    #expect(theta1.isNaN && expected1.isNaN, "θ1τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta1 == expected1, "θ1τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta2.isNaN || expected2.isNaN {
                    #expect(theta2.isNaN && expected2.isNaN, "θ2τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta2 == expected2, "θ2τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta3.isNaN || expected3.isNaN {
                    #expect(theta3.isNaN && expected3.isNaN, "θ3τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta3 == expected3, "θ3τ mismatch at x=\(x), tau=\(tau)")
                }
                if theta4.isNaN || expected4.isNaN {
                    #expect(theta4.isNaN && expected4.isNaN, "θ4τ NaN mismatch at x=\(x), tau=\(tau)")
                } else {
                    #expect(theta4 == expected4, "θ4τ mismatch at x=\(x), tau=\(tau)")
                }
            }
        }
    }
    #endif

    @Test("Rejects invalid nome q")
    func rejectsInvalidNome() {
        let invalidQs: [Double] = [-0.5, 0.0, 1.0, 1.5, .infinity, -.infinity]
        for q in invalidQs {
            #expect(throws: SpecialFunctionError<Double>.parameterNotInDomain(name: "q", value: q)) {
                _ = try SpecialFunctions.jacobiTheta1(0.3, q: q)
            }
            #expect(throws: SpecialFunctionError<Double>.parameterNotInDomain(name: "q", value: q)) {
                _ = try SpecialFunctions.jacobiTheta2(0.3, q: q)
            }
            #expect(throws: SpecialFunctionError<Double>.parameterNotInDomain(name: "q", value: q)) {
                _ = try SpecialFunctions.jacobiTheta3(0.3, q: q)
            }
            #expect(throws: SpecialFunctionError<Double>.parameterNotInDomain(name: "q", value: q)) {
                _ = try SpecialFunctions.jacobiTheta4(0.3, q: q)
            }
        }
        do {
            _ = try SpecialFunctions.jacobiTheta2(0.3, q: .nan)
            #expect(Bool(false), "Expected parameterNotInDomain for q=NaN")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotInDomain(name: name, value: value) = error {
                #expect(name == "q")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }

    @Test("Rejects non-finite x or τ")
    func rejectsNonFiniteInputs() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotFinite(name: "x", value: .infinity)) {
            _ = try SpecialFunctions.jacobiTheta1(.infinity, q: 0.5)
        }
        do {
            _ = try SpecialFunctions.jacobiTheta1Tau(0.2, tau: .nan)
            #expect(Bool(false), "Expected parameterNotFinite for τ")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: name, value: value) = error {
                #expect(name == "tau")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }
}
