import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Jacobi Elliptic Functions")
struct JacobiEllipticTests {

    @Test("Double wrappers match bridge")
    func doubleMatches() throws {
        let ks: [Double] = [0.0, 0.5, 0.99]
        let thetas: [Double] = [-1.0, -0.2, 0.0, 0.7, 1.3]
        for k in ks {
            for theta in thetas {
                let trio = try SpecialFunctions.jacobiElliptic(k, theta: theta)
                #expect(try SpecialFunctions.jacobiEllipticSn(k, theta: theta) == bs_jacobi_elliptic_sn_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCn(k, theta: theta) == bs_jacobi_elliptic_cn_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDn(k, theta: theta) == bs_jacobi_elliptic_dn_d(k, theta))
                #expect(trio.sn == bs_jacobi_elliptic_sn_d(k, theta))
                #expect(trio.cn == bs_jacobi_elliptic_cn_d(k, theta))
                #expect(trio.dn == bs_jacobi_elliptic_dn_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticSd(k, theta: theta) == bs_jacobi_elliptic_sd_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticSc(k, theta: theta) == bs_jacobi_elliptic_sc_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNs(k, theta: theta) == bs_jacobi_elliptic_ns_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNd(k, theta: theta) == bs_jacobi_elliptic_nd_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNc(k, theta: theta) == bs_jacobi_elliptic_nc_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDs(k, theta: theta) == bs_jacobi_elliptic_ds_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDc(k, theta: theta) == bs_jacobi_elliptic_dc_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCs(k, theta: theta) == bs_jacobi_elliptic_cs_d(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCd(k, theta: theta) == bs_jacobi_elliptic_cd_d(k, theta))
            }
        }
    }

    @Test("Float wrappers match bridge")
    func floatMatches() throws {
        let ks: [Float] = [0.1, 0.6, 0.95]
        let thetas: [Float] = [-0.9, -0.3, 0.0, 0.5, 1.0]
        for k in ks {
            for theta in thetas {
                let trio = try SpecialFunctions.jacobiElliptic(k, theta: theta)
                #expect(try SpecialFunctions.jacobiEllipticSn(k, theta: theta) == bs_jacobi_elliptic_sn_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCn(k, theta: theta) == bs_jacobi_elliptic_cn_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDn(k, theta: theta) == bs_jacobi_elliptic_dn_f(k, theta))
                #expect(trio.sn == bs_jacobi_elliptic_sn_f(k, theta))
                #expect(trio.cn == bs_jacobi_elliptic_cn_f(k, theta))
                #expect(trio.dn == bs_jacobi_elliptic_dn_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticSd(k, theta: theta) == bs_jacobi_elliptic_sd_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticSc(k, theta: theta) == bs_jacobi_elliptic_sc_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNs(k, theta: theta) == bs_jacobi_elliptic_ns_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNd(k, theta: theta) == bs_jacobi_elliptic_nd_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNc(k, theta: theta) == bs_jacobi_elliptic_nc_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDs(k, theta: theta) == bs_jacobi_elliptic_ds_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDc(k, theta: theta) == bs_jacobi_elliptic_dc_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCs(k, theta: theta) == bs_jacobi_elliptic_cs_f(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCd(k, theta: theta) == bs_jacobi_elliptic_cd_f(k, theta))
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 wrappers match bridge")
    func float80Matches() throws {
        let ks: [Float80] = [0.0, 0.7, 0.99]
        let thetas: [Float80] = [-0.8, -0.25, 0.0, 0.6, 1.1]
        for k in ks {
            for theta in thetas {
                let trio = try SpecialFunctions.jacobiElliptic(k, theta: theta)
                #expect(try SpecialFunctions.jacobiEllipticSn(k, theta: theta) == bs_jacobi_elliptic_sn_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCn(k, theta: theta) == bs_jacobi_elliptic_cn_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDn(k, theta: theta) == bs_jacobi_elliptic_dn_l(k, theta))
                #expect(trio.sn == bs_jacobi_elliptic_sn_l(k, theta))
                #expect(trio.cn == bs_jacobi_elliptic_cn_l(k, theta))
                #expect(trio.dn == bs_jacobi_elliptic_dn_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticSd(k, theta: theta) == bs_jacobi_elliptic_sd_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticSc(k, theta: theta) == bs_jacobi_elliptic_sc_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNs(k, theta: theta) == bs_jacobi_elliptic_ns_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNd(k, theta: theta) == bs_jacobi_elliptic_nd_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticNc(k, theta: theta) == bs_jacobi_elliptic_nc_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDs(k, theta: theta) == bs_jacobi_elliptic_ds_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticDc(k, theta: theta) == bs_jacobi_elliptic_dc_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCs(k, theta: theta) == bs_jacobi_elliptic_cs_l(k, theta))
                #expect(try SpecialFunctions.jacobiEllipticCd(k, theta: theta) == bs_jacobi_elliptic_cd_l(k, theta))
            }
        }
    }
    #endif

    @Test("Rejects non-finite inputs")
    func rejectsInvalid() {
        do {
            _ = try SpecialFunctions.jacobiEllipticSn(Double.nan, theta: 0.4) as Double
            #expect(Bool(false), "Expected parameterNotFinite for k")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: name, value: value) = error {
                #expect(name == "k")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        #expect(throws: SpecialFunctionError<Double>.parameterNotFinite(name: "theta", value: .infinity)) {
            _ = try SpecialFunctions.jacobiEllipticCn(0.3, theta: .infinity) as Double
        }
    }
}
